"""Experimental two-level MX visibility blocks; NOT a complete MS archive format.

See devnotes/microscaling.md. Standard finite element encodings and UE8M0
micro-scales are combined with a custom float32 outer scale, adaptive format
tags, exact exceptions, and a bounded independently decodable record.
"""

from dataclasses import dataclass
import struct
import zlib

import numpy as np

from .codec import pack_codes, unpack_codes
from .visibility import protected_mask

HEADER = struct.Struct("<4sBBBBHHfIIII")
MAX_RAW = 16384
# ID: (name, total bits, exponent bits, fraction bits, exponent bias, max code)
FORMATS = {
    1: ("E2M1", 4, 2, 1, 1, 7),
    2: ("E2M3", 6, 2, 3, 1, 31),
    3: ("E3M2", 6, 3, 2, 3, 31),
    4: ("E4M3", 8, 4, 3, 7, 126),
    5: ("E5M2", 8, 5, 2, 15, 123),
    6: ("custom-scaled-FP16", 16, 5, 10, 15, 31743),
}


def positive_levels(fid):
    _, _, _, mantissa, bias, maximum = FORMATS[fid]
    code = np.arange(maximum + 1, dtype=np.uint16)
    exponent, fraction = code >> mantissa, code & ((1 << mantissa) - 1)
    return np.ldexp(np.where(exponent == 0, 0., 1.) + fraction / (1 << mantissa),
                    np.maximum(exponent.astype(int), 1) - bias).astype(np.float32)


LEVELS = {fid: positive_levels(fid) for fid in FORMATS}
MIDPOINTS = {fid: (v[:-1].astype(np.float64) + v[1:]) / 2 for fid, v in LEVELS.items()}


def quantize_elements(values, fid):
    """Finite inputs: nearest, ties to even low code bit, including subnormals."""
    if fid not in FORMATS or not np.all(np.isfinite(values)):
        raise ValueError("Invalid format or nonfinite unprotected element")
    magnitude = np.abs(values)
    if np.any(magnitude > LEVELS[fid][-1]):
        raise ValueError("Element outside supported range; use exceptions")
    if fid == 6:
        return values.astype(np.float16).view(np.uint16)
    mids = MIDPOINTS[fid]
    code = np.searchsorted(mids, magnitude, side="left").astype(np.uint16)
    tie = (code < len(mids)) & (magnitude == mids[np.minimum(code, len(mids)-1)])
    code += (tie & ((code & 1) != 0)).astype(np.uint16)
    return code | (np.signbit(values).astype(np.uint16) << (FORMATS[fid][1] - 1))


def element_values(codes, fid):
    bits = FORMATS[fid][1]
    code = np.asarray(codes, dtype=np.uint16)
    magnitude = code & ((1 << (bits - 1)) - 1)
    if np.any(code >= 2**bits) or np.any(magnitude >= len(LEVELS[fid])):
        raise ValueError("Reserved/nonfinite element code")
    return np.copysign(LEVELS[fid][magnitude], np.where(code >> (bits-1), -1., 1.)).astype(np.float32)


def _micro(values):
    # [block, baseline, channel] complex -> [block, channel*2+component, baseline].
    b, r, c = values.shape
    return np.ascontiguousarray(values).view(np.float32).reshape(b, r, c, 2).transpose(0, 2, 3, 1).reshape(b, 2*c, r)


def _complex(values):
    b, m, r = values.shape
    return np.ascontiguousarray(values.reshape(b, m//2, 2, r).transpose(0, 3, 1, 2)).view(np.complex64).reshape(b, r, m//2)


def _reconstruct(codes, exponents, fid, outer):
    # Decoder order: float64 ldexp, float64 outer multiply, one float32 rounding.
    with np.errstate(over="ignore", invalid="ignore"):
        return (np.ldexp(element_values(codes, fid).astype(np.float64), exponents[..., None])
                * np.asarray(outer, dtype=np.float64)[..., None, None]).astype(np.float32)


def _candidate(source, protected, outer, fid):
    x = _micro(np.where(protected, np.complex64(0), source)).astype(np.float64)
    maxima = np.max(np.abs(x), axis=-1)
    ratio = maxima / outer[:, None].astype(np.float64) / float(LEVELS[fid][-1])
    with np.errstate(divide="ignore"):
        exponent = np.where(ratio > 0, np.ceil(np.log2(ratio)), 0).clip(-127, 127).astype(np.int16)
    denominator = np.ldexp(outer[:, None].astype(np.float64), exponent)
    normalized = x / denominator[..., None]
    overflow = np.abs(normalized) > LEVELS[fid][-1]
    # Floating log2 can round down at an exact power boundary: raise that scale.
    raise_scale = overflow.any(axis=-1) & (exponent < 127)
    exponent += raise_scale.astype(np.int16)
    denominator = np.ldexp(outer[:, None].astype(np.float64), exponent)
    normalized = x / denominator[..., None]
    overflow = np.abs(normalized) > LEVELS[fid][-1]
    codes = quantize_elements(np.where(overflow, 0, normalized), fid)
    decoded = _complex(_reconstruct(codes, exponent, fid, outer))
    overflow_complex = overflow.reshape(len(source), source.shape[2], 2, 32).any(axis=2).transpose(0, 2, 1)
    exceptions = protected | overflow_complex | ~np.isfinite(decoded)
    decoded[exceptions] = source[exceptions]
    return codes, (exponent + 127).astype(np.uint8), decoded, exceptions


def _passing(source, decoded, valid, bounds):
    error = np.zeros(source.shape, dtype=np.complex128)
    np.subtract(decoded, source, out=error, where=valid, dtype=np.complex128)
    square = error.real**2 + error.imag**2
    count = np.maximum(1, valid.sum(axis=(1, 2)))
    metrics = {
        "max_component_mse": square.sum(axis=(1, 2))/(2*count),
        "max_baseline_component_mse": (square.sum(axis=2)/(2*np.maximum(1, valid.sum(axis=2)))).max(axis=1),
        "max_complex_error": np.sqrt(square.max(axis=(1, 2))),
        "max_complex_bias": np.abs(error.sum(axis=(1, 2)) / count),
    }
    passed = np.all(np.isfinite(decoded) | ~valid, axis=(1, 2))
    for key, value in bounds.items():
        if value is not None:
            if key not in metrics:
                raise ValueError("Experimental MX does not implement " + key)
            passed &= metrics[key] <= value
    return passed


def _frame(raw, rows, channels, outer, fid, exceptions, compress):
    zipped = zlib.compress(raw, 1) if compress else raw
    use_zip = compress and len(zipped) < len(raw)
    payload = zipped if use_zip else raw
    head = HEADER.pack(b"MXB1", int(use_zip), fid, 0, 0, rows, channels, outer,
                       exceptions, len(raw), len(payload), 0)
    checksum = zlib.crc32(payload, zlib.crc32(head[:-4]))
    return head[:-4] + struct.pack("<I", checksum) + payload


def _serialize(source, codes, scales, tags, outer, exceptions, rows, compress):
    channels = source.shape[1]
    unique = np.unique(tags)
    fid = int(unique[0]) if len(unique) == 1 and unique[0] != 0 else 255
    chunks = [tags.astype(np.uint8).tobytes()] if fid == 255 else []
    chunks.append(scales.astype(np.uint8).tobytes())
    components = _micro(source[None])[0]
    code_bytes = 0
    # Format banks appear by increasing ID; microblocks retain channel/component order.
    for tag in sorted(unique):
        selected = tags == tag
        body = (components[selected].astype("<f4").tobytes() if tag == 0 else
                pack_codes(codes[selected], FORMATS[int(tag)][1]))
        code_bytes += len(body)
        chunks.append(body)
    indices = np.flatnonzero(exceptions[:rows].ravel()).astype("<u2")
    chunks.extend([indices.tobytes(), source[:rows].ravel()[indices].astype("<c8").tobytes()])
    raw = b"".join(chunks)
    blob = _frame(raw, rows, channels, outer, fid, len(indices), compress)
    return blob, {"header": HEADER.size, "format_tags": len(tags) if fid == 255 else 0,
                  "micro_scales": len(scales), "codes": code_bytes,
                  "exceptions": len(indices)*10, "raw_body": len(raw),
                  "stored_body": len(blob)-HEADER.size}


@dataclass
class Encoded:
    blob: bytes
    accounting: dict
    formats: dict


def encode_batch(values, flags, rows, bounds, *, outer_mode="max", granularity="micro", compress=True):
    """Encode padded [n,32,channels<=24] arrays. Rows records actual presence.

    Uniform-format candidates compete with exact storage. Micro mode additionally
    tests a conservative local-error allocation, then checks complete megablocks.
    This heuristic is not a globally optimal precision allocator.
    """
    values = np.asarray(values)
    flags = np.asarray(flags)
    rows = np.asarray(rows)
    if (values.dtype != np.complex64 or values.ndim != 3 or values.shape[1] != 32
            or not 1 <= values.shape[2] <= 24 or flags.dtype != bool or flags.shape != values.shape
            or rows.shape != (len(values),) or np.any(rows < 1) or np.any(rows > 32)
            or outer_mode not in ("unit", "max") or granularity not in ("mega", "micro")):
        raise ValueError("Unsupported experimental MX input/configuration")
    present = np.arange(32)[None, :, None] < rows[:, None, None]
    source = np.where(present, values, np.complex64(0))
    protected = protected_mask(source, flags) & present
    valid = ~flags & np.isfinite(source) & present
    safe = np.where(protected, np.complex64(0), source)
    outer = np.max(np.abs(_micro(safe)), axis=(1, 2)) if outer_mode == "max" else np.ones(len(values), dtype=np.float32)
    outer = np.where(outer == 0, np.float32(1), outer).astype(np.float32)
    chosen = []
    for value, count in zip(source, rows):
        raw = value[:count].astype("<c8").tobytes()
        blob = _frame(raw, int(count), values.shape[2], 1., 0, 0, compress)
        chosen.append(Encoded(blob, {"header": HEADER.size, "format_tags": 0, "micro_scales": 0,
                                      "codes": len(raw), "exceptions": 0, "raw_body": len(raw),
                                      "stored_body": len(blob)-HEADER.size}, {"exact32": int(count)*values.shape[2]*2}))
    candidates, diagnostics = {}, {}
    # Select micro formats using equal absolute per-component error budgets.
    valid_micro = np.repeat(valid[:, :, :, None], 2, axis=3).transpose(0, 2, 3, 1).reshape(len(values), -1, 32)
    real_source = _micro(source).astype(np.float64)
    tags = np.zeros(real_source.shape[:2], dtype=np.uint8)
    micro_best = np.full(tags.shape, np.inf)
    for fid in FORMATS:
        codes, scales, decoded, exceptions = _candidate(source, protected, outer, fid)
        candidates[fid] = (codes, scales, decoded, exceptions)
        passing = _passing(source, decoded, valid, bounds)
        diagnostics[FORMATS[fid][0]] = int(passing.sum())
        for i in np.flatnonzero(passing):
            uniform_tags = np.full(tags.shape[1], fid, dtype=np.uint8)
            blob, account = _serialize(source[i], codes[i], scales[i], uniform_tags, float(outer[i]), exceptions[i], int(rows[i]), compress)
            if len(blob) < len(chosen[i].blob):
                chosen[i] = Encoded(blob, account, {FORMATS[fid][0]: int(rows[i])*values.shape[2]*2})
        if granularity == "micro":
            error = np.zeros(real_source.shape, dtype=np.float64)
            np.subtract(_micro(decoded), real_source, out=error, where=valid_micro, dtype=np.float64)
            n = np.maximum(1, valid_micro.sum(axis=-1))
            mse = np.sum(error**2, axis=-1)/n
            okay = np.all(np.isfinite(error), axis=-1)
            if bounds.get("max_component_mse") is not None:
                okay &= mse <= bounds["max_component_mse"]
            if bounds.get("max_complex_error") is not None:
                okay &= np.max(np.abs(error), axis=-1) <= bounds["max_complex_error"]/np.sqrt(2)
            if bounds.get("max_complex_bias") is not None:
                okay &= np.abs(error.sum(axis=-1)/n) <= bounds["max_complex_bias"]/np.sqrt(2)
            bits = FORMATS[fid][1]
            previous_bits = np.asarray([32, 4, 6, 6, 8, 8, 16])[tags]
            take = okay & ((bits < previous_bits) | ((bits == previous_bits) & (mse < micro_best)))
            tags[take], micro_best[take] = fid, mse[take]
    if granularity == "micro":
        mixed_codes = np.zeros(real_source.shape, dtype=np.uint16)
        mixed_scales = np.full(tags.shape, 127, dtype=np.uint8)
        mixed_real = _micro(source).copy()
        mixed_exceptions = protected.copy()
        for fid, (codes, scales, decoded, exceptions) in candidates.items():
            select = tags == fid
            mixed_codes[select], mixed_scales[select] = codes[select], scales[select]
            mixed_real[select] = _micro(decoded)[select]
            exception_micro = np.repeat(exceptions[:, :, :, None], 2, axis=3).transpose(0, 2, 3, 1).reshape(len(values), -1, 32)
            selected_exception = exception_micro & select[..., None]
            mixed_exceptions |= selected_exception.reshape(len(values), values.shape[2], 2, 32).any(axis=2).transpose(0, 2, 1)
        decoded = _complex(mixed_real)
        decoded[mixed_exceptions] = source[mixed_exceptions]
        passing = _passing(source, decoded, valid, bounds)
        diagnostics["mixed"] = int(passing.sum())
        for i in np.flatnonzero(passing):
            blob, account = _serialize(source[i], mixed_codes[i], mixed_scales[i], tags[i], float(outer[i]), mixed_exceptions[i], int(rows[i]), compress)
            if len(blob) < len(chosen[i].blob):
                ids, counts = np.unique(tags[i], return_counts=True)
                histogram = {"exact32" if f == 0 else FORMATS[int(f)][0]: int(n)*int(rows[i]) for f, n in zip(ids, counts)}
                chosen[i] = Encoded(blob, account, histogram)
    return chosen, diagnostics


def decode_block(blob):
    if len(blob) < HEADER.size:
        raise ValueError("Truncated MX header")
    magic, zipped, fid, layout, reserved, rows, channels, outer, count, rawlen, storedlen, crc = HEADER.unpack_from(blob)
    payload = blob[HEADER.size:]
    if (magic != b"MXB1" or zipped not in (0, 1) or fid not in (0, 255, *FORMATS)
            or layout or reserved or not 1 <= rows <= 32 or not 1 <= channels <= 24
            or not np.isfinite(outer) or outer <= 0 or count > rows*channels
            or rawlen > MAX_RAW or storedlen > MAX_RAW or len(payload) != storedlen):
        raise ValueError("Invalid MX header/length")
    if zlib.crc32(payload, zlib.crc32(blob[:HEADER.size-4])) != crc:
        raise ValueError("MX checksum mismatch")
    if zipped:
        stream = zlib.decompressobj()
        raw = stream.decompress(payload, rawlen+1)
        if not stream.eof or stream.unused_data or stream.unconsumed_tail:
            raise ValueError("Invalid compressed MX body")
    else:
        raw = payload
    if len(raw) != rawlen:
        raise ValueError("MX decoded length mismatch")
    if fid == 0:
        if count or outer != 1 or len(raw) != rows*channels*8:
            raise ValueError("Invalid exact MX block")
        return np.frombuffer(raw, dtype="<c8").reshape(rows, channels).copy()
    micros = channels*2
    pos = micros if fid == 255 else 0
    tags = np.frombuffer(raw[:pos], dtype=np.uint8) if fid == 255 else np.full(micros, fid, dtype=np.uint8)
    if len(tags) != micros or any(int(tag) not in (0, *FORMATS) for tag in tags):
        raise ValueError("Invalid micro format tags")
    scales = np.frombuffer(raw[pos:pos+micros], dtype=np.uint8)
    if len(scales) != micros or np.any(scales == 255):
        raise ValueError("Invalid micro scales")
    pos += micros
    decoded = np.empty((micros, 32), dtype=np.float32)
    for tag in sorted(np.unique(tags)):
        selected = tags == tag
        n = int(selected.sum())*32
        width = 32 if tag == 0 else FORMATS[int(tag)][1]
        length = n*width//8
        if pos+length > len(raw):
            raise ValueError("Truncated MX codes")
        if tag == 0:
            decoded[selected] = np.frombuffer(raw[pos:pos+length], dtype="<f4").reshape(-1, 32)
        else:
            codes = unpack_codes(raw[pos:pos+length], width, n).reshape(1, -1, 32)
            exponent = scales[selected].astype(np.int16)[None] - 127
            decoded[selected] = _reconstruct(codes, exponent, int(tag), np.asarray([outer]))[0]
        pos += length
    if pos+count*10 != len(raw):
        raise ValueError("Invalid MX exception byte count")
    result = _complex(decoded[None])[0, :rows].copy()
    positions = np.frombuffer(raw[pos:pos+count*2], dtype="<u2")
    if np.any(positions >= rows*channels) or (len(positions)>1 and np.any(positions[1:] <= positions[:-1])):
        raise ValueError("Invalid MX exception positions")
    result.ravel()[positions] = np.frombuffer(raw[pos+count*2:], dtype="<c8")
    return result


def decode_component(path):
    """Decode a complete experimental DATA component using only its contents."""
    import h5py
    from .codec import decode_block as uniform_decode
    with h5py.File(path, 'r') as h5:
        if (h5.attrs.get('format') != 'varquanms.experimental-mx-component/1'
                or not h5.attrs.get('complete') or not h5.attrs.get('full_observation')
                or h5.attrs.get('variant') not in ('uniform','unit-mega','max-mega','unit-micro','max-micro')):
            raise ValueError('Not a complete experimental visibility component')
        shape = tuple(int(x) for x in h5.attrs['shape'])
        if len(shape) != 3 or any(x <= 0 for x in shape) or np.prod(shape, dtype=object)*8 > 2*1024**3:
            raise ValueError('Invalid component shape')
        permutations = [h5[f'permutations/{p}'][:] for p in range(shape[2])]
        for permutation in permutations:
            if not np.array_equal(np.sort(permutation), np.arange(shape[0])):
                raise ValueError('Invalid component permutation')
        records, offsets = h5['records'][:], h5['offsets'][:]
        payload = h5['payload'][:].tobytes()
        if (records.ndim != 2 or records.shape[1] != 5 or len(offsets) != len(records)+1
                or offsets[0] != 0 or offsets[-1] != len(payload) or np.any(offsets[1:] <= offsets[:-1])):
            raise ValueError('Invalid component index')
        result = np.empty(shape, dtype=np.complex64)
        present = np.zeros(shape, dtype=bool)
        for record, left, right in zip(records, offsets[:-1], offsets[1:]):
            p,a,b,lo,hi = map(int,record)
            if not (0 <= p < shape[2] and 0 <= a < b <= shape[0] and 0 <= lo < hi <= shape[1]):
                raise ValueError('Invalid component block bounds')
            rows = permutations[p][a:b]
            if np.any(present[rows,lo:hi,p]):
                raise ValueError('Overlapping component records')
            blob = payload[int(left):int(right)]
            block = uniform_decode(blob,(b-a,hi-lo)) if h5.attrs['variant']=='uniform' else decode_block(blob)
            if block.shape != (b-a,hi-lo):
                raise ValueError('Component block shape mismatch')
            result[rows,lo:hi,p] = block
            present[rows,lo:hi,p] = True
        if not np.all(present):
            raise ValueError('Missing component records')
        return result
