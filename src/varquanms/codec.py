"""Versioned uniform codecs: packing, exceptions and exact escape blocks.

Header <BBHIIf: component bits (0 = exact complex64), body compression
(0 = none, 1 = zlib), mapping uint16 (0 = symmetric, 1 = midrise), sample count uint32,
exception count uint32, positive float32 quantization step. All little endian.
The raw body is either complex64 bits, or packed component codes followed by
uint32 exception positions and original complex64 bits. See docs/format.md.
"""

from collections import Counter
from dataclasses import dataclass, field
import math
import struct
import time
import zlib
import sys

import numpy as np

from .config import CODEBOOK_IDS, SUPPORTED_BITS
from .metadata import portable

HEADER = struct.Struct("<BBHIIf")
MAX_SAMPLES = 262144
ALLOWED_BITS = SUPPORTED_BITS
try:
    from . import _native
except ImportError:
    _native = None
if sys.byteorder != "little":
    _native = None


def pack_codes(codes, bits):
    if bits not in ALLOWED_BITS:
        raise ValueError("Unsupported packing width")
    a = np.asarray(codes)
    if a.dtype.kind not in "iu" or np.any(a < 0) or np.any(a >= 2**bits):
        raise ValueError("Code is outside its bit width")
    a = a.astype(np.uint16, copy=False).ravel()
    if bits == 8:
        return a.astype(np.uint8).tobytes()
    if bits == 16:
        return a.astype("<u2", copy=False).tobytes()
    if bits in (2, 4):
        per_byte = 8 // bits
        padded = np.zeros(((a.size + per_byte - 1) // per_byte) * per_byte, dtype=np.uint8)
        padded[:a.size] = a
        result = padded[::per_byte].copy()
        for j in range(1, per_byte):
            result |= padded[j::per_byte] << (bits * j)
        return result.tobytes()
    # 3/6/12-bit streams repeat after eight/four/two codes respectively. Work in
    # bytes, without the old samples-by-bits temporary array.
    count = (a.size * bits + 7) // 8
    result = np.zeros(count, dtype=np.uint8)
    per_group = 8 // math.gcd(bits, 8)
    for j in range(per_group):
        values = a[j::per_group]
        offset, shift = divmod(j * bits, 8)
        stride = per_group * bits // 8
        result[offset:offset + len(values)*stride:stride] |= (values << shift).astype(np.uint8)
        if shift + bits > 8:
            high = values >> (8 - shift)
            target = result[offset+1:offset+1+len(values)*stride:stride]
            target |= high[:len(target)].astype(np.uint8)
    return result.tobytes()


def unpack_codes(payload, bits, count):
    if bits not in ALLOWED_BITS or type(count) is not int or count < 0 or count > 2 * MAX_SAMPLES:
        raise ValueError("Invalid packing width/count")
    expected = (count * bits + 7) // 8
    if len(payload) != expected:
        raise ValueError("Packed code length mismatch")
    remainder = count * bits % 8
    if remainder and payload[-1] >> remainder:
        raise ValueError("Nonzero code padding bits")
    if bits == 8:
        return np.frombuffer(payload, dtype=np.uint8).copy()
    if bits == 16:
        return np.frombuffer(payload, dtype="<u2").astype(np.uint16, copy=True)
    source = np.frombuffer(payload, dtype=np.uint8)
    result = np.empty(count, dtype=np.uint8 if bits <= 8 else np.uint16)
    per_group = 8 // math.gcd(bits, 8)
    stride = per_group * bits // 8
    for j in range(per_group):
        target = result[j::per_group]
        offset, shift = divmod(j * bits, 8)
        values = source[offset:offset+len(target)*stride:stride].astype(np.uint16) >> shift
        if shift + bits > 8:
            high = source[offset+1:offset+1+len(target)*stride:stride].astype(np.uint16)
            values[:len(high)] |= high << (8-shift)
        target[:] = values & (2**bits - 1)
    return result


def error_metrics(original, decoded, valid, percentile=True):
    if original.shape != decoded.shape or original.shape != valid.shape:
        raise ValueError("Error metric shape mismatch")
    valid = valid & np.isfinite(original)
    if np.any(~np.isfinite(decoded[valid])):
        raise ValueError("Nonfinite reconstruction")
    r = np.zeros(original.shape, dtype=np.float64)
    i = np.zeros_like(r)
    np.subtract(decoded.real, original.real, out=r, where=valid, dtype=np.float64)
    np.subtract(decoded.imag, original.imag, out=i, where=valid, dtype=np.float64)
    count = int(np.count_nonzero(valid))
    squared = r * r + i * i
    absolute = np.hypot(r[valid], i[valid])
    by_row_count = valid.sum(axis=1)
    by_row = squared.sum(axis=1) / (2 * np.maximum(1, by_row_count))
    mean_r, mean_i = float(r.sum() / max(1, count)), float(i.sum() / max(1, count))
    return {
        "valid_samples": count,
        "component_mse": float(squared.sum() / (2 * max(1, count))),
        "real_mse": float(np.sum(r * r) / max(1, count)),
        "imag_mse": float(np.sum(i * i) / max(1, count)),
        "bias_real": mean_r, "bias_imag": mean_i,
        "complex_bias_abs": math.hypot(mean_r, mean_i),
        "max_complex_error": float(absolute.max()) if count else 0.0,
        "p99_error": float(np.quantile(absolute, 0.99)) if count and percentile else 0.0,
        "max_baseline_component_mse": float(np.max(by_row)) if by_row.size else 0.0,
    }


def within_bounds(metrics, bounds):
    fields = {"max_component_mse": "component_mse", "max_complex_error": "max_complex_error",
              "max_complex_bias": "complex_bias_abs", "max_baseline_component_mse": "max_baseline_component_mse",
              "max_p99_error": "p99_error"}
    return all(limit is None or metrics[fields[name]] <= limit for name, limit in bounds.items())


def _compress(raw, config):
    if config.compression == "gzip":
        compressed = zlib.compress(raw, config.compression_level)
        if len(compressed) < len(raw):
            return 1, compressed
    return 0, raw


def _body(payload, compressed, expected):
    if compressed == 0:
        raw = payload
    elif compressed == 1:
        decoder = zlib.decompressobj()
        raw = decoder.decompress(payload, expected + 1)
        if not decoder.eof or decoder.unused_data or decoder.unconsumed_tail:
            raise ValueError("Truncated, oversized or trailing compressed block")
    else:
        raise ValueError("Unknown block compression identifier")
    if len(raw) != expected:
        raise ValueError("Decoded block byte count mismatch")
    return raw


def decode_block(blob, shape, codebook=None):
    if len(blob) < HEADER.size:
        raise ValueError("Truncated block header")
    bits, compressed, mapping, count, exceptions, step = HEADER.unpack_from(blob)
    if mapping not in CODEBOOK_IDS.values() or count != math.prod(shape) or not 0 < count <= MAX_SAMPLES:
        raise ValueError("Invalid block sample count/header")
    if codebook is not None and (codebook not in CODEBOOK_IDS or (bits and mapping != CODEBOOK_IDS[codebook])):
        raise ValueError("Block mapping disagrees with declared codebook")
    if mapping == 2:
        from .power import decode_power
        return decode_power(blob, shape)
    if bits == 0:
        if mapping or exceptions or step != 0:
            raise ValueError("Invalid exact block header")
        return np.frombuffer(_body(blob[HEADER.size:], compressed, count * 8), dtype="<c8").copy().reshape(shape)
    if bits not in ALLOWED_BITS or not math.isfinite(step) or step < 0 or exceptions > count:
        raise ValueError("Invalid quantized block header")
    packed_bytes = (2 * count * bits + 7) // 8
    raw = _body(blob[HEADER.size:], compressed, packed_bytes + exceptions * 12)
    codes = unpack_codes(raw[:packed_bytes], bits, 2 * count)
    if mapping == 0 and np.any(codes == 2**bits - 1):
        raise ValueError("Reserved uniform code encountered")
    center = (2**bits - 1) / 2 if mapping else 2**(bits - 1) - 1
    coefficients = codes.astype(np.float32) - np.float32(center)
    with np.errstate(over="ignore", invalid="ignore"):
        components = coefficients * np.float32(step)
    if not np.all(np.isfinite(components)):
        raise ValueError("Quantized reconstruction overflow")
    result = components.astype("<f4", copy=False).view("<c8").reshape(-1)
    if exceptions:
        positions = np.frombuffer(raw, dtype="<u4", count=exceptions, offset=packed_bytes)
        if np.any(positions >= count) or np.any(positions[1:] <= positions[:-1]):
            raise ValueError("Invalid exception positions")
        values = np.frombuffer(raw, dtype="<c8", count=exceptions, offset=packed_bytes + exceptions * 4)
        result[positions] = values
    return result.reshape(shape)


@dataclass
class EncodedBlock:
    blob: bytes
    decoded: np.ndarray
    bits: int
    exceptions: int
    metrics: dict
    accounting: dict
    timings: dict = field(default_factory=dict)
    parameters: dict = field(default_factory=dict)


def _native_encode_block(original, flags, config, bounds):
    started = stamp = time.perf_counter()
    mapping = CODEBOOK_IDS[config.codebook]
    valid = ~flags & np.isfinite(original)
    signed_zero = ((original.real == 0) & np.signbit(original.real)) | ((original.imag == 0) & np.signbit(original.imag))
    mandatory = flags | ~np.isfinite(original) | signed_zero
    amplitude = np.maximum(np.abs(original.real), np.abs(original.imag))[valid & ~signed_zero].astype(np.float64)
    cutoffs = [float(amplitude.max()) if q == 1 else float(np.quantile(amplitude, q))
               for q in config.scale_quantiles] if amplitude.size else [0.] * len(config.scale_quantiles)
    widths, steps = [], []
    for bits in sorted(config.bits):
        maximum = (2**bits-1)/2 if mapping else 2**(bits-1)-1
        for cutoff in cutoffs:
            with np.errstate(over="ignore"):
                step = np.float32(cutoff/maximum)
                if float(step)*maximum < cutoff:
                    step = np.nextafter(step, np.float32(np.inf))
            widths.append(bits); steps.append(step)
    if len(widths) > 4096:
        return None
    limits = [bounds[name] if bounds[name] is not None else np.inf for name in
              ("max_component_mse", "max_complex_error", "max_complex_bias", "max_baseline_component_mse")]
    timings = {"quantization_and_scale_estimation": time.perf_counter()-stamp}
    stamp = time.perf_counter()
    blob, reconstructed, quantization_cpu, compression_cpu = _native.search(original, np.ascontiguousarray(flags, dtype=bool),
        np.asarray(widths, dtype="u1"), np.asarray(steps, dtype="<f4"), np.asarray(limits, dtype="<f8"),
        original.shape[1], config.compression_level if config.compression == "gzip" else -1,
        *([mapping] if mapping else []))
    timings["native_candidate_search"] = time.perf_counter()-stamp
    timings["native_search_quantization_cpu"] = quantization_cpu
    timings["native_search_compression_cpu"] = compression_cpu
    stamp = time.perf_counter()
    decoded = decode_block(blob, original.shape, config.codebook)
    if decoded.tobytes() != reconstructed:
        raise ValueError("Native candidate differs from the reference decoder")
    metrics = error_metrics(original, decoded, valid)
    timings["candidate_reconstruction_and_error_checks"] = time.perf_counter()-stamp
    if not within_bounds(metrics, bounds):
        # A native reduction very near a threshold can differ from NumPy.
        # Repeat the complete reference selection, keeping all bounds strict.
        return None
    bits, _, _, count, exceptions, _ = HEADER.unpack_from(blob)
    accounting = {"lossless_component_bytes": count*8 if bits == 0 else 0,
                  "packed_code_bytes": (2*count*bits+7)//8,
                  "exception_position_bytes": exceptions*4, "exception_value_bytes": exceptions*8,
                  "header_bytes": HEADER.size,
                  "flagged_exception_samples": int(flags.sum()) if bits else 0,
                  "nonfinite_exception_samples": int(np.count_nonzero(~flags & ~np.isfinite(original))) if bits else 0,
                  "signed_zero_exception_samples": int(np.count_nonzero(~flags & np.isfinite(original) & signed_zero)) if bits else 0,
                  "outlier_exception_samples": exceptions-int(mandatory.sum()) if bits else 0}
    timings["codec_total"] = time.perf_counter()-started
    return EncodedBlock(blob, decoded, bits, exceptions, metrics, accounting, timings)


def encode_block(original, flags, config, bounds, stream=(), *, global_std=None):
    total_start = time.perf_counter()
    times = Counter()
    original = portable(original)
    if original.dtype != np.dtype("<c8") or original.ndim != 2 or not 0 < original.size <= MAX_SAMPLES or flags.shape != original.shape:
        raise ValueError("Codec expects bounded complex64[baseline,channel] blocks")
    from .error_budget import configured_bounds
    bounds = configured_bounds(original, flags, config, bounds, global_std)
    mapping = CODEBOOK_IDS[config.codebook]
    native_supported = (_native is not None and (not mapping or getattr(_native, "MIDRISE_V1", 0) == 1)
                        and (3 not in config.bits or getattr(_native, "THREE_BIT_V1", 0) == 1))
    if mapping == 2:
        from .power import encode_power
        return encode_power(original, flags, config, bounds)
    if native_supported and config.rounding == "nearest-even" and bounds["max_p99_error"] is None:
        accelerated = _native_encode_block(original, flags, config, bounds)
        if accelerated is not None:
            return accelerated
    valid = ~flags & np.isfinite(original)
    raw = original.tobytes()
    stamp = time.perf_counter()
    compression, body = _compress(raw, config)
    times["packing_and_secondary_compression"] += time.perf_counter() - stamp
    zero_metrics = {name: 0.0 for name in ("component_mse", "real_mse", "imag_mse", "bias_real", "bias_imag", "complex_bias_abs", "max_complex_error", "p99_error", "max_baseline_component_mse")}
    zero_metrics["valid_samples"] = int(np.count_nonzero(valid))
    exact = EncodedBlock(HEADER.pack(0, compression, 0, original.size, 0, 0.0) + body,
                         original.copy(), 0, 0, zero_metrics,
                         {"lossless_component_bytes": len(raw), "packed_code_bytes": 0,
                          "exception_position_bytes": 0, "exception_value_bytes": 0,
                          "header_bytes": HEADER.size, "flagged_exception_samples": 0,
                          "nonfinite_exception_samples": 0, "signed_zero_exception_samples": 0,
                          "outlier_exception_samples": 0})
    best = exact
    stamp = time.perf_counter()
    components = np.stack([original.real, original.imag], axis=-1)
    absolute_components = np.abs(components.astype(np.float64))
    signed_zero = np.any((components == 0) & np.signbit(components), axis=-1)
    mandatory = flags | ~np.isfinite(original) | signed_zero
    estimable = valid & ~signed_zero
    magnitudes = np.max(absolute_components, axis=-1)[estimable]
    cutoffs = [float(magnitudes.max()) if q == 1.0 else float(np.quantile(magnitudes, q))
               for q in config.scale_quantiles] if magnitudes.size else [0.0] * len(config.scale_quantiles)
    times["quantization_and_scale_estimation"] += time.perf_counter() - stamp
    for bits in sorted(config.bits):
        maximum = (2**bits - 1)/2 if mapping else 2**(bits - 1) - 1
        for q_index, quantile in enumerate(config.scale_quantiles):
            stamp = time.perf_counter()
            cutoff = cutoffs[q_index]
            step = np.float32(cutoff / maximum)
            if float(step) * maximum < cutoff:
                step = np.nextafter(step, np.float32(np.inf))
            if not np.isfinite(step):
                times["quantization_and_scale_estimation"] += time.perf_counter() - stamp
                continue
            if native_supported and config.rounding == "nearest-even":
                candidate = _native.candidate(original, np.ascontiguousarray(flags, dtype=bool), bits, float(step), original.shape[1],
                                              *([mapping] if mapping else []))
                times["quantization_and_scale_estimation"] += time.perf_counter() - stamp
                if candidate is None:
                    continue
                decoded_raw, packed, position_raw, value_raw, mse, maximum_error, row_mse = candidate
                # Positive-sum native prechecks have a generous rounding margin.
                # Exact NumPy checks below decide every candidate that can win;
                # bias/percentile are never rejected by this prefilter.
                if any(limit is not None and measured > limit * (1 + 1e-9)
                       for measured, limit in ((mse, bounds["max_component_mse"]),
                                               (maximum_error, bounds["max_complex_error"]),
                                               (row_mse, bounds["max_baseline_component_mse"]))):
                    continue
                stamp = time.perf_counter()
                compression, body = _compress(packed + position_raw + value_raw, config)
                exception_count = len(position_raw) // 4
                blob = HEADER.pack(bits, compression, mapping, original.size, exception_count, step) + body
                times["packing_and_secondary_compression"] += time.perf_counter() - stamp
                if len(blob) >= len(best.blob):
                    continue
                stamp = time.perf_counter()
                decoded = np.frombuffer(decoded_raw, dtype="<c8").reshape(original.shape)
                metrics = error_metrics(original, decoded, valid, percentile=bounds["max_p99_error"] is not None)
                times["candidate_reconstruction_and_error_checks"] += time.perf_counter() - stamp
                if not within_bounds(metrics, bounds):
                    continue
                mandatory_count = int(np.count_nonzero(mandatory))
                best = EncodedBlock(blob, decoded, bits, exception_count, metrics, {
                    "lossless_component_bytes": 0, "packed_code_bytes": len(packed),
                    "exception_position_bytes": len(position_raw), "exception_value_bytes": len(value_raw),
                    "header_bytes": HEADER.size, "flagged_exception_samples": int(np.count_nonzero(flags)),
                    "nonfinite_exception_samples": int(np.count_nonzero(~flags & ~np.isfinite(original))),
                    "signed_zero_exception_samples": int(np.count_nonzero(~flags & np.isfinite(original) & signed_zero)),
                    "outlier_exception_samples": exception_count - mandatory_count})
                continue
            limit = float(step) * maximum
            outside = np.any(absolute_components > limit, axis=-1)
            exception_mask = mandatory | outside
            safe = np.where(exception_mask[..., None], np.float32(0), components)
            scaled = safe.astype(np.float64) / float(step) if step else np.zeros_like(safe, dtype=np.float64)
            if mapping:
                scaled += maximum
            if step == 0:
                # Degenerate blocks have one exact +0 value; stochastic
                # half-integer rounding must not manufacture negative zeros.
                rounded = np.full_like(scaled, 2**(bits - 1) if mapping else 0)
            elif config.rounding == "nearest-even":
                rounded = np.rint(scaled)
            else:
                rounded = np.empty_like(scaled)
                for component in range(2):
                    seed_words = [config.seed & 0xFFFFFFFF, config.seed >> 32, *stream, bits, q_index, component]
                    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed_words)))
                    low = np.floor(scaled[..., component])
                    rounded[..., component] = low + (rng.random(low.shape) < scaled[..., component] - low)
            if (mapping and (np.any(rounded < 0) or np.any(rounded > 2**bits - 1))) or (not mapping and np.any(np.abs(rounded) > maximum)):
                times["quantization_and_scale_estimation"] += time.perf_counter() - stamp
                continue
            codes = (rounded if mapping else rounded + maximum).astype(np.uint16).ravel()
            positions = np.flatnonzero(exception_mask).astype("<u4")
            values = original.ravel()[positions]
            times["quantization_and_scale_estimation"] += time.perf_counter() - stamp
            # Evaluate the decoder's float32 coefficients and multiplication.
            stamp = time.perf_counter()
            with np.errstate(over="ignore", invalid="ignore"):
                coefficients = rounded.astype(np.float32) if mapping else rounded.astype(np.int32).astype(np.float32)
                if mapping:
                    coefficients -= np.float32(maximum)
                reconstructed = coefficients * step
            if not np.all(np.isfinite(reconstructed)):
                times["candidate_reconstruction_and_error_checks"] += time.perf_counter() - stamp
                continue
            decoded = reconstructed.astype("<f4", copy=False).view("<c8").reshape(original.shape)
            decoded.ravel()[positions] = values
            metrics = error_metrics(original, decoded, valid, percentile=bounds["max_p99_error"] is not None)
            times["candidate_reconstruction_and_error_checks"] += time.perf_counter() - stamp
            if not within_bounds(metrics, bounds):
                continue
            stamp = time.perf_counter()
            packed = pack_codes(codes, bits)
            raw_body = packed + positions.tobytes() + values.tobytes()
            compression, body = _compress(raw_body, config)
            blob = HEADER.pack(bits, compression, mapping, original.size, len(positions), step) + body
            times["packing_and_secondary_compression"] += time.perf_counter() - stamp
            if len(blob) >= len(best.blob):
                continue
            best = EncodedBlock(blob, decoded, bits, len(positions), metrics, {
                "lossless_component_bytes": 0, "packed_code_bytes": len(packed),
                "exception_position_bytes": len(positions) * 4, "exception_value_bytes": len(positions) * 8,
                "header_bytes": HEADER.size,
                "flagged_exception_samples": int(np.count_nonzero(flags)),
                "nonfinite_exception_samples": int(np.count_nonzero(~flags & ~np.isfinite(original))),
                "signed_zero_exception_samples": int(np.count_nonzero(~flags & np.isfinite(original) & signed_zero)),
                "outlier_exception_samples": int(np.count_nonzero(outside & ~mandatory)),
            })
    if best.bits and bounds["max_p99_error"] is None:
        stamp = time.perf_counter()
        best.metrics = error_metrics(original, best.decoded, valid)
        times["candidate_reconstruction_and_error_checks"] += time.perf_counter() - stamp
    # Check the selected serialized representation with the public decoder.
    # The early error check must never conceal packing/serialization mistakes.
    if best.bits:
        stamp = time.perf_counter()
        if decode_block(best.blob, original.shape, config.codebook).tobytes() != best.decoded.tobytes():
            raise ValueError("Selected encoding differs from candidate reconstruction")
        times["candidate_reconstruction_and_error_checks"] += time.perf_counter() - stamp
    times["codec_total"] = time.perf_counter() - total_start
    best.timings = dict(times)
    return best
