"""uniform-symmetric-v1: specified packing, exceptions and exact escape blocks.

Header <BBHIIf: component bits (0 = exact complex64), body compression
(0 = none, 1 = zlib), reserved uint16 = 0, complex sample count uint32,
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

import numpy as np

from .config import SUPPORTED_BITS
from .metadata import portable

HEADER = struct.Struct("<BBHIIf")
MAX_SAMPLES = 262144
ALLOWED_BITS = SUPPORTED_BITS


def pack_codes(codes, bits):
    if bits not in ALLOWED_BITS:
        raise ValueError("Unsupported packing width")
    a = np.asarray(codes)
    if a.dtype.kind not in "iu" or np.any(a < 0) or np.any(a >= 2**bits):
        raise ValueError("Code is outside its bit width")
    a = a.astype(np.uint16, copy=False).ravel()
    unpacked = ((a[:, None] >> np.arange(bits, dtype=np.uint16)) & 1).ravel()
    return np.packbits(unpacked, bitorder="little").tobytes()


def unpack_codes(payload, bits, count):
    if bits not in ALLOWED_BITS or type(count) is not int or count < 0 or count > 2 * MAX_SAMPLES:
        raise ValueError("Invalid packing width/count")
    expected = (count * bits + 7) // 8
    if len(payload) != expected:
        raise ValueError("Packed code length mismatch")
    remainder = count * bits % 8
    if remainder and payload[-1] >> remainder:
        raise ValueError("Nonzero code padding bits")
    values = np.unpackbits(np.frombuffer(payload, dtype=np.uint8), bitorder="little", count=count * bits).reshape(count, bits)
    dtype = np.uint8 if bits <= 8 else np.uint16
    return np.sum(values.astype(np.uint16) << np.arange(bits, dtype=np.uint16), axis=1).astype(dtype)


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


def decode_block(blob, shape):
    if len(blob) < HEADER.size:
        raise ValueError("Truncated block header")
    bits, compressed, reserved, count, exceptions, step = HEADER.unpack_from(blob)
    if reserved != 0 or count != math.prod(shape) or not 0 < count <= MAX_SAMPLES:
        raise ValueError("Invalid block sample count/header")
    if bits == 0:
        if exceptions or step != 0:
            raise ValueError("Invalid exact block header")
        return np.frombuffer(_body(blob[HEADER.size:], compressed, count * 8), dtype="<c8").copy().reshape(shape)
    if bits not in ALLOWED_BITS or not math.isfinite(step) or step < 0 or exceptions > count:
        raise ValueError("Invalid quantized block header")
    packed_bytes = (2 * count * bits + 7) // 8
    raw = _body(blob[HEADER.size:], compressed, packed_bytes + exceptions * 12)
    codes = unpack_codes(raw[:packed_bytes], bits, 2 * count)
    if np.any(codes == 2**bits - 1):
        raise ValueError("Reserved uniform code encountered")
    maximum = 2**(bits - 1) - 1
    signed = codes.astype(np.int32) - maximum
    with np.errstate(over="ignore", invalid="ignore"):
        components = signed.astype(np.float32) * np.float32(step)
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


def encode_block(original, flags, config, bounds, stream=()):
    total_start = time.perf_counter()
    times = Counter()
    original = portable(original)
    if original.dtype != np.dtype("<c8") or original.ndim != 2 or not 0 < original.size <= MAX_SAMPLES or flags.shape != original.shape:
        raise ValueError("Codec expects bounded complex64[baseline,channel] blocks")
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
    signed_zero = np.any((components == 0) & np.signbit(components), axis=-1)
    mandatory = flags | ~np.isfinite(original) | signed_zero
    estimable = valid & ~signed_zero
    magnitudes = np.max(np.abs(components.astype(np.float64)), axis=-1)[estimable]
    cutoffs = [float(magnitudes.max()) if q == 1.0 else float(np.quantile(magnitudes, q))
               for q in config.scale_quantiles] if magnitudes.size else [0.0] * len(config.scale_quantiles)
    times["quantization_and_scale_estimation"] += time.perf_counter() - stamp
    for bits in sorted(config.bits):
        maximum = 2**(bits - 1) - 1
        for q_index, quantile in enumerate(config.scale_quantiles):
            stamp = time.perf_counter()
            cutoff = cutoffs[q_index]
            step = np.float32(cutoff / maximum)
            if float(step) * maximum < cutoff:
                step = np.nextafter(step, np.float32(np.inf))
            if not np.isfinite(step):
                times["quantization_and_scale_estimation"] += time.perf_counter() - stamp
                continue
            limit = float(step) * maximum
            outside = np.any(np.abs(components.astype(np.float64)) > limit, axis=-1)
            exception_mask = mandatory | outside
            safe = np.where(exception_mask[..., None], np.float32(0), components)
            scaled = safe.astype(np.float64) / float(step) if step else np.zeros_like(safe, dtype=np.float64)
            if config.rounding == "nearest-even":
                rounded = np.rint(scaled)
            else:
                rounded = np.empty_like(scaled)
                for component in range(2):
                    seed_words = [config.seed & 0xFFFFFFFF, config.seed >> 32, *stream, bits, q_index, component]
                    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed_words)))
                    low = np.floor(scaled[..., component])
                    rounded[..., component] = low + (rng.random(low.shape) < scaled[..., component] - low)
            if np.any(np.abs(rounded) > maximum):
                times["quantization_and_scale_estimation"] += time.perf_counter() - stamp
                continue
            codes = (rounded.astype(np.int32) + maximum).astype(np.uint16).ravel()
            positions = np.flatnonzero(exception_mask).astype("<u4")
            values = original.ravel()[positions]
            times["quantization_and_scale_estimation"] += time.perf_counter() - stamp
            stamp = time.perf_counter()
            packed = pack_codes(codes, bits)
            raw_body = packed + positions.tobytes() + values.tobytes()
            compression, body = _compress(raw_body, config)
            blob = HEADER.pack(bits, compression, 0, original.size, len(positions), step) + body
            times["packing_and_secondary_compression"] += time.perf_counter() - stamp
            if len(blob) >= len(best.blob):
                continue
            try:
                stamp = time.perf_counter()
                decoded = decode_block(blob, original.shape)
                metrics = error_metrics(original, decoded, valid, percentile=bounds["max_p99_error"] is not None)
            except ValueError:
                continue
            finally:
                times["candidate_reconstruction_and_error_checks"] += time.perf_counter() - stamp
            if not within_bounds(metrics, bounds):
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
    times["codec_total"] = time.perf_counter() - total_start
    best.timings = dict(times)
    return best
