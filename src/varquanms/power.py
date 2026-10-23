"""power-midrise-v1: median-centered signed power levels, with exact exceptions.

Mapping 2 uses the base header's last float32 as residual range t, followed by
<fff: real median, imaginary median, gamma. Gamma is exactly 1, 1/2 or 1/4.
Decoder powers are fixed binary64 products, without libm pow or fused arithmetic.
"""
import math
import struct
import time

import numpy as np

from . import codec

PARAMETERS = struct.Struct("<fff")
HEADER_BYTES = codec.HEADER.size + PARAMETERS.size
GAMMAS = (1., .5, .25)


def reconstruct(codes, bits, radius, center, gamma):
    """IEEE binary64 divide, products, add; one final float32 rounding."""
    if radius == 0:
        shape = np.broadcast_shapes(np.shape(codes), np.shape(center))
        return np.broadcast_to(np.asarray(center, dtype=np.float32), shape).copy()
    u = (2 * np.asarray(codes, dtype=np.float64) - (2**bits-1)) / (2**bits-1)
    magnitude = np.abs(u)
    if gamma != 1:
        magnitude = magnitude * magnitude
        if gamma == .25:
            magnitude = magnitude * magnitude
    residual = np.copysign(magnitude, u) * float(radius)
    with np.errstate(over="ignore", invalid="ignore"):
        return (residual + np.asarray(center, dtype=np.float64)).astype(np.float32)


def decode_power(blob, shape):
    if len(blob) < HEADER_BYTES:
        raise ValueError("Truncated power block metadata")
    bits, compressed, mapping, count, exceptions, radius = codec.HEADER.unpack_from(blob)
    cr, ci, gamma = PARAMETERS.unpack_from(blob, codec.HEADER.size)
    if mapping != 2 or bits not in codec.ALLOWED_BITS or count != math.prod(shape) or not 0 < count <= codec.MAX_SAMPLES or exceptions > count or not math.isfinite(radius) or radius < 0 or not math.isfinite(cr) or not math.isfinite(ci) or gamma not in GAMMAS:
        raise ValueError("Invalid power block header/parameters")
    packed_bytes = (2*count*bits+7)//8
    raw = codec._body(blob[HEADER_BYTES:], compressed, packed_bytes+12*exceptions)
    codes = codec.unpack_codes(raw[:packed_bytes], bits, 2*count).reshape(-1, 2)
    components = reconstruct(codes, bits, radius, [cr, ci], gamma)
    if not np.all(np.isfinite(components)):
        raise ValueError("Power reconstruction overflow")
    result = components.astype("<f4", copy=False).view("<c8").reshape(-1)
    if exceptions:
        positions = np.frombuffer(raw, dtype="<u4", count=exceptions, offset=packed_bytes)
        if np.any(positions >= count) or np.any(positions[1:] <= positions[:-1]):
            raise ValueError("Invalid power exception positions")
        values = np.frombuffer(raw, dtype="<c8", count=exceptions, offset=packed_bytes+4*exceptions)
        result[positions] = values
    return result.reshape(shape)


def _parameters(original, flags, config):
    components = np.stack([original.real, original.imag], axis=-1).astype(np.float64)
    signed_zero = np.any((components == 0) & np.signbit(components), axis=-1)
    mandatory = flags | ~np.isfinite(original) | signed_zero
    eligible = components[~mandatory]
    center = (np.median(eligible, axis=0).astype(np.float32) if eligible.size and config.power_center == "median"
              else np.zeros(2, dtype=np.float32))
    # Evaluate residuals against the already rounded/stored median.
    residuals = eligible-center.astype(np.float64)
    ranges = np.max(np.abs(residuals), axis=-1)
    cutoffs = [float(ranges.max()) if q == 1 else float(np.quantile(ranges, q))
               for q in config.scale_quantiles] if ranges.size else [0.] * len(config.scale_quantiles)
    widths, parameters = [], []
    for bits in sorted(config.bits):
        for cutoff in cutoffs:
            with np.errstate(over="ignore"):
                radius = np.float32(cutoff)
                if float(radius) < cutoff:
                    radius = np.nextafter(radius, np.float32(np.inf))
            if not np.isfinite(radius):
                continue
            for gamma in config.power_gammas:
                widths.append(bits); parameters.append([radius, *center, gamma])
    return components, mandatory, signed_zero, np.asarray(widths, dtype="u1"), np.asarray(parameters, dtype="<f4").reshape(-1,4)


def _record(blob, original, flags, mandatory, signed_zero, timings):
    decoded = codec.decode_block(blob, original.shape, "power-midrise-v1")
    metrics = codec.error_metrics(original, decoded, ~flags & np.isfinite(original))
    bits, _, _, count, exceptions, radius = codec.HEADER.unpack_from(blob)
    parameters = {}
    if bits:
        cr, ci, gamma = PARAMETERS.unpack_from(blob, codec.HEADER.size)
        parameters = {"center_real": cr, "center_imag": ci, "radius": radius, "gamma": gamma}
    accounting = {"lossless_component_bytes": count*8 if not bits else 0,
                  "packed_code_bytes": (2*count*bits+7)//8,
                  "exception_position_bytes": exceptions*4, "exception_value_bytes": exceptions*8,
                  "header_bytes": HEADER_BYTES if bits else codec.HEADER.size,
                  "flagged_exception_samples": int(flags.sum()) if bits else 0,
                  "nonfinite_exception_samples": int(np.count_nonzero(~flags & ~np.isfinite(original))) if bits else 0,
                  "signed_zero_exception_samples": int(np.count_nonzero(~flags & np.isfinite(original) & signed_zero)) if bits else 0,
                  "outlier_exception_samples": exceptions-int(mandatory.sum()) if bits else 0}
    return codec.EncodedBlock(blob, decoded, bits, exceptions, metrics, accounting, timings, parameters)


def encode_power(original, flags, config, bounds):
    started = time.perf_counter()
    if config.rounding != "nearest-even":
        raise ValueError("Power codec supports nearest reconstructed level, with even-code ties")
    components, mandatory, signed_zero, widths, parameters = _parameters(original, flags, config)
    timings = {"quantization_and_scale_estimation": time.perf_counter()-started}
    native = codec._native
    if native is not None and getattr(native, "POWER_V1", 0) and 0 < len(widths) <= 4096 and bounds["max_p99_error"] is None:
        limits = [bounds[k] if bounds[k] is not None else np.inf for k in
                  ("max_component_mse", "max_complex_error", "max_complex_bias", "max_baseline_component_mse")]
        stamp = time.perf_counter()
        blob, reconstructed, quant_cpu, compress_cpu = native.power_search(original, np.ascontiguousarray(flags, dtype=bool),
            widths, parameters, np.asarray(limits, dtype="<f8"), original.shape[1],
            config.compression_level if config.compression == "gzip" else -1)
        timings.update(native_candidate_search=time.perf_counter()-stamp,
                       native_search_quantization_cpu=quant_cpu, native_search_compression_cpu=compress_cpu)
        candidate = _record(blob, original, flags, mandatory, signed_zero, timings)
        if candidate.decoded.tobytes() != reconstructed:
            raise ValueError("Native power reconstruction differs from reference")
        if codec.within_bounds(candidate.metrics, bounds):
            timings["codec_total"] = time.perf_counter()-started
            return candidate
    stamp = time.perf_counter()
    compression, body = codec._compress(original.tobytes(), config)
    best_blob = codec.HEADER.pack(0, compression, 0, original.size, 0, 0.) + body
    valid = ~flags & np.isfinite(original)
    for bits, (radius, cr, ci, gamma) in zip(widths, parameters):
        bits = int(bits); maximum = 2**bits-1
        center = np.asarray([cr, ci], dtype=np.float64)
        with np.errstate(invalid="ignore"):
            residuals = components-center
            outside = np.any(np.abs(residuals) > float(radius), axis=-1)
        exception_mask = mandatory | outside
        safe = np.where(exception_mask[...,None], 0., residuals)
        u = safe/float(radius) if radius else np.zeros_like(safe)
        magnitude = np.abs(u)
        if gamma != 1:
            magnitude = np.sqrt(magnitude)
            if gamma == .25:
                magnitude = np.sqrt(magnitude)
        coordinate = (np.copysign(magnitude, u)+1)*(maximum/2)
        low = np.floor(coordinate).astype(np.int64)
        # Range exceptions happen before this; clipping only protects endpoints
        # from one-ulp coordinate arithmetic and cannot saturate input outliers.
        low = np.clip(low, 0, maximum); high = np.minimum(low+1, maximum)
        lower = reconstruct(low, bits, radius, center, gamma)
        upper = reconstruct(high, bits, radius, center, gamma)
        target = np.where(exception_mask[...,None], center, components)
        dl = np.abs(lower.astype(np.float64)-target); dh = np.abs(upper.astype(np.float64)-target)
        choose_high = (dh < dl) | ((dh == dl) & (high % 2 == 0))
        codes = np.where(choose_high, high, low).astype(np.uint16)
        if not radius:
            codes.fill(2**(bits-1))
        decoded_components = reconstruct(codes, bits, radius, center, gamma)
        if not np.all(np.isfinite(decoded_components)):
            continue
        decoded = decoded_components.astype("<f4", copy=False).view("<c8").reshape(original.shape)
        positions = np.flatnonzero(exception_mask).astype("<u4")
        values = original.ravel()[positions]; decoded.ravel()[positions] = values
        metrics = codec.error_metrics(original, decoded, valid, percentile=bounds["max_p99_error"] is not None)
        if not codec.within_bounds(metrics, bounds):
            continue
        raw = codec.pack_codes(codes.ravel(), bits)+positions.tobytes()+values.tobytes()
        compression, body = codec._compress(raw, config)
        blob = codec.HEADER.pack(bits, compression, 2, original.size, len(positions), radius) + PARAMETERS.pack(cr,ci,gamma) + body
        if len(blob) < len(best_blob):
            best_blob = blob
    timings["reference_candidate_search"] = time.perf_counter()-stamp
    best = _record(best_blob, original, flags, mandatory, signed_zero, timings)
    if not codec.within_bounds(best.metrics, bounds):
        raise ValueError("Selected power reconstruction exceeds error bounds")
    timings["codec_total"] = time.perf_counter()-started
    return best
