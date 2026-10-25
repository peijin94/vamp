"""One-product quantization blocks, invertible orderings and measured quality."""

from collections import Counter, deque
from concurrent.futures import ThreadPoolExecutor
import math
import time
import zlib

import numpy as np

from .codec import _native, decode_block, encode_block, error_metrics, within_bounds
from .config import CODEBOOK_IDS, Config
from .error_budget import configured_bounds
from .metadata import equal
from .ordering import all_orders


def _flags(t):
    flags = np.asarray(t.getcol("FLAG"), dtype=bool)
    flags |= np.asarray(t.getcol("FLAG_ROW"), dtype=bool)[:, None, None]
    return flags


def protected_mask(values, flags):
    return flags | ~np.isfinite(values) | ((values.real == 0) & np.signbit(values.real)) | ((values.imag == 0) & np.signbit(values.imag))


def quality(original, decoded, flags, names, method=None, workers=1):
    if method not in (None, "numpy-v1", "native-moments-v2"):
        raise ValueError("Unsupported quality calculation")
    if _native is not None and method != "numpy-v1" and 0 < original.shape[1] <= 65536 and original.shape[2] <= 32:
        return _native_quality(original, decoded, flags, names, workers)
    return _numpy_quality(original, decoded, flags, names)


def _native_quality(original, decoded, flags, names, workers=1):
    rows, channels, products = original.shape
    if original.dtype != np.dtype("<c8") or decoded.dtype != original.dtype or decoded.shape != original.shape or flags.shape != original.shape:
        raise ValueError("Quality expects matching complex64 visibility and flag arrays")
    def calculate(start):
        stop = min(rows, start+1024)
        return _native.quality(np.ascontiguousarray(original[start:stop]), np.ascontiguousarray(decoded[start:stop]),
                               np.ascontiguousarray(flags[start:stop], dtype=bool), stop-start, channels, products)
    starts = range(0, rows, 1024)
    if workers == 1:
        pieces = list(map(calculate, starts))
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            pieces = list(pool.map(calculate, starts))
    # Fixed 1024-row partitions and ordered merges make reductions independent
    # of worker count. Retain exactly the valid amplitudes for global p99.
    piece_moments = [np.frombuffer(piece[0], dtype="<f8").reshape(products, 7) for piece in pieces]
    moments = np.sum(piece_moments, axis=0)
    moments[:, 4:6] = np.max(np.stack(piece_moments)[:, :, 4:6], axis=0)
    raw = [None] * 11
    raw[0] = moments.tobytes()
    for key in (1, 2):
        raw[key] = b"".join(piece[key] for piece in pieces)
    for key in (3, 4, 5, 6, 7, 8):
        dtype = "<u8" if key in (4, 6, 8) else "<f8"
        raw[key] = np.sum([np.frombuffer(piece[key], dtype=dtype) for piece in pieces], axis=0).tobytes()
    raw[10] = sum(piece[10] for piece in pieces)
    moments = np.frombuffer(raw[0], dtype="<f8").reshape(products, 7)
    def percentile(p):
        values = np.concatenate([np.frombuffer(piece[9], dtype="<f8").reshape(products, -1)[p, :int(pm[p,6])]
                                 for piece, pm in zip(pieces, piece_moments)])
        return float(np.quantile(values, .99, overwrite_input=True)) if len(values) else 0.0
    if workers == 1:
        percentiles = list(map(percentile, range(products)))
    else:
        with ThreadPoolExecutor(max_workers=min(workers, products)) as pool:
            percentiles = list(pool.map(percentile, range(products)))
    metrics = {}
    for p, name in enumerate(names):
        r, i, rr, ii, maximum, max_row, valid_count = moments[p]
        count = int(valid_count); divisor = max(1, count)
        metrics[name] = {"valid_samples": count, "component_mse": float((rr+ii)/(2*divisor)),
                         "real_mse": float(rr/divisor), "imag_mse": float(ii/divisor),
                         "bias_real": float(r/divisor), "bias_imag": float(i/divisor),
                         "complex_bias_abs": math.hypot(r/divisor, i/divisor),
                         "max_complex_error": float(maximum), "max_baseline_component_mse": float(max_row),
                         "p99_error": percentiles[p]}

    def covariance(sums, counts, shape):
        s = np.frombuffer(sums, dtype="<f8").reshape(*shape, 6)
        n = np.frombuffer(counts, dtype="<u8").reshape(shape)
        denominator = np.maximum(1, n)
        x = s[..., 2] + 1j*s[..., 3]; y = s[..., 4] + 1j*s[..., 5]
        return (s[..., 0]+1j*s[..., 1])/denominator - x*np.conj(y)/denominator**2, n

    adjacent, adjacent_counts = covariance(raw[5], raw[6], (max(0, channels-1), products))
    cross, cross_counts = covariance(raw[7], raw[8], (products, products))
    channel_sum = np.frombuffer(raw[3], dtype="<f8").reshape(channels, products)
    channel_count = np.frombuffer(raw[4], dtype="<u8").reshape(channels, products)
    arrays = {"baseline_component_mse": np.frombuffer(raw[1], dtype="<f8").reshape(rows, products),
              "baseline_max_complex_error": np.frombuffer(raw[2], dtype="<f8").reshape(rows, products),
              "channel_component_mse": channel_sum/(2*np.maximum(1, channel_count)),
              "cross_product_error_covariance": cross, "cross_product_paired_counts": cross_counts,
              "adjacent_channel_error_covariance": adjacent, "adjacent_channel_paired_counts": adjacent_counts}
    return {"by_product": metrics, "protected_sample_count": raw[10],
            "protected_sample_bits": "exact", "metric_domain": "unflagged finite input samples",
            "covariance_definition": "E[e_x conj(e_y)] - E[e_x] conj(E[e_y]), restricted to paired valid samples",
            "science_validation": "not performed", "calculation": "native-moments-v2"}, arrays


def _numpy_quality(original, decoded, flags, names):
    rows, channels, products = original.shape
    valid = ~flags & np.isfinite(original)
    protected = protected_mask(original, flags)
    if original[protected].tobytes() != decoded[protected].tobytes():
        raise ValueError("Flagged/nonfinite/signed-zero sample bits changed")
    metrics = {}
    baseline_mse = np.zeros((rows, products), dtype=np.float64)
    baseline_max = np.zeros_like(baseline_mse)
    channel_mse = np.zeros((channels, products), dtype=np.float64)
    adjacent_covariance = np.zeros((max(0, channels - 1), products), dtype=np.complex128)
    adjacent_counts = np.zeros(adjacent_covariance.shape, dtype=np.uint64)
    for p, name in enumerate(names):
        mask = valid[:, :, p]
        metrics[name] = error_metrics(original[:, :, p], decoded[:, :, p], mask)
        error = np.zeros((rows, channels), dtype=np.complex128)
        error[mask] = decoded[:, :, p][mask].astype(np.complex128) - original[:, :, p][mask].astype(np.complex128)
        squared = error.real**2 + error.imag**2
        baseline_mse[:, p] = squared.sum(axis=1) / (2 * np.maximum(1, mask.sum(axis=1)))
        baseline_max[:, p] = np.max(np.abs(error), axis=1)
        channel_mse[:, p] = squared.sum(axis=0) / (2 * np.maximum(1, mask.sum(axis=0)))
        if channels > 1:
            paired = mask[:, :-1] & mask[:, 1:]
            counts = paired.sum(axis=0)
            x = np.where(paired, error[:, :-1], 0)
            y = np.where(paired, error[:, 1:], 0)
            adjacent_covariance[:, p] = np.sum(x * np.conj(y), axis=0) / np.maximum(1, counts) - x.sum(axis=0) * np.conj(y.sum(axis=0)) / np.maximum(1, counts)**2
            adjacent_counts[:, p] = counts
    # Paired, centered covariance, evaluated in small row slabs to bound RAM.
    cross_sum = np.zeros((products, products), dtype=np.complex128)
    x_sum = np.zeros_like(cross_sum)
    y_sum = np.zeros_like(cross_sum)
    cross_counts = np.zeros((products, products), dtype=np.uint64)
    for start in range(0, rows, 1024):
        stop = min(rows, start + 1024)
        mask = valid[start:stop]
        error = np.zeros(mask.shape, dtype=np.complex128)
        error[mask] = decoded[start:stop][mask].astype(np.complex128) - original[start:stop][mask].astype(np.complex128)
        for p in range(products):
            for q in range(products):
                paired = mask[:, :, p] & mask[:, :, q]
                x, y = error[:, :, p][paired], error[:, :, q][paired]
                cross_sum[p, q] += np.sum(x * np.conj(y))
                x_sum[p, q] += np.sum(x)
                y_sum[p, q] += np.sum(y)
                cross_counts[p, q] += len(x)
    count = np.maximum(1, cross_counts)
    covariance = cross_sum / count - x_sum * np.conj(y_sum) / count**2
    arrays = {"baseline_component_mse": baseline_mse, "baseline_max_complex_error": baseline_max,
              "channel_component_mse": channel_mse,
              "cross_product_error_covariance": covariance,
              "cross_product_paired_counts": cross_counts,
              "adjacent_channel_error_covariance": adjacent_covariance,
              "adjacent_channel_paired_counts": adjacent_counts}
    return {"by_product": metrics, "protected_sample_count": int(protected.sum()),
            "protected_sample_bits": "exact", "metric_domain": "unflagged finite input samples",
            "covariance_definition": "E[e_x conj(e_y)] - E[e_x] conj(E[e_y]), restricted to paired valid samples",
            "science_validation": "not performed"}, arrays


def _encoded_blocks(values, flags, permutation, boundaries, product, column, name, config,
                    references=None, autocorrelation=None):
    def jobs():
        for g, (start, stop) in enumerate(zip(boundaries[:-1], boundaries[1:])):
            rows = permutation[start:stop]
            for segment, lo in enumerate(range(0, values.shape[1], config.channel_segment)):
                yield rows, lo, min(values.shape[1], lo + config.channel_segment), g, segment

    def run(job):
        rows, lo, hi, g, segment = job
        stamp = time.perf_counter()
        block, flagged = values[rows, lo:hi, product], flags[rows, lo:hi, product]
        gathering = time.perf_counter() - stamp
        global_std = None
        if config.normalized_bounds is not None:
            if np.any(autocorrelation[rows] != autocorrelation[rows[0]]):
                raise ValueError('Normalized blocks cannot mix auto/cross classes')
            label = 'auto' if autocorrelation[rows[0]] else 'cross'
            global_std = references[name][label]['component_std']
        chosen = encode_block(block, flagged, config, config.bounds(name),
                              stream=(config.columns.index(column), product, g, segment), global_std=global_std)
        chosen.timings["memory_gathering"] = gathering
        return rows, lo, hi, chosen

    if config.workers == 1:
        yield from map(run, jobs())
    else:
        # Bound outstanding work; preserve submission order and deterministic
        # stream IDs. Only the parent thread ever touches HDF5/casacore.
        with ThreadPoolExecutor(max_workers=config.workers) as pool:
            iterator, pending = iter(jobs()), deque()
            for job in iterator:
                pending.append(pool.submit(run, job))
                if len(pending) >= 2 * config.workers:
                    yield pending.popleft().result()
            while pending:
                yield pending.popleft().result()


def write_visibility(writer, t, column, prefix, config, inspection, *, distribution_references=None):
    timings = Counter()
    stamp = time.perf_counter()
    values = np.asarray(t.getcol(column))
    flags = _flags(t)
    pairs = np.column_stack([t.getcol("ANTENNA1"), t.getcol("ANTENNA2")])
    uvw = np.asarray(t.getcol("UVW"))
    timings["input_reading"] += time.perf_counter() - stamp
    if flags.shape != values.shape:
        raise ValueError("FLAG and selected visibility shape mismatch")
    stamp = time.perf_counter()
    orders, statistics, counts = all_orders(values, flags, pairs, uvw, config, timings)
    rows, channels, products = values.shape
    names = inspection["layout"]["configurations"][0]["correlation_names"]
    references = None
    if config.normalized_bounds is not None:
        from .references import column_references
        references = column_references(t,column,distribution_references)
        if set(references) != set(names):
            raise ValueError('Distribution reference correlation identities differ')
    reconstructed = np.empty_like(values)
    spec = {"encoding": "quantized", "row_count": rows, "cell_shape": [channels, products],
            "dtype": "<c8", "codec": config.codebook, "channel_segment": config.channel_segment,
            "block_header": config.block_header,
            "product_names": names, "products": [],
            "statistics": writer.put(prefix + "/amplitude_medians", statistics, "quantization_statistics"),
            "valid_counts": writer.put(prefix + "/valid_counts", counts, "quantization_statistics")}
    shared_permutation_path = None
    if references is not None:
        spec['distribution_references'] = references
        from .references import key, table_context
        spec['distribution_reference_key'] = key(table_context(t))
    worker_timings = Counter()
    for p, order in enumerate(orders):
        product_prefix = prefix + f"/p{p}"
        permutation = order["permutation"]
        permutation_path = None
        if not np.array_equal(permutation, np.arange(rows, dtype=np.uint32)):
            if config.shared_order and shared_permutation_path is not None:
                permutation_path = shared_permutation_path
            else:
                permutation_path = writer.put(product_prefix + "/permutation", permutation, "permutations")
                if config.shared_order:
                    shared_permutation_path = permutation_path
        boundaries = order["group_boundaries"]
        boundaries_path = writer.put(product_prefix + "/group_boundaries", boundaries, "visibility_indexes")
        payload_path = product_prefix + "/payload"
        payload = writer.create_stream(payload_path, "visibility_block_payload")
        offsets, crc = [0], []
        histogram, sample_histogram, accounting = Counter(), Counter(), Counter()
        gamma_histogram = Counter()
        buffered, buffered_size, written = [], 0, 0
        pipeline_start = time.perf_counter()
        previous_writes = timings["payload_writing"]
        for source_rows, lo, hi, chosen in _encoded_blocks(values, flags, permutation, boundaries, p, column, names[p], config,
                                                         references, pairs[:,0] == pairs[:,1]):
            worker_timings.update(chosen.timings)
            if config.workers == 1:
                timings.update(chosen.timings)
            reconstructed[source_rows, lo:hi, p] = chosen.decoded
            offsets.append(offsets[-1] + len(chosen.blob))
            crc.append(zlib.crc32(chosen.blob))
            label = "float32" if chosen.bits == 0 else str(chosen.bits)
            histogram[label] += 1
            sample_histogram[label] += chosen.decoded.size
            accounting.update(chosen.accounting)
            if "gamma" in chosen.parameters:
                gamma_histogram[str(chosen.parameters["gamma"])] += 1
            buffered.append(chosen.blob)
            buffered_size += len(chosen.blob)
            if buffered_size >= 1024 * 1024:
                stamp = time.perf_counter()
                payload.resize((written + buffered_size,))
                payload[written:written + buffered_size] = np.frombuffer(b"".join(buffered), dtype=np.uint8)
                written += buffered_size
                buffered, buffered_size = [], 0
                timings["payload_writing"] += time.perf_counter() - stamp
        if config.workers > 1:
            timings["parallel_block_pipeline_wall"] += time.perf_counter() - pipeline_start - (timings["payload_writing"] - previous_writes)
        if buffered:
            stamp = time.perf_counter()
            payload.resize((written + buffered_size,))
            payload[written:written + buffered_size] = np.frombuffer(b"".join(buffered), dtype=np.uint8)
            timings["payload_writing"] += time.perf_counter() - stamp
        offsets_path = writer.put(product_prefix + "/offsets", np.asarray(offsets, dtype="<u8"), "visibility_indexes", compress=False)
        crc_path = writer.put(product_prefix + "/crc32", np.asarray(crc, dtype="<u4"), "visibility_indexes", compress=False)
        spec["products"].append({"permutation": permutation_path, "group_boundaries": boundaries_path,
                                 "payload": payload_path, "offsets": offsets_path, "crc32": crc_path,
                                 "block_count": len(crc), "component_bit_width_histogram": dict(histogram),
                                 "complex_sample_bit_width_histogram": dict(sample_histogram),
                                 "logical_bytes_before_secondary_compression": dict(accounting),
                                 "serialized_block_bytes": offsets[-1], "index_bytes": len(offsets) * 8 + len(crc) * 4,
                                 "bucket_reference": order["bucket_reference"], "buckets": order["buckets"],
                                 "grouping": order["grouping"],
                                 "permutation_direction": order["permutation_direction"]})
        if config.codebook == "power-midrise-v1":
            spec["products"][-1]["power_gamma_block_histogram"] = dict(gamma_histogram)
    stamp = time.perf_counter()
    summary, arrays = quality(values, reconstructed, flags, names, workers=config.workers)
    if config.keep_calib_only:
        summary["protected_sample_bits"] = "exact after replacing source-flagged values with complex positive zero; original flagged values discarded"
    for name, metrics in summary["by_product"].items():
        if not within_bounds(metrics, config.bounds(name)):
            raise ValueError(f"Whole-product error criteria failed for {name}; no complete archive was published")
    spec["quality_arrays"] = {name: writer.put(prefix + "/quality/" + name, array, "visibility_quality") for name, array in arrays.items()}
    spec["quality_summary"] = summary
    timings["quality_accounting"] += time.perf_counter() - stamp
    spec["timings_seconds"] = dict(timings)
    spec["worker_timings_seconds"] = dict(worker_timings)
    spec["worker_count"] = config.workers
    spec["worker_timer_scope"] = "*_cpu fields sum thread CPU seconds. Other fields sum elapsed task subregions, overlapping across workers; not additive to pipeline wall time. codec_total/native_candidate_search are nested."
    spec["byte_accounting_note"] = "Header/body serialization and index bytes are exact. Code/exception byte counts are before joint secondary compression; their compressed physical bytes cannot be attributed separately. HDF5 allocation and all metadata are counted in archive totals."
    return spec


def _product_layout(reader, spec, p):
    if spec.get("codec") not in CODEBOOK_IDS or not 0 <= p < len(spec["products"]):
        raise ValueError("Unsupported codec/product")
    product = spec["products"][p]
    rows = spec["row_count"]
    boundaries = reader.h5[product["group_boundaries"]][...]
    if len(boundaries) < 2 or boundaries[0] != 0 or boundaries[-1] != rows or np.any(boundaries[1:] <= boundaries[:-1]):
        raise ValueError("Invalid group boundaries")
    permutation = np.arange(rows, dtype=np.uint32) if product["permutation"] is None else reader.h5[product["permutation"]][...]
    if permutation.shape != (rows,) or not np.array_equal(np.sort(permutation), np.arange(rows)):
        raise ValueError("Invalid baseline permutation")
    segments = math.ceil(spec["cell_shape"][0] / spec["channel_segment"])
    count = (len(boundaries) - 1) * segments
    offsets = reader.h5[product["offsets"]][...]
    checksums = reader.h5[product["crc32"]][...]
    if not hasattr(reader, "block_datasets"):
        reader.block_datasets = {}
    if product["payload"] not in reader.block_datasets:
        # Keep the dataset handle/chunk cache alive across independent block reads.
        reader.block_datasets[product["payload"]] = reader.h5[product["payload"]]
    payload = reader.block_datasets[product["payload"]]
    if len(offsets) != count + 1 or len(checksums) != count or offsets[0] != 0 or offsets[-1] != payload.size or np.any(offsets[1:] <= offsets[:-1]):
        raise ValueError("Invalid block index")
    return product, boundaries, permutation, segments, offsets, payload, checksums


def read_visibility_block(reader, spec, product_index, block_index, layout=None):
    product, boundaries, permutation, segments, offsets, payload, checksums = layout or _product_layout(reader, spec, product_index)
    if not 0 <= block_index < len(checksums):
        raise ValueError("Block index out of range")
    group, segment = divmod(block_index, segments)
    rows = permutation[boundaries[group]:boundaries[group + 1]]
    lo = segment * spec["channel_segment"]
    hi = min(spec["cell_shape"][0], lo + spec["channel_segment"])
    blob = payload[int(offsets[block_index]):int(offsets[block_index + 1])].tobytes()
    if zlib.crc32(blob) != checksums[block_index]:
        raise ValueError("Visibility block CRC32 mismatch")
    header = spec.get("block_header", "literal-bits-v1")
    if header not in ("literal-bits-v1", "precision-nibble-v1") or not blob:
        raise ValueError("Unsupported visibility block header")
    if (blob[0] & 0xF0 == 0xA0) != (header == "precision-nibble-v1"):
        raise ValueError("Block precision header disagrees with manifest")
    return rows, lo, hi, decode_block(blob, (len(rows), hi - lo), spec["codec"])


def read_visibility(reader, spec):
    if not hasattr(reader, "visibility_cache"):
        reader.visibility_cache = {}
    key = spec["products"][0]["payload"]
    if key in reader.visibility_cache:
        return reader.visibility_cache[key]
    shape = (spec["row_count"], *spec["cell_shape"])
    if math.prod(shape) * 8 > reader.manifest["config"]["max_visibility_bytes"]:
        raise ValueError("Decoded visibility exceeds configured prototype RAM limit")
    values = np.empty(shape, dtype="<c8")
    for p in range(shape[2]):
        layout = _product_layout(reader, spec, p)
        for index in range(len(layout[-1])):
            rows, lo, hi, block = read_visibility_block(reader, spec, p, index, layout)
            values[rows, lo:hi, p] = block
    reader.visibility_cache[key] = values
    return values


def verify_visibility(reader, spec, t, column, configuration):
    config = Config(**configuration)
    config.validate()
    decoded = read_visibility(reader, spec)
    original, flags = np.asarray(t.getcol(column)), _flags(t)
    autocorrelation = np.asarray(t.getcol('ANTENNA1')) == np.asarray(t.getcol('ANTENNA2'))
    references = None
    if config.normalized_bounds is not None:
        from .references import column_references
        references = column_references(t,column,reader.manifest.get('distribution_references'))
        if references != spec.get('distribution_references'):
            raise ValueError('Column distribution references differ from archive inventory')
    method = spec["quality_summary"].get("calculation", "numpy-v1")
    summary, arrays = quality(original, decoded, flags, spec["product_names"], method, config.workers)
    if config.keep_calib_only:
        summary["protected_sample_bits"] = "exact after replacing source-flagged values with complex positive zero; original flagged values discarded"
    numerical_fallback = method == "native-moments-v2" and _native is None
    if numerical_fallback:
        # The portable verification path has a different floating reduction
        # order. Counts/protected values stay exact; numerical moments are
        # compared with a documented rounding tolerance. Bounds stay strict.
        assert_quality_close(spec["quality_summary"], summary, spec["quality_arrays"], arrays, reader)
    elif summary != spec["quality_summary"]:
        raise ValueError("Reconstructed errors differ from encode-time statistics")
    for name, path in spec["quality_arrays"].items():
        if not numerical_fallback and not equal(arrays[name], reader.h5[path][...]):
            raise ValueError(f"Reconstructed error accounting differs: {name}")
    for p, name in enumerate(spec["product_names"]):
        layout = _product_layout(reader, spec, p)
        for index in range(len(layout[-1])):
            rows, lo, hi, block = read_visibility_block(reader, spec, p, index, layout)
            valid = ~flags[rows, lo:hi, p] & np.isfinite(original[rows, lo:hi, p])
            measured = error_metrics(original[rows, lo:hi, p], block, valid)
            global_std = None
            if references is not None:
                if np.any(autocorrelation[rows] != autocorrelation[rows[0]]):
                    raise ValueError('Normalized block mixes auto/cross classes')
                label = 'auto' if autocorrelation[rows[0]] else 'cross'
                global_std = references[name][label]['component_std']
            bounds = configured_bounds(original[rows, lo:hi, p], flags[rows, lo:hi, p],
                                       config, config.bounds(name), global_std)
            if not within_bounds(measured, bounds):
                raise ValueError(f"Block failed measured error budget: product {name}, block {index}")
        if not within_bounds(summary["by_product"][name], config.bounds(name)):
            raise ValueError(f"Whole-product measured error criteria failed: {name}")
    return summary


def assert_quality_close(expected, actual, paths, arrays, reader):
    if expected["protected_sample_count"] != actual["protected_sample_count"]:
        raise ValueError("Protected sample counts differ")
    scale = max((v["component_mse"] for v in actual["by_product"].values()), default=0)
    for product, metrics in actual["by_product"].items():
        for name, value in metrics.items():
            reference = expected["by_product"][product][name]
            atol = 1e-12 * (scale if "mse" in name else math.sqrt(scale))
            if (name == "valid_samples" and value != reference) or not math.isclose(value, reference, rel_tol=1e-10, abs_tol=atol):
                raise ValueError(f"Reconstructed error statistic differs: {product}/{name}")
    for name, path in paths.items():
        reference = reader.h5[path][...]
        array = arrays[name]
        atol = 1e-12 * (math.sqrt(scale) if "max_complex_error" in name else scale)
        if (array.dtype.kind in "iu" and not equal(array, reference)) or not np.allclose(array, reference, rtol=1e-10, atol=atol):
            raise ValueError(f"Reconstructed error accounting differs: {name}")
