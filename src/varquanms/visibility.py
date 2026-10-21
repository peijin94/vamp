"""One-product quantization blocks, invertible orderings and measured quality."""

from collections import Counter
import math
import time
import zlib

import numpy as np

from .codec import decode_block, encode_block, error_metrics, within_bounds
from .config import Config
from .metadata import equal
from .ordering import all_orders


def _flags(t):
    flags = np.asarray(t.getcol("FLAG"), dtype=bool)
    flags |= np.asarray(t.getcol("FLAG_ROW"), dtype=bool)[:, None, None]
    return flags


def protected_mask(values, flags):
    return flags | ~np.isfinite(values) | ((values.real == 0) & np.signbit(values.real)) | ((values.imag == 0) & np.signbit(values.imag))


def quality(original, decoded, flags, names):
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


def write_visibility(writer, t, column, prefix, config, inspection):
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
    reconstructed = np.empty_like(values)
    spec = {"encoding": "quantized", "row_count": rows, "cell_shape": [channels, products],
            "dtype": "<c8", "codec": "uniform-symmetric-v1", "channel_segment": config.channel_segment,
            "product_names": names, "products": [],
            "statistics": writer.put(prefix + "/amplitude_medians", statistics, "quantization_statistics"),
            "valid_counts": writer.put(prefix + "/valid_counts", counts, "quantization_statistics")}
    shared_permutation_path = None
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
        buffered, buffered_size, written = [], 0, 0
        for g, (start, stop) in enumerate(zip(boundaries[:-1], boundaries[1:])):
            source_rows = permutation[start:stop]
            for segment, lo in enumerate(range(0, channels, config.channel_segment)):
                hi = min(channels, lo + config.channel_segment)
                stamp = time.perf_counter()
                block = values[source_rows, lo:hi, p]
                flagged = flags[source_rows, lo:hi, p]
                timings["memory_gathering"] += time.perf_counter() - stamp
                stamp = time.perf_counter()
                chosen = encode_block(block, flagged, config, config.bounds(names[p]),
                                      stream=(config.columns.index(column), p, g, segment))
                timings.update(chosen.timings)
                reconstructed[source_rows, lo:hi, p] = chosen.decoded
                offsets.append(offsets[-1] + len(chosen.blob))
                crc.append(zlib.crc32(chosen.blob))
                label = "float32" if chosen.bits == 0 else str(chosen.bits)
                histogram[label] += 1
                sample_histogram[label] += block.size
                accounting.update(chosen.accounting)
                buffered.append(chosen.blob)
                buffered_size += len(chosen.blob)
                if buffered_size >= 1024 * 1024:
                    stamp = time.perf_counter()
                    payload.resize((written + buffered_size,))
                    payload[written:written + buffered_size] = np.frombuffer(b"".join(buffered), dtype=np.uint8)
                    written += buffered_size
                    buffered, buffered_size = [], 0
                    timings["payload_writing"] += time.perf_counter() - stamp
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
                                 "permutation_direction": order["permutation_direction"]})
    stamp = time.perf_counter()
    summary, arrays = quality(values, reconstructed, flags, names)
    for name, metrics in summary["by_product"].items():
        if not within_bounds(metrics, config.bounds(name)):
            raise ValueError(f"Whole-product error criteria failed for {name}; no complete archive was published")
    spec["quality_arrays"] = {name: writer.put(prefix + "/quality/" + name, array, "visibility_quality") for name, array in arrays.items()}
    spec["quality_summary"] = summary
    timings["quality_accounting"] += time.perf_counter() - stamp
    spec["timings_seconds"] = dict(timings)
    spec["byte_accounting_note"] = "Header/body serialization and index bytes are exact. Code/exception byte counts are before joint secondary compression; their compressed physical bytes cannot be attributed separately. HDF5 allocation and all metadata are counted in archive totals."
    return spec


def _product_layout(reader, spec, p):
    if spec.get("codec") != "uniform-symmetric-v1" or not 0 <= p < len(spec["products"]):
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
    return rows, lo, hi, decode_block(blob, (len(rows), hi - lo))


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
    summary, arrays = quality(original, decoded, flags, spec["product_names"])
    if summary != spec["quality_summary"]:
        raise ValueError("Reconstructed errors differ from encode-time statistics")
    for name, path in spec["quality_arrays"].items():
        if not equal(arrays[name], reader.h5[path][...]):
            raise ValueError(f"Reconstructed error accounting differs: {name}")
    for p, name in enumerate(spec["product_names"]):
        layout = _product_layout(reader, spec, p)
        for index in range(len(layout[-1])):
            rows, lo, hi, block = read_visibility_block(reader, spec, p, index, layout)
            valid = ~flags[rows, lo:hi, p] & np.isfinite(original[rows, lo:hi, p])
            measured = error_metrics(original[rows, lo:hi, p], block, valid)
            if not within_bounds(measured, config.bounds(name)):
                raise ValueError(f"Block failed measured error budget: product {name}, block {index}")
        if not within_bounds(summary["by_product"][name], config.bounds(name)):
            raise ValueError(f"Whole-product measured error criteria failed: {name}")
    return summary
