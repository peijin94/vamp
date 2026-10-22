"""Stable baseline permutations; spectra and ordered antenna IDs remain intact."""

import math
import time

import numpy as np


def amplitude_statistics(values, flags):
    nrow, _, products = values.shape
    statistics = np.zeros((nrow, products), dtype=np.float64)
    counts = np.zeros((nrow, products), dtype=np.uint32)
    for product in range(products):
        data = values[:, :, product]
        valid = ~flags[:, :, product] & np.isfinite(data)
        counts[:, product] = valid.sum(axis=1)
        amplitude = np.hypot(data.real.astype(np.float64), data.imag.astype(np.float64))
        amplitude[~valid] = np.nan
        absent = counts[:, product] == 0
        amplitude[absent, 0] = 0.0
        statistics[:, product] = np.nanmedian(amplitude, axis=1)
        statistics[absent, product] = 0.0
    return statistics, counts


def make_order(statistic, valid_counts, autocorrelation, lengths, config, *, uvw=None, pairs=None):
    rows = len(statistic)
    original = np.arange(rows, dtype=np.uint32)
    no_valid = valid_counts == 0
    bucket = np.zeros(rows, dtype=np.int32)
    zero = statistic == 0
    reference = None
    extra_labels = []
    grouping = {"method": config.ordering,
                "classes": "autocorrelations first, then cross-correlations; no-valid samples separate"}
    if config.ordering == "buckets":
        positive = ~no_valid & ~zero
        reference = config.bucket_reference or (float(np.median(statistic[positive])) if np.any(positive) else 1.0)
        bucket[positive] = np.floor(np.log2(statistic[positive]) - math.log2(reference)).astype(np.int32)
    if config.ordering == "original":
        permutation = original
    elif config.ordering == "class-original":
        permutation = np.lexsort((original, no_valid, ~autocorrelation)).astype(np.uint32)
    elif config.ordering == "antenna":
        if pairs is None or np.shape(pairs) != (rows, 2):
            raise ValueError("Antenna grouping requires ordered ANTENNA1/ANTENNA2 identities")
        endpoint = np.asarray(pairs)[:, config.antenna_endpoint - 1]
        permutation = np.lexsort((original, endpoint, no_valid, ~autocorrelation)).astype(np.uint32)
        extra_labels.append(endpoint)
        grouping.update(endpoint_column=f"ANTENNA{config.antenna_endpoint}", within_endpoint="original row order")
    elif config.ordering in ("orientation", "length-orientation"):
        if uvw is None or np.shape(uvw) != (rows, 3):
            raise ValueError("Orientation grouping requires original UVW in metres")
        uvw = np.asarray(uvw)
        invalid = ~np.all(np.isfinite(uvw), axis=1) | ~np.isfinite(lengths)
        zero_projected = (uvw[:, 0] == 0) & (uvw[:, 1] == 0)
        special = np.where(invalid, 2, np.where(zero_projected, 1, 0))
        angle = np.zeros(rows, dtype=np.float64)
        ordinary = special == 0
        angle[ordinary] = np.mod(np.arctan2(uvw[ordinary, 1], uvw[ordinary, 0]), 2 * np.pi)
        # Ordered UVW bearing, not geographic antenna bearing. No pi folding,
        # antenna swapping, conjugation, or changes to the archived coordinates.
        angular_bin = np.minimum((angle * config.orientation_bins / (2 * np.pi)).astype(np.int32), config.orientation_bins - 1)
        extra_labels.extend([special, angular_bin])
        if config.ordering == "length-orientation":
            radial_bin = np.zeros(rows, dtype=np.int32)
            positive = ~invalid & (lengths > 0)
            radial_bin[positive] = np.floor(np.log2(lengths[positive]) - math.log2(config.length_reference_m)).astype(np.int32)
            extra_labels.append(radial_bin)
            permutation = np.lexsort((original, angular_bin, radial_bin, special, no_valid, ~autocorrelation)).astype(np.uint32)
            grouping.update(length_reference_m=config.length_reference_m, radial_bins="floor(log2(norm(UVW)/length_reference_m)); factor of two")
        else:
            permutation = np.lexsort((original, angle, special, no_valid, ~autocorrelation)).astype(np.uint32)
        grouping.update(orientation="atan2(V,U) modulo 2*pi, original ordered metre UVW", orientation_bins=config.orientation_bins,
                        bin_width_degrees=360.0 / config.orientation_bins,
                        special_geometry="zero projected UV and nonfinite geometry form separate groups")
    else:
        score = statistic
        if config.ordering == "length":
            score = np.where(np.isfinite(lengths), lengths, np.inf)
            extra_labels.append(~np.isfinite(lengths))
            grouping["length"] = "hypot(hypot(U,V),W), original metre UVW"
        if config.ordering == "buckets":
            permutation = np.lexsort((original, bucket, ~zero, no_valid, ~autocorrelation)).astype(np.uint32)
        else:
            permutation = np.lexsort((original, score, no_valid, ~autocorrelation)).astype(np.uint32)
    # Labels stop blocks from crossing auto/cross, no-valid-data or bucket edges.
    labels = np.column_stack([autocorrelation, no_valid])
    if extra_labels:
        labels = np.column_stack([labels, *extra_labels])
    if config.ordering == "buckets":
        labels = np.column_stack([labels, zero, bucket])
    ordered_labels = labels[permutation]
    breaks = np.r_[0, np.flatnonzero(np.any(ordered_labels[1:] != ordered_labels[:-1], axis=1)) + 1, rows]
    boundaries = [0]
    bucket_records = []
    for start, stop in zip(breaks[:-1], breaks[1:]):
        boundaries.extend(range(int(start) + config.baseline_group, int(stop), config.baseline_group))
        boundaries.append(int(stop))
        if config.ordering == "buckets":
            row = permutation[start]
            kind = "no-valid-samples" if no_valid[row] else "zero" if zero[row] else "positive"
            item = {"start": int(start), "stop": int(stop), "class": "auto" if autocorrelation[row] else "cross", "kind": kind}
            if kind == "positive":
                item.update(bucket_index=int(bucket[row]), lower=math.ldexp(reference, int(bucket[row])), upper=math.ldexp(reference, int(bucket[row] + 1)))
            bucket_records.append(item)
    sizes, number = np.unique(np.diff(boundaries), return_counts=True)
    grouping["group_size_histogram"] = {str(int(size)): int(n) for size, n in zip(sizes, number)}
    grouping["group_count"] = len(boundaries) - 1
    grouping["partial_group_count"] = int(np.count_nonzero(np.diff(boundaries) < config.baseline_group))
    return {"permutation": permutation, "group_boundaries": np.asarray(boundaries, dtype="<u4"), "grouping": grouping,
            "bucket_reference": reference, "buckets": bucket_records,
            "permutation_direction": "perm[j] is original row index at archived baseline position j; restored[perm,:]=decoded"}


def all_orders(values, flags, pairs, uvw, config, timings=None):
    stamp = time.perf_counter()
    statistics, counts = amplitude_statistics(values, flags)
    if timings is not None:
        timings["amplitude_statistics"] += time.perf_counter() - stamp
    stamp = time.perf_counter()
    auto = pairs[:, 0] == pairs[:, 1]
    lengths = np.hypot(np.hypot(uvw[:, 0], uvw[:, 1]), uvw[:, 2])
    shared = None
    if config.shared_order:
        masked = np.where(counts > 0, statistics, np.nan)
        absent = np.all(counts == 0, axis=1)
        masked[absent, 0] = 0
        shared = np.nanmedian(masked, axis=1)
        shared[absent] = 0
    orders = []
    for product in range(values.shape[2]):
        score = shared if shared is not None else statistics[:, product]
        valid = counts.sum(axis=1) if shared is not None else counts[:, product]
        orders.append(make_order(score, valid, auto, lengths, config, uvw=uvw, pairs=pairs))
    if timings is not None:
        timings["baseline_sorting_and_grouping"] += time.perf_counter() - stamp
    return orders, statistics, counts
