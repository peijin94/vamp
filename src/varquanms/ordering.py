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


def make_order(statistic, valid_counts, autocorrelation, lengths, config):
    rows = len(statistic)
    original = np.arange(rows, dtype=np.uint32)
    no_valid = valid_counts == 0
    bucket = np.zeros(rows, dtype=np.int32)
    zero = statistic == 0
    reference = None
    if config.ordering == "buckets":
        positive = ~no_valid & ~zero
        reference = config.bucket_reference or (float(np.median(statistic[positive])) if np.any(positive) else 1.0)
        bucket[positive] = np.floor(np.log2(statistic[positive]) - math.log2(reference)).astype(np.int32)
    if config.ordering == "original":
        permutation = original
    else:
        score = statistic
        if config.ordering == "length":
            score = np.where(np.isfinite(lengths), lengths, np.inf)
        if config.ordering == "buckets":
            permutation = np.lexsort((original, bucket, ~zero, no_valid, ~autocorrelation)).astype(np.uint32)
        else:
            permutation = np.lexsort((original, score, no_valid, ~autocorrelation)).astype(np.uint32)
    # Labels stop blocks from crossing auto/cross, no-valid-data or bucket edges.
    labels = np.column_stack([autocorrelation, no_valid])
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
    return {"permutation": permutation, "group_boundaries": np.asarray(boundaries, dtype="<u4"),
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
        orders.append(make_order(score, valid, auto, lengths, config))
    if timings is not None:
        timings["baseline_sorting_and_grouping"] += time.perf_counter() - stamp
    return orders, statistics, counts
