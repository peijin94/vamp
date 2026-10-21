"""Read every declared column without modifying the Measurement Set.

This module produces an inspection inventory, not an archival representation.
Numeric metadata arrays carry their dtype, shape and original bytes. MAIN
visibility arrays are scanned in bounded row batches, rather than embedded.
"""

from __future__ import annotations

import base64
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import time

import casacore
from casacore.tables import table
import numpy as np

from . import __version__

VISIBILITY_COLUMNS = ("DATA", "FLOAT_DATA", "CORRECTED_DATA", "MODEL_DATA")
FLAG_WEIGHT_COLUMNS = (
    "FLAG", "FLAG_ROW", "FLAG_CATEGORY", "WEIGHT", "SIGMA",
    "WEIGHT_SPECTRUM", "SIGMA_SPECTRUM",
)
ROW_KEYS = (
    "ANTENNA1", "ANTENNA2", "FEED1", "FEED2", "TIME", "TIME_CENTROID",
    "DATA_DESC_ID", "FIELD_ID", "SCAN_NUMBER", "ARRAY_ID", "OBSERVATION_ID",
    "PROCESSOR_ID", "STATE_ID", "INTERVAL", "EXPOSURE", "TIME_EXTRA_PREC",
    "ANTENNA3", "FEED3", "PULSAR_BIN", "PULSAR_GATE_ID", "PHASE_ID", "BASELINE_REF",
)
CORRELATIONS = {
    1: "I", 2: "Q", 3: "U", 4: "V", 5: "RR", 6: "RL", 7: "LR", 8: "LL",
    9: "XX", 10: "XY", 11: "YX", 12: "YY",
}


def json_value(value):
    """Represent numeric metadata without discarding NaN payloads or -0 bits."""
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, np.ndarray):
        if value.dtype.kind in "biufc":
            a = np.ascontiguousarray(value)
            result = {"dtype": a.dtype.str, "shape": list(value.shape),
                      "data_base64": base64.b64encode(a.tobytes()).decode("ascii")}
            if a.size <= 1024:
                result["values"] = json_value(value.tolist())
            return result
        return {"dtype": value.dtype.str, "shape": list(value.shape),
                "values": json_value(value.tolist())}
    if isinstance(value, np.generic):
        return json_value(np.asarray(value))
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return {"float64_bits_hex": np.asarray(value, dtype="<f8").tobytes().hex()}
    if isinstance(value, complex):
        return {"complex128_bits_hex": np.asarray(value, dtype="<c16").tobytes().hex()}
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f"Unsupported metadata value: {type(value).__name__}")


def unique_bits(values):
    """Deduplicate complete records by bytes, preserving first-occurrence order."""
    a = np.asarray(values)
    if a.dtype.hasobject:
        raise TypeError("Object arrays have no portable bitwise record identity")
    if len(a) == 0:
        return a.copy(), np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    width = a.dtype.itemsize * math.prod(a.shape[1:])
    if width == 0:
        return a[:1].copy(), np.zeros(len(a), dtype=np.int64), np.asarray([len(a)], dtype=np.int64)
    keys = np.ascontiguousarray(a).reshape(len(a), -1).view(np.dtype((np.void, width))).reshape(-1)
    _, first, inverse, counts = np.unique(
        keys, return_index=True, return_inverse=True, return_counts=True)
    order = np.argsort(first)
    old_to_new = np.empty(len(order), dtype=np.int64)
    old_to_new[order] = np.arange(len(order))
    return a[first[order]], old_to_new[inverse], counts[order]


def _open(path):
    # A frozen local copy is expected. usernoread avoids lock-file modifications.
    return table(str(path), readonly=True, lockoptions="usernoread", ack=False)


def _batches(t, column, chunk_rows):
    for start in range(0, t.nrows(), chunk_rows):
        count = min(chunk_rows, t.nrows() - start)
        try:
            values = np.asarray(t.getcol(column, startrow=start, nrow=count))
        except RuntimeError:
            # Variable shapes and undefined cells cannot be read with getcol.
            # Reading each defined cell also distinguishes absent from empty.
            for row in range(start, start + count):
                if t.iscelldefined(column, row):
                    cell = np.asarray(t.getcell(column, row))
                    yield np.asarray([row]), np.expand_dims(cell, 0)
                else:
                    yield np.asarray([row]), None
        else:
            yield np.arange(start, start + count), values


def _logical_bytes(a):
    if a.dtype.kind in "USO":
        return sum(len(str(v).encode("utf-8")) for v in a.ravel())
    return a.nbytes


def numeric_counts(a):
    """Count components independently, including signed zero components."""
    if a.dtype.kind not in "fc":
        return {}
    components = (a.real, a.imag) if a.dtype.kind == "c" else (a,)
    return {
        "nonfinite_values": int(np.count_nonzero(~np.isfinite(a))),
        "nan_components": sum(int(np.count_nonzero(np.isnan(v))) for v in components),
        "positive_inf_components": sum(int(np.count_nonzero(np.isposinf(v))) for v in components),
        "negative_inf_components": sum(int(np.count_nonzero(np.isneginf(v))) for v in components),
        "negative_zero_components": sum(int(np.count_nonzero((v == 0) & np.signbit(v))) for v in components),
        "zero_values": int(np.count_nonzero(a == 0)),
    }


def _scan_column(t, column, chunk_rows):
    descriptor = t.getcoldesc(column)
    shapes = Counter()
    dtypes = set()
    stats = Counter()
    byte_count = 0
    digest = hashlib.sha256()
    undefined = 0
    for rows, values in _batches(t, column, chunk_rows):
        if values is None:
            undefined += len(rows)
            continue
        shapes[str(list(values.shape[1:]))] += len(rows)
        dtypes.add(values.dtype.str)
        byte_count += _logical_bytes(values)
        stats["elements"] += values.size
        if values.dtype.kind == "b":
            stats["true_values"] += int(np.count_nonzero(values))
        stats.update(numeric_counts(values))
        if values.dtype.kind in "USO":
            for value in values.ravel():
                encoded = str(value).encode("utf-8")
                digest.update(len(encoded).to_bytes(8, "little"))
                digest.update(encoded)
        else:
            digest.update(np.ascontiguousarray(values).tobytes())
    return {
        "descriptor": json_value(descriptor),
        "keywords": json_value(t.getcolkeywords(column)),
        "units": json_value(t.getcolkeywords(column).get("QuantumUnits")),
        "measure": json_value(t.getcolkeywords(column).get("MEASINFO")),
        "storage_manager_type": descriptor.get("dataManagerType"),
        "storage_manager_group": descriptor.get("dataManagerGroup"),
        "cell_shapes": dict(shapes), "numpy_dtypes": sorted(dtypes),
        "undefined_cells": undefined, "logical_array_bytes": int(byte_count),
        "statistics": dict(stats), "defined_value_bytes_sha256": digest.hexdigest(),
        "digest_note": "Defined cell values in row/C order; shapes and undefined counts are separate. Strings are UTF-8 with uint64 little-endian lengths.",
    }


def _scan_table(t, chunk_rows, include_rows=False):
    result = {
        "row_count": t.nrows(), "info": json_value(t.info()),
        "keywords": json_value(t.getkeywords()),
        "descriptor": json_value(t.getdesc()),
        "storage_managers": json_value(t.getdminfo()),
        "columns": {c: _scan_column(t, c, chunk_rows) for c in t.colnames()},
    }
    result["logical_array_bytes"] = sum(c["logical_array_bytes"] for c in result["columns"].values())
    if include_rows:
        result["rows"] = [
            {c: json_value(np.asarray(t.getcell(c, row))) if t.iscelldefined(c, row)
             else {"undefined": True} for c in t.colnames()}
            for row in range(t.nrows())
        ]
    return result


def _physical_inventory(root, managers):
    files = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise ValueError(f"Symlink encountered in MS inventory: {path}")
        if not path.is_file():
            continue
        stat = path.stat()
        relative = path.relative_to(root)
        table_name = str(relative.parent) if relative.parent != Path(".") else "MAIN"
        match = re.match(r"table\.f(\d+)(?:\D.*)?$", path.name)
        candidates = [m for m in managers.get(table_name, {}).values()
                      if match and int(m.get("SEQNR", -1)) == int(match.group(1))]
        files.append({
            "path": relative.as_posix(), "bytes": stat.st_size,
            "allocated_bytes": getattr(stat, "st_blocks", 0) * 512,
            "table": table_name,
            "candidate_storage_manager_columns": [str(c) for m in candidates for c in m.get("COLUMNS", [])],
            "attribution": "Shared storage file; no exact per-column attribution" if candidates else "Table overhead or unassigned file",
        })
    return {
        "files": files, "file_bytes": sum(f["bytes"] for f in files),
        "allocated_file_bytes": sum(f["allocated_bytes"] for f in files),
        "note": "File lengths and allocated blocks measured on the local filesystem, excluding directory inodes. Storage-manager file assignment is by table.f<SEQNR> naming; shared files cannot be allocated exactly to individual columns.",
    }


def _key_matrix(columns, names):
    arrays = [np.ascontiguousarray(columns[c]).view(np.uint8).reshape(len(columns[c]), -1)
              for c in names]
    return np.concatenate(arrays, axis=1)


def _equal_cells(t, row_a, row_b):
    for column in t.colnames():
        defined = t.iscelldefined(column, row_a), t.iscelldefined(column, row_b)
        if defined[0] != defined[1]:
            return False
        if not defined[0]:
            continue
        a, b = np.asarray(t.getcell(column, row_a)), np.asarray(t.getcell(column, row_b))
        if a.dtype != b.dtype or a.shape != b.shape:
            return False
        if a.dtype.kind in "USO":
            if not np.array_equal(a, b):
                return False
        elif a.tobytes() != b.tobytes():
            return False
    return True


def _row_inventory(t):
    names = [c for c in ROW_KEYS if c in t.colnames()]
    columns = {c: np.asarray(t.getcol(c)) for c in names}
    if any(v.ndim != 1 or v.dtype.kind not in "biuf" for v in columns.values()):
        raise ValueError("A required row-identity column is not a numeric scalar")
    keys = _key_matrix(columns, names)
    _, inverse, counts = unique_bits(keys)
    duplicate_groups = []
    exact_duplicates = 0
    for group in np.flatnonzero(counts > 1):
        rows = np.flatnonzero(inverse == group).tolist()
        duplicate_groups.append(rows)
        representatives = []
        for row in rows:
            if any(_equal_cells(t, other, row) for other in representatives):
                exact_duplicates += 1
            else:
                representatives.append(row)
    distinct = {}
    for name, values in columns.items():
        unique, _, frequencies = unique_bits(values)
        distinct[name] = {"count": len(unique), "values": json_value(unique),
                          "row_counts": frequencies.tolist()}
    pairs = np.column_stack([columns["ANTENNA1"], columns["ANTENNA2"]])
    unique_pairs, _, pair_counts = unique_bits(pairs)
    context_names = [c for c in names if c not in ("ANTENNA1", "ANTENNA2", "DATA_DESC_ID")]
    _, contexts, _ = unique_bits(_key_matrix(columns, context_names))
    return {
        "identity_columns": names, "distinct": distinct,
        "row_key_duplicate_count": int(np.sum(counts - 1)),
        "duplicate_row_groups": duplicate_groups,
        "exact_duplicate_rows_beyond_first": exact_duplicates,
        "ordered_antenna_pairs": json_value(unique_pairs),
        "ordered_pair_row_counts": pair_counts.tolist(),
        "autocorrelation_rows": int(np.count_nonzero(pairs[:, 0] == pairs[:, 1])),
        "antenna1_greater_than_antenna2_rows": int(np.count_nonzero(pairs[:, 0] > pairs[:, 1])),
        "context_count_excluding_pair_and_data_description": len(np.unique(contexts)),
    }, columns


def _visibility_inventory(t, selected, chunk_rows):
    result = {}
    for column in selected:
        totals = Counter()
        shape_groups = {}
        for rows, values in _batches(t, column, chunk_rows):
            if values is None:
                continue
            if values.ndim != 3:
                raise ValueError(f"{column} must have channel/correlation cells; got {values.shape[1:]}")
            flags = np.stack([t.getcell("FLAG", int(row)) for row in rows])
            if flags.shape != values.shape:
                raise ValueError(f"FLAG shape does not match {column} at row {rows[0]}")
            row_flags = np.asarray([t.getcell("FLAG_ROW", int(row)) for row in rows], dtype=bool)
            flagged = flags | row_flags[:, None, None]
            finite = np.isfinite(values)
            valid = ~flagged & finite
            totals["complex_or_float_samples"] += values.size
            totals["flagged_samples_including_flag_row"] += int(np.count_nonzero(flagged))
            totals["all_flagged_rows"] += int(np.count_nonzero(np.all(flagged, axis=(1, 2))))
            totals["nonfinite_unflagged_samples"] += int(np.count_nonzero(~finite & ~flagged))
            totals["nonfinite_flagged_samples"] += int(np.count_nonzero(~finite & flagged))
            key = str(list(values.shape[1:]))
            if key not in shape_groups:
                shape_groups[key] = {"row_count": 0,
                                     "valid_sample_count_by_channel_product": np.zeros(values.shape[1:], dtype=np.int64)}
            group = shape_groups[key]
            group["row_count"] += len(rows)
            group["valid_sample_count_by_channel_product"] += np.sum(valid, axis=0)
        for group in shape_groups.values():
            counts = group["valid_sample_count_by_channel_product"]
            group["all_invalid_channel_product_positions"] = np.argwhere(counts == 0).tolist()
            group["valid_sample_count_by_product"] = counts.sum(axis=0).tolist()
            group["valid_sample_count_by_channel_product"] = counts.tolist()
        result[column] = {**dict(totals), "by_cell_shape": shape_groups,
                          "note": "Flagged samples and nonfinite samples are excluded from valid counts; original values have not been changed."}
    return result


def _layout(t, rows, columns, subtable_paths, selected, main_inventory):
    issues = []
    required = ("ANTENNA", "DATA_DESCRIPTION", "SPECTRAL_WINDOW", "POLARIZATION", "FIELD")
    absent = [name for name in required if name not in subtable_paths]
    if absent:
        raise ValueError(f"Missing required subtables: {', '.join(absent)}")
    with _open(subtable_paths["ANTENNA"]) as antennas:
        antenna_count = antennas.nrows()
    with _open(subtable_paths["DATA_DESCRIPTION"]) as dd, _open(subtable_paths["SPECTRAL_WINDOW"]) as spw, _open(subtable_paths["POLARIZATION"]) as pol:
        configs = []
        for dd_id in np.unique(columns["DATA_DESC_ID"]):
            dd_id = int(dd_id)
            if not 0 <= dd_id < dd.nrows():
                raise ValueError(f"DATA_DESC_ID out of range: {dd_id}")
            spw_id = int(dd.getcell("SPECTRAL_WINDOW_ID", dd_id))
            pol_id = int(dd.getcell("POLARIZATION_ID", dd_id))
            if not 0 <= spw_id < spw.nrows() or not 0 <= pol_id < pol.nrows():
                raise ValueError(f"Invalid spectral-window/polarization reference for DATA_DESC_ID {dd_id}")
            nchan, ncorr = int(spw.getcell("NUM_CHAN", spw_id)), int(pol.getcell("NUM_CORR", pol_id))
            corr_types = np.asarray(pol.getcell("CORR_TYPE", pol_id))
            frequencies = np.asarray(spw.getcell("CHAN_FREQ", spw_id))
            if len(frequencies) != nchan or len(corr_types) != ncorr:
                raise ValueError("Spectral/polarization metadata lengths are inconsistent")
            row_indices = np.flatnonzero(columns["DATA_DESC_ID"] == dd_id)
            pairs = np.column_stack([columns["ANTENNA1"][row_indices], columns["ANTENNA2"][row_indices]])
            if np.any(pairs < 0) or np.any(pairs >= antenna_count):
                raise ValueError("Antenna reference is outside the ANTENNA table")
            unique_pairs, _, _ = unique_bits(pairs)
            upper = bool(np.all(pairs[:, 0] <= pairs[:, 1]))
            lower = bool(np.all(pairs[:, 0] >= pairs[:, 1]))
            orientation = "upper_triangle" if upper else "lower_triangle" if lower else "mixed_ordered_pairs"
            auto_pairs = int(np.count_nonzero(unique_pairs[:, 0] == unique_pairs[:, 1]))
            observed_cross = len(unique_pairs) - auto_pairs
            missing_cross = antenna_count * (antenna_count - 1) // 2 - observed_cross if upper or lower else None
            for name in selected:
                if t.getcoldesc(name).get("valueType") != "complex":
                    issues.append(f"{name}: initial quantization target requires complex64")
                for row in row_indices:
                    if not t.iscelldefined(name, int(row)):
                        issues.append(f"{name}: undefined visibility cells require explicit presence maps")
                        break
                    if np.asarray(t.getcell(name, int(row))).shape != (nchan, ncorr):
                        issues.append(f"{name}: cell shape differs from its DATA_DESCRIPTION at row {row}")
                        break
            configs.append({
                "data_description_id": dd_id, "spectral_window_id": spw_id,
                "polarization_id": pol_id, "row_count": len(row_indices),
                "channel_count": nchan, "correlation_count": ncorr,
                "correlation_types": corr_types.tolist(),
                "correlation_names": [CORRELATIONS.get(int(c), f"CORR_TYPE_{c}") for c in corr_types],
                "correlation_products": json_value(np.asarray(pol.getcell("CORR_PRODUCT", pol_id))),
                "channel_frequencies": json_value(frequencies),
                "channel_widths": json_value(np.asarray(spw.getcell("CHAN_WIDTH", spw_id))),
                "frequency_reference_code": int(spw.getcell("MEAS_FREQ_REF", spw_id)),
                "ordered_pair_count": len(unique_pairs), "antenna_table_rows": antenna_count,
                "pair_orientation": orientation,
                "missing_cross_pairs_relative_to_full_antenna_table": missing_cross,
                "missing_autocorrelation_pairs_relative_to_full_antenna_table": antenna_count - auto_pairs,
                "missing_pairs_note": "Counts use the full ANTENNA table and observed triangular orientation. Omitted pairs may be intentional; no absent records are replaced by flags. Mixed orientations have no inferred full domain.",
            })
    if rows["distinct"]["TIME"]["count"] != 1:
        issues.append("Multiple TIME values require later multi-integration support")
    if rows["context_count_excluding_pair_and_data_description"] != 1:
        issues.append("Mixed row contexts require separate partitions")
    if rows["row_key_duplicate_count"]:
        issues.append("Duplicate row keys require explicit reconstruction maps")
    if len(configs) != 1:
        issues.append("Multiple data descriptions require explicit partitions/channel maps")
    for name in selected:
        item = main_inventory["columns"][name]
        if len(item["cell_shapes"]) != 1 or item["undefined_cells"]:
            issues.append(f"{name}: heterogeneous or undefined cells require explicit mappings")
    if len(configs) == 1 and configs[0]["channel_count"] > 1:
        organization = "channels_inside_cells"
    elif configs and all(c["channel_count"] == 1 for c in configs):
        organization = "single_channel_data_description_records"
    else:
        organization = "mixed_channel_configurations"
    return {"organization": organization, "configurations": configs,
            "regular_single_configuration_candidate": not issues,
            "initial_regular_profile_issues": sorted(set(issues)),
            "archive_implementation_available": True,
            "quantized_profile_supported": not issues and organization == "channels_inside_cells",
            "channel_presence_note": "Per-cell dimensions were compared with spectral-window metadata. Unrecorded channels beyond that metadata cannot be inferred.",
            "validation_scope": "Inspection only. This is not an MS conformance or CASA compatibility certification."}


def inspect_ms(source, selected=("DATA",), chunk_rows=1024):
    """Return a strict JSON-compatible inspection report, or raise on read failure."""
    source = Path(source).resolve()
    selected = tuple(selected)
    if chunk_rows < 1:
        raise ValueError("chunk_rows must be positive")
    if not selected or len(selected) != len(set(selected)):
        raise ValueError("Selected columns must be nonempty and unique")
    if not (source / "table.dat").is_file():
        raise ValueError(f"Not a casacore table directory: {source}")
    started = time.perf_counter()
    subtables, subtable_paths, managers = {}, {}, {}
    with _open(source) as t:
        if t.nrows() == 0:
            raise ValueError("Empty MAIN table is not supported by this inspection profile")
        missing = [c for c in selected if c not in t.colnames()]
        required_columns = ("UVW", "FLAG", "FLAG_ROW", "ANTENNA1", "ANTENNA2", "TIME", "TIME_CENTROID", "DATA_DESC_ID")
        missing.extend(c for c in required_columns if c not in t.colnames())
        if missing:
            raise ValueError(f"Missing required or selected columns: {', '.join(sorted(set(missing)))}")
        if any(c not in VISIBILITY_COLUMNS for c in selected):
            raise ValueError("Select DATA, FLOAT_DATA, CORRECTED_DATA or MODEL_DATA explicitly")
        for name, value in t.getkeywords().items():
            if not isinstance(value, str) or not value.startswith("Table: "):
                continue
            path = Path(value[7:])
            path = (source / path).resolve() if not path.is_absolute() else path.resolve()
            if not path.is_relative_to(source):
                raise ValueError(f"External subtable reference is unsupported: {name} -> {path}")
            subtable_paths[name] = path
        for name, path in subtable_paths.items():
            with _open(path) as sub:
                subtables[name] = _scan_table(sub, chunk_rows, include_rows=True)
                managers[str(path.relative_to(source))] = sub.getdminfo()
        main = _scan_table(t, chunk_rows)
        managers["MAIN"] = t.getdminfo()
        row_inventory, row_columns = _row_inventory(t)
        layout = _layout(t, row_inventory, row_columns, subtable_paths, selected, main)
        visibility = _visibility_inventory(t, selected, chunk_rows)
        uvw = np.asarray(t.getcol("UVW"))
        if uvw.shape != (t.nrows(), 3) or uvw.dtype != np.dtype("float64"):
            raise ValueError(f"UVW must be float64[row,3], found {uvw.dtype} {uvw.shape}")
        unique_uvw, _, _ = unique_bits(uvw)
        uvw_inventory = {
            "dtype": uvw.dtype.str, "shape": list(uvw.shape),
            "units": json_value(t.getcolkeywords("UVW").get("QuantumUnits")),
            "measure": json_value(t.getcolkeywords("UVW").get("MEASINFO")),
            "bitwise_distinct_geometry_records": len(unique_uvw),
            "logical_array_bytes": uvw.nbytes,
            "dictionary_plus_uint32_row_map_logical_bytes": unique_uvw.nbytes + len(uvw) * 4,
            "note": "Logical dictionary estimate only; no measured archive savings. Only identical bytes are shareable. No UVW regeneration or wavelength conversion was performed.",
        }
    report = {
        "schema": "varquanms.inspection/1", "status": "complete",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": str(source), "selected_visibility_columns": list(selected),
        "other_visibility_columns_preserved_by_future_archive": [c for c in VISIBILITY_COLUMNS if c in main["columns"] and c not in selected],
        "flags_weights_columns_present": [c for c in FLAG_WEIGHT_COLUMNS if c in main["columns"]],
        "software": {"varquanMS": __version__, "python": platform.python_version(),
                     "numpy": np.__version__, "python_casacore": casacore.__version__},
        "main": main, "subtables": subtables, "row_identity": row_inventory,
        "layout": layout, "uvw": uvw_inventory, "visibility": visibility,
        "storage": _physical_inventory(source, managers),
        "policies": {"source_access": "readonly/usernoread; requires a frozen input copy",
                     "inspection_scan": "All defined values, all rows; no sampling",
                     "numeric_metadata_encoding": "dtype + shape + base64 raw C-order bytes; values are a human-readable preview",
                     "logical_string_bytes": "UTF-8 value lengths; excludes Python or container overhead",
                     "noise_model": "None validated; WEIGHT/SIGMA do not imply a thermal-noise guarantee",
                     "archive_status": "Lossless reconstruction and regular-layout uniform quantization are implemented; scientific acceptance is separate"},
        "elapsed_seconds": time.perf_counter() - started,
    }
    # Validate portability before any output file is finalized.
    json.dumps(report, allow_nan=False)
    return report


def human_summary(report):
    main, layout = report["main"], report["layout"]
    lines = [f"Inspection: {Path(report['source']).name}",
             f"Rows: {main['row_count']:,}; organization: {layout['organization']}",
             f"Selected visibility columns: {', '.join(report['selected_visibility_columns'])}"]
    for c in layout["configurations"]:
        lines.append(f"DDID {c['data_description_id']}: {c['channel_count']} channels, {c['correlation_count']} products ({', '.join(c['correlation_names'])}), {c['ordered_pair_count']:,} ordered antenna pairs")
    for name in ("TIME", "TIME_CENTROID"):
        distinct = report["row_identity"]["distinct"][name]
        lines.append(f"{name}: {distinct['count']} bitwise-distinct value(s); {distinct['values'].get('values', 'see report')}; reference {main['columns'][name]['measure']}")
    lines += [
        f"Autocorrelation rows: {report['row_identity']['autocorrelation_rows']:,}; duplicate row keys: {report['row_identity']['row_key_duplicate_count']:,}",
        f"UVW: {report['uvw']['bitwise_distinct_geometry_records']:,} bitwise-distinct vectors; {report['uvw']['logical_array_bytes']:,} logical bytes; {report['uvw']['measure']}",
        f"Physical file lengths: {report['storage']['file_bytes']:,} bytes; allocated file blocks: {report['storage']['allocated_file_bytes']:,} bytes",
        f"MAIN logical arrays: {main['logical_array_bytes']:,} bytes; subtables: {sum(s['logical_array_bytes'] for s in report['subtables'].values()):,} bytes",
    ]
    for name, vis in report["visibility"].items():
        stats = main["columns"][name]["statistics"]
        lines.append(f"{name}: {main['columns'][name]['logical_array_bytes']:,} logical bytes; {vis.get('flagged_samples_including_flag_row', 0):,} flagged samples; {stats.get('nonfinite_values', 0):,} nonfinite samples; {stats.get('negative_zero_components', 0):,} signed-zero components")
    if layout["initial_regular_profile_issues"]:
        lines.append("Initial regular-profile issues: " + "; ".join(layout["initial_regular_profile_issues"]))
    else:
        lines.append("Layout meets the initial regular quantization profile; acceptance still requires encode/decode/source verification.")
    lines.append("Storage files may contain multiple columns; logical repetition does not establish physical savings.")
    lines.append("This inspection alone establishes no compression ratio, noise guarantee, scientific fidelity or CASA workflow compatibility.")
    return "\n".join(lines) + "\n"
