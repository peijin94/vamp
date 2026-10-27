"""Explicitly destructive MAIN-flag policy; no visibility values are discarded."""

import json

import numpy as np

FLAG_COLUMNS = ("FLAG", "FLAG_CATEGORY", "FLAG_ROW")
ENCODING = "burned-flags-v1"
POLICY_KEYWORD = "VARQUANMS_FLAG_POLICY"


def export_marker(columns):
    return {"policy": ENCODING, "columns": list(columns),
            "reconstruction": "Defined MAIN flag cells reset to false; original flag values discarded",
            "reflagging_performed": False}


def burned_columns(manifest):
    """Validate the opt-in against every column specification, including old archives."""
    enabled = manifest["config"].get("burn_flag", False)
    if type(enabled) is not bool:
        raise ValueError("burn_flag must be boolean")
    names = []
    for table in manifest["tables"]:
        for name, spec in table["columns"].items():
            expected = enabled and table["path"] == "." and name in FLAG_COLUMNS
            if (spec["encoding"] == ENCODING) != expected:
                raise ValueError("Burned flag encoding disagrees with explicit MAIN flag policy")
            if expected:
                names.append(name)
    return names


def policy_summary(manifest):
    columns = burned_columns(manifest)
    enabled = manifest["config"].get("burn_flag", False)
    return {"burn_flag": enabled, "discarded_main_columns": columns,
            "original_main_flag_values_preserved": not enabled,
            "export": "all defined MAIN flag cells false; undefined cells remain undefined" if enabled else "original flag values",
            "reflagging_performed": False, "subtable_flags": "preserved exactly",
            "flagged_visibility_values": ("discarded; CORRECTED_DATA reconstructs as complex positive zero"
                                          if manifest["config"].get("keep_calib_only", False)
                                          else "preserved bitwise using source flags during encoding"),
            "warning": "Destructive flag policy; reflag before scientific use. Recomputed flags may differ from the original." if enabled else None}


def write_burned_flags(writer, t, column, prefix, info, config):
    """Store only flag-cell shapes/presence, never their original Boolean values."""
    if t.getcoldesc(column)["valueType"] != "boolean":
        raise ValueError(f"burn_flag requires Boolean column {column}")
    spec = {"encoding": ENCODING, "row_count": t.nrows(), "dtype": "|b1"}
    if t.nrows() == 0:
        return dict(spec, cell_shape=[])
    if info["undefined_cells"] == 0 and len(info["cell_shapes"]) == 1:
        return dict(spec, cell_shape=json.loads(next(iter(info["cell_shapes"]))))
    chunks = []
    for start in range(0, t.nrows(), config.io_chunk_rows):
        stop = min(t.nrows(), start + config.io_chunk_rows)
        shapes = [list(np.asarray(t.getcell(column, row)).shape) if t.iscelldefined(column, row) else None
                  for row in range(start, stop)]
        path = writer.put_json(prefix + f"/shapes_{start}", shapes, "reconstruction")
        chunks.append({"start": start, "stop": stop, "path": path})
    return dict(spec, shape_chunks=chunks)


def _zeros(shape, rows=None):
    if not isinstance(shape, list) or any(type(n) is not int or n < 0 for n in shape):
        raise ValueError("Invalid burned flag cell shape")
    if np.prod(shape, dtype=object) > 64 * 1024**2:
        raise ValueError("Burned flag cell exceeds the 64 MiB prototype limit")
    return np.zeros(shape if rows is None else [rows, *shape], dtype=bool)


def read_burned_flags(reader, spec, start, stop):
    if spec.get("dtype") != "|b1" or ("cell_shape" in spec) == ("shape_chunks" in spec):
        raise ValueError("Invalid burned flag reconstruction specification")
    if "cell_shape" in spec:
        return _zeros(spec["cell_shape"], stop - start)
    values = []
    cursor = start
    for chunk in spec["shape_chunks"]:
        if chunk["stop"] <= start or chunk["start"] >= stop:
            continue
        lo, hi = max(start, chunk["start"]), min(stop, chunk["stop"])
        shapes = reader.json(chunk["path"])
        if lo != cursor or len(shapes) != chunk["stop"] - chunk["start"]:
            raise ValueError("Burned flag shape chunks have missing/overlapping rows")
        values.extend(None if shape is None else _zeros(shape)
                      for shape in shapes[lo - chunk["start"]:hi - chunk["start"]])
        cursor = hi
    if cursor != stop:
        raise ValueError("Burned flag shape chunks have missing rows")
    return values
