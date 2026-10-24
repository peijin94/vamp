"""Lossless MS values, relocatable metadata and deterministic MS reconstruction."""

from collections import Counter
from datetime import datetime, timezone
import json
import hashlib
import os
from pathlib import Path
import platform
import resource
import tempfile
import time

import casacore
from casacore.tables import table
import h5py
import numpy as np

from . import __version__
from .config import Config
from .calibrated import CalibratedView, OMITTED, descriptor as calibrated_descriptor, marker as calibrated_marker
from .calibrated import DATA_KEYWORD, data_marker, materialize_ms_v2_data
from .weights import COLUMNS as WEIGHT_COLUMNS, ENCODING as WEIGHT_ENCODING, read_weights, write_weights, verify_weights
from .container import Reader, Writer, new_output
from .flags import (ENCODING as BURNED_FLAGS, FLAG_COLUMNS, POLICY_KEYWORD,
                    burned_columns, export_marker, policy_summary,
                    read_burned_flags, write_burned_flags)
from .inspection import _batches, _open, inspect_ms, unique_bits
from .metadata import casacore_keywords, dumps, equal, logical_descriptor, pack, portable, references, relocate, unpack

SUPPORTED_TYPES = {"boolean", "uchar", "short", "ushort", "int", "uint", "int64", "float", "double", "complex", "dcomplex", "string"}


def discover_tables(source):
    source = Path(source).resolve()
    pending = [source]
    discovered = {}
    while pending:
        path = pending.pop(0)
        if path in discovered:
            continue
        if not path.is_relative_to(source) or path.is_symlink():
            raise ValueError(f"Unsupported external/symlink table: {path}")
        with _open(path) as t:
            description = t.getdesc()
            for name, value in description.get("_private_keywords_", {}).items():
                # Standard casacore hypercolumn declarations are descriptive,
                # and the table constructor restores them with StandardStMan.
                # Unknown private extensions still fail before archive creation.
                if (not name.startswith("Hypercolumn_") or not isinstance(value, dict)
                        or set(value) != {"ndim", "data", "coord", "id"}
                        or type(value["ndim"]) is not int or value["ndim"] < 1
                        or any(not isinstance(value[k], (list, np.ndarray))
                               or any(c not in t.colnames() for c in value[k])
                               for k in ("data", "coord", "id"))):
                    raise ValueError("Unsupported private casacore table keyword: " + name)
            for column in t.colnames():
                if t.getcoldesc(column)["valueType"] not in SUPPORTED_TYPES:
                    raise ValueError(f"Unsupported column type: {path.name}/{column}")
            normalized = relocate(description, source)
            schema = {"path": path.relative_to(source).as_posix(), "row_count": t.nrows(),
                      "column_names": t.colnames(), "descriptor": pack(normalized),
                      "info": pack(t.info()), "original_storage_managers": pack(t.getdminfo())}
            # Validate every schema before any incomplete archive is created.
            dumps(schema)
            discovered[path] = schema
            for reference in references(description):
                target = Path(reference)
                target = (source / target).resolve() if not target.is_absolute() else target.resolve()
                if not target.is_relative_to(source):
                    raise ValueError(f"External table reference: {reference}")
                pending.append(target)
    expected = {p.parent.resolve() for p in source.rglob("table.dat")}
    if set(discovered) != expected:
        raise ValueError("Unreferenced nested casacore tables require explicit support")
    for path in source.rglob("*"):
        if path.is_symlink():
            raise ValueError(f"Symlink inside source MS: {path}")
        if path.is_file() and (path.parent.resolve() not in discovered or not path.name.startswith("table.")):
            raise ValueError(f"Non-table auxiliary file requires explicit support: {path}")
    return discovered


def _exact_role(table_path, column, selected):
    if table_path != ".":
        return "subtables"
    if column in selected:
        return "visibility_lossless"
    if column in ("FLAG", "FLAG_ROW", "FLAG_CATEGORY", "WEIGHT", "SIGMA", "WEIGHT_SPECTRUM", "SIGMA_SPECTRUM"):
        return "flags_weights"
    if column in ("UVW", "TIME", "TIME_CENTROID", "ANTENNA1", "ANTENNA2"):
        return "coordinates"
    return "auxiliary_columns"


def _write_exact_column(writer, t, column, prefix, info, config, role, deduplicate):
    rows = t.nrows()
    if rows == 0:
        return {"encoding": "empty", "row_count": 0}
    # Numeric homogeneous columns can be represented as one typed array.
    if info["undefined_cells"] == 0 and len(info["cell_shapes"]) == 1 and len(info["numpy_dtypes"]) == 1 and np.dtype(info["numpy_dtypes"][0]).kind in "biufc":
        shape = tuple(json.loads(next(iter(info["cell_shapes"]))))
        dtype = np.dtype(info["numpy_dtypes"][0]).newbyteorder("<")
        spec = {"encoding": "dense", "row_count": rows, "cell_shape": list(shape), "dtype": dtype.str}
        if deduplicate and rows > 1 and info["logical_array_bytes"] <= config.max_visibility_bytes:
            values = portable(t.getcol(column))
            if values[0].size:
                bytes_by_row = values.view(np.uint8).reshape(rows, -1)
                if np.all(bytes_by_row == bytes_by_row[0]):
                    spec.update(encoding="constant", values=writer.put(prefix + "/values", values[:1], role))
                    return spec
                distinct, inverse, _ = unique_bits(values)
                if len(distinct) <= np.iinfo(np.uint32).max and distinct.nbytes + rows * 4 < values.nbytes:
                    spec.update(encoding="dictionary",
                                values=writer.put(prefix + "/values", distinct, role),
                                row_map=writer.put(prefix + "/row_map", inverse.astype("<u4"), "reconstruction"))
                    return spec
        dataset = writer.create(prefix + "/values", (rows, *shape), dtype, role)
        for start in range(0, rows, config.io_chunk_rows):
            stop = min(rows, start + config.io_chunk_rows)
            dataset[start:stop] = portable(t.getcol(column, startrow=start, nrow=stop - start))
        spec["values"] = prefix + "/values"
        return spec
    # Explicit cell records retain undefined versus empty and arbitrary shapes.
    chunks = []
    for start in range(0, rows, config.io_chunk_rows):
        stop = min(rows, start + config.io_chunk_rows)
        records = [pack(np.asarray(t.getcell(column, row))) if t.iscelldefined(column, row) else ["undefined"]
                   for row in range(start, stop)]
        path = prefix + f"/cells_{start}"
        writer.put_json(path, records, role)
        chunks.append({"start": start, "stop": stop, "path": path})
    return {"encoding": "cells", "row_count": rows, "chunks": chunks}


def read_column(reader, spec, start, stop):
    """Return a numeric row slab, or explicit cells (None means undefined)."""
    if not 0 <= start <= stop <= spec["row_count"]:
        raise ValueError("Column row bounds out of range")
    encoding = spec["encoding"]
    if encoding == BURNED_FLAGS:
        return read_burned_flags(reader, spec, start, stop)
    if encoding == WEIGHT_ENCODING:
        return read_weights(reader, spec, start, stop)
    if encoding == "empty":
        return []
    if encoding == "dense":
        return reader.h5[spec["values"]][start:stop]
    if encoding in ("constant", "dictionary"):
        if not hasattr(reader, "column_cache"):
            reader.column_cache = {}
        if spec["values"] not in reader.column_cache:
            reader.column_cache[spec["values"]] = reader.h5[spec["values"]][...]
        values = reader.column_cache[spec["values"]]
        if encoding == "constant":
            return np.repeat(values, stop - start, axis=0)
        indices = reader.h5[spec["row_map"]][start:stop]
        if np.any(indices >= len(values)):
            raise ValueError("Dictionary reconstruction index out of range")
        return values[indices]
    if encoding == "cells":
        result = []
        for chunk in spec["chunks"]:
            if chunk["stop"] <= start or chunk["start"] >= stop:
                continue
            nodes = reader.json(chunk["path"])
            if len(nodes) != chunk["stop"] - chunk["start"]:
                raise ValueError("Cell-record chunk length mismatch")
            lo, hi = max(start, chunk["start"]) - chunk["start"], min(stop, chunk["stop"]) - chunk["start"]
            result.extend(None if v == ["undefined"] else unpack(v) for v in nodes[lo:hi])
        if len(result) != stop - start:
            raise ValueError("Cell-record reconstruction has missing rows")
        return result
    if encoding == "quantized":
        from .visibility import read_visibility
        return read_visibility(reader, spec)[start:stop]
    raise ValueError(f"Unknown column encoding: {encoding}")


def _environment():
    from .codec import _native
    return {"varquanMS": __version__, "python": platform.python_version(),
            "numpy": np.__version__, "python_casacore": casacore.__version__, "h5py": h5py.__version__,
            "platform": platform.platform(), "machine": platform.machine(),
            "cpu_count": os.cpu_count(), "threads_environment": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
            "native_acceleration": _native is not None,
            "implementation_files_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(Path(__file__).parent.iterdir()) if p.is_file() and p.suffix in (".py", ".c", ".so", ".pyd")}}


def encode(source, output, config=None, inspection=None, *, distribution_references=None):
    config = config or Config()
    config.validate()
    source = Path(source).resolve()
    output = new_output(output, source)
    started = time.perf_counter()
    report = inspection or inspect_ms(source, config.columns, config.io_chunk_rows)
    if report["source"] != str(source) or report["selected_visibility_columns"] != config.columns:
        raise ValueError("Inspection does not match the requested source/columns")
    schemas = discover_tables(source)
    removed_columns = []
    if config.keep_calib_only:
        main_schema = schemas[source]
        removed_columns = [c for c in main_schema["column_names"] if c in OMITTED]
        main_schema["column_names"] = [c for c in main_schema["column_names"] if c not in OMITTED]
        main_schema["descriptor"] = pack(calibrated_descriptor(unpack(main_schema["descriptor"]), config.to_dict()))
        if config.weight_bits == 4:
            for name in WEIGHT_COLUMNS:
                info = report["main"]["columns"].get(name)
                if info is not None and (info["undefined_cells"] or len(info["cell_shapes"]) != 1 or info["numpy_dtypes"] != [np.dtype("float32").str]):
                    raise ValueError("Four-bit weights require homogeneous defined float32 cells; use weight_bits=32 for " + name)
    if config.burn_flag:
        main = report["main"]
        columns = [name for name in main["columns"] if name in FLAG_COLUMNS]
        for name in columns:
            if main["descriptor"][name]["valueType"] != "boolean":
                raise ValueError(f"burn_flag requires Boolean column {name}")
        existing = unpack(schemas[source]["descriptor"]).get("_keywords_", {}).get(POLICY_KEYWORD)
        if existing is not None and not equal(existing, export_marker(columns)):
            raise ValueError(f"Cannot overwrite existing {POLICY_KEYWORD} provenance keyword")
    if config.mode != "lossless":
        issues = report["layout"]["initial_regular_profile_issues"]
        if issues or report["layout"]["organization"] != "channels_inside_cells":
            raise ValueError("Quantized profile requires one regular configuration: " + "; ".join(issues or [report["layout"]["organization"]]))
        names = report["layout"]["configurations"][0]["correlation_names"]
        if set(config.product_bounds) - set(names):
            raise ValueError("Per-product bounds name a product absent from POLARIZATION")
        for c in config.columns:
            if report["main"]["columns"][c]["numpy_dtypes"] != [np.dtype("complex64").str]:
                raise ValueError("Quantized profile currently requires complex64 cells")
            if report["main"]["columns"][c]["logical_array_bytes"] > config.max_visibility_bytes:
                raise ValueError("Selected visibility array exceeds max_visibility_bytes prototype limit")
    if config.normalized_bounds is not None and distribution_references is None:
        from .references import sample_references
        distribution_references = sample_references(source, config)
    preflight_seconds = time.perf_counter() - started
    writer = Writer(output, config)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(), "config": config.to_dict(),
        "source": {"name": source.name, "physical_file_bytes": report["storage"]["file_bytes"],
                   "selected_visibility_logical_bytes": sum(report["main"]["columns"][c]["logical_array_bytes"] for c in config.columns),
                   "layout": report["layout"], "main_row_count": report["main"]["row_count"],
                   "omitted_visibility_columns": removed_columns},
        "environment": _environment(), "tables": [], "timings_seconds": {"inspection_preflight": preflight_seconds},
        "provenance": {"source_modified": False, "source_loss_history": "Unknown; original storage-manager assignments are recorded per table",
                       "compression_generation": 1, "noise_model": config.noise,
                       "science_validation": "Not performed; absolute numeric error bounds are not scientific acceptance",
                       "storage_manager_layout_preserved": False,
                       "source_snapshot_requirement": "Frozen input; readonly usernoread access"},
    }
    exact_started = time.perf_counter()
    if config.normalized_bounds is not None:
        manifest['distribution_references'] = distribution_references
    try:
        for table_index, (path, schema) in enumerate(schemas.items()):
            item = dict(schema, columns={})
            table_name = schema["path"]
            scan = report["main"] if table_name == "." else report["subtables"].get(path.name)
            if scan is None:
                raise ValueError("Nested subtable is not covered by the inspection profile")
            with _open(path) as t:
                if table_name == "." and config.keep_calib_only:
                    t = CalibratedView(t, config.to_dict())
                for column_index, column in enumerate(t.colnames()):
                    prefix = f"tables/t{table_index}/c{column_index}"
                    if table_name == "." and config.keep_calib_only and config.weight_bits == 4 and column in WEIGHT_COLUMNS:
                        item["columns"][column] = write_weights(writer, t, column, prefix, config, scan["columns"][column])
                    elif table_name == "." and config.burn_flag and column in FLAG_COLUMNS:
                        item["columns"][column] = write_burned_flags(
                            writer, t, column, prefix, scan["columns"][column], config)
                    elif table_name == "." and column in config.columns and config.mode != "lossless":
                        from .visibility import write_visibility
                        item["columns"][column] = write_visibility(writer, t, column, prefix, config, report,
                                                                   distribution_references=distribution_references)
                    else:
                        role = _exact_role(table_name, column, config.columns)
                        item["columns"][column] = _write_exact_column(
                            writer, t, column, prefix, scan["columns"][column], config,
                            role, config.deduplicate and role != "visibility_lossless")
            manifest["tables"].append(item)
        manifest["provenance"]["flag_policy"] = policy_summary(manifest)
        manifest["provenance"]["calibrated_only_policy"] = calibrated_marker() if config.keep_calib_only else None
        manifest["timings_seconds"]["column_encoding_and_writing"] = time.perf_counter() - exact_started
        manifest["peak_process_rss_bytes"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024
        manifest["timings_seconds"]["before_integrity_finalization"] = time.perf_counter() - started
        manifest = writer.finish(manifest)
    finally:
        writer.close()
    return archive_summary(output, manifest, time.perf_counter() - started)


def archive_summary(path, manifest, elapsed=None):
    categories = dict(manifest["dataset_stored_bytes_by_category"])
    total = Path(path).stat().st_size
    categories["container_and_manifest_overhead"] = total - sum(categories.values())
    original = manifest["source"]["physical_file_bytes"]
    visibility_logical = manifest["source"]["selected_visibility_logical_bytes"]
    visibility_stored = sum(v for k, v in categories.items() if k.startswith("visibility_") or k in ("permutations", "quantization_statistics"))
    return {"archive": str(path), "archive_bytes": total, "source_physical_file_bytes": original,
            "whole_archive_ratio": original / total if total else None,
            "selected_visibility_logical_bytes": visibility_logical,
            "visibility_datasets_stored_bytes": visibility_stored,
            "visibility_dataset_ratio_including_maps_statistics": visibility_logical / visibility_stored if visibility_stored else None,
            "stored_bytes_by_category": categories, "elapsed_seconds": elapsed,
            "mode": manifest["config"]["mode"], "deduplicate": manifest["config"]["deduplicate"],
            "flag_policy": policy_summary(manifest),
            "keep_calib_only": manifest["config"].get("keep_calib_only", False),
            "omitted_visibility_columns": manifest["source"].get("omitted_visibility_columns", []),
            "weights_quality": {c: s["quality"] for t in manifest["tables"] if t["path"] == "." for c,s in t["columns"].items() if s["encoding"] == WEIGHT_ENCODING},
            "exact_ms_round_trip": manifest["config"]["mode"] == "lossless" and not manifest["config"].get("burn_flag", False) and not manifest["config"].get("keep_calib_only", False),
            "timings_seconds": manifest["timings_seconds"],
            "peak_process_rss_bytes": manifest["peak_process_rss_bytes"],
            "visibility_quality": {c: s.get("quality_summary") for t in manifest["tables"] if t["path"] == "." for c, s in t["columns"].items() if s["encoding"] == "quantized"},
            "visibility_block_statistics": {c: {"products": s["products"], "timings_seconds": s["timings_seconds"],
                                                "worker_timings_seconds": s.get("worker_timings_seconds", {}),
                                                "worker_count": s.get("worker_count", 1)}
                                           for t in manifest["tables"] if t["path"] == "." for c, s in t["columns"].items() if s["encoding"] == "quantized"},
            "note": "Visibility dataset ratio includes indexes, mappings, counts and quality datasets, but container-wide overhead is reported only in whole-archive bytes. No scientific or CASA workflow acceptance is implied."}


def decode(archive, output, ms_v2_data=None, *, group=None):
    output = new_output(output, archive)
    started = time.perf_counter()
    with Reader(archive, group=group) as reader:
        if ms_v2_data is None:
            ms_v2_data = reader.manifest['config'].get('keep_calib_only', False)
        if ms_v2_data and not reader.manifest['config'].get('keep_calib_only', False):
            raise ValueError('--ms-v2-data is only for a calibrated-only archive')
        staged = Path(tempfile.mkdtemp(prefix=output.name + ".", suffix=".partial", dir=output.parent))
        schemas = reader.manifest["tables"]
        # Create the entire table graph before restoring any table references.
        for spec in schemas:
            relative = Path(spec["path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise ValueError("Unsafe reconstruction table path")
            target = staged / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            description = unpack(spec["descriptor"])
            description.pop("_keywords_", None)
            for name in spec["column_names"]:
                description[name] = dict(description[name])
                description[name].pop("keywords", None)
                description[name]["dataManagerType"] = "StandardStMan"
                description[name]["dataManagerGroup"] = "StandardStMan"
            with table(str(target), tabledesc=description, nrow=spec["row_count"], readonly=False, ack=False) as t:
                info = unpack(spec["info"])
                if info.get("readme", "").endswith("\n"):
                    # CTDS putinfo appends a line terminator to a nonempty readme.
                    info["readme"] = info["readme"][:-1]
                t.putinfo(info)
                for name, column in spec["columns"].items():
                    if column["encoding"] == "quantized":
                        values = read_column(reader, column, 0, spec["row_count"])
                        for start in range(0, spec["row_count"], 1024):
                            chunk = values[start:start + 1024]
                            t.putcol(name, chunk, startrow=start, nrow=len(chunk))
                    else:
                        for start in range(0, spec["row_count"], 1024):
                            values = read_column(reader, column, start, min(spec["row_count"], start + 1024))
                            if isinstance(values, np.ndarray):
                                if values.size:
                                    t.putcol(name, values, startrow=start, nrow=len(values))
                                else:
                                    for index, value in enumerate(values):
                                        t.putcell(name, start + index, value)
                            else:
                                for index, value in enumerate(values):
                                    if value is not None:
                                        cell = value.item() if isinstance(value, np.ndarray) and value.ndim == 0 else value
                                        if isinstance(cell, np.ndarray) and cell.size == 0 and description[name]["valueType"] == "string":
                                            # casacore returns an untyped [] for empty string
                                            # vectors; np.asarray([]) is float64. Restore the
                                            # declared column type when crossing its API.
                                            cell = np.asarray([], dtype="U1").reshape(cell.shape)
                                        t.putcell(name, start + index, cell)
        for spec in schemas:
            with table(str(staged / spec["path"]), readonly=False, ack=False) as t:
                description = relocate(unpack(spec["descriptor"]), staged, restore=True)
                t.putkeywords(casacore_keywords(description.get("_keywords_", {})))
                for name in spec["column_names"]:
                    t.putcolkeywords(name, casacore_keywords(description[name].get("keywords", {})))
                if spec["path"] == "." and reader.manifest["config"].get("burn_flag", False):
                    t.putkeyword(POLICY_KEYWORD, export_marker(burned_columns(reader.manifest)))
        if ms_v2_data:
            materialize_ms_v2_data(staged)
        # Reopen every table prior to publication; reconstruction references must resolve.
        for spec in schemas:
            with _open(staged / spec["path"]) as t:
                if t.nrows() != spec["row_count"]:
                    raise ValueError("Export row count mismatch")
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite {output}")
        staged.rename(output)
    return {"output": str(output), "status": "exported", "elapsed_seconds": time.perf_counter() - started,
            "derived_calibrated_DATA": ms_v2_data,
            "python_casacore_reopen": "passed", "CASA_workflow": "not run",
            "flag_policy": policy_summary(reader.manifest)}


def _check_column(reader, spec, t, column, shapes_only=False):
    """Check exact reconstruction, or only shapes/presence for discarded source flags."""
    def matches(actual, expected):
        if shapes_only:
            return actual.shape == expected.shape and actual.dtype == expected.dtype
        return equal(actual, expected)

    for start in range(0, t.nrows(), 1024):
        stop = min(t.nrows(), start + 1024)
        values = read_column(reader, spec, start, stop)
        if isinstance(values, np.ndarray):
            actual = np.asarray(t.getcol(column, startrow=start, nrow=stop - start))
            if not matches(actual, values):
                raise ValueError(f"Column reconstruction mismatch: {column}, rows {start}:{stop}")
        else:
            for index, value in enumerate(values):
                row = start + index
                defined = t.iscelldefined(column, row)
                if defined != (value is not None) or defined and not matches(np.asarray(t.getcell(column, row)), value):
                    raise ValueError(f"Column reconstruction mismatch: {column}, row {row}")


def verify(archive, source=None, restored=None):
    started = time.perf_counter()
    checks = []
    with Reader(archive) as reader:
        for spec in reader.manifest["tables"]:
            if source is None:
                continue
            root = Path(source).resolve()
            with _open(root / spec["path"]) as t:
                if spec["path"] == "." and reader.manifest["config"].get("keep_calib_only", False):
                    t = CalibratedView(t, reader.manifest["config"])
                expected = unpack(spec["descriptor"])
                actual = relocate(t.getdesc(), root)
                if t.nrows() != spec["row_count"] or t.colnames() != spec["column_names"] or not equal(logical_descriptor(expected), logical_descriptor(actual)) or not equal(unpack(spec["info"]), t.info()):
                    raise ValueError(f"Source schema/metadata mismatch: {spec['path']}")
                for column, encoding in spec["columns"].items():
                    if encoding["encoding"] == WEIGHT_ENCODING:
                        quality = verify_weights(reader, encoding, t, column, reader.manifest["config"])
                        checks.append({"table": spec["path"], "column": column,
                                       "status": "bitwise_exact" if quality["original_bits_exact"] else "weight_error_bound_passed", "quality": quality})
                        continue
                    if encoding["encoding"] == "quantized":
                        from .visibility import verify_visibility
                        quality = verify_visibility(reader, encoding, t, column, reader.manifest["config"])
                        checks.append({"table": spec["path"], "column": column, "status": "measured_error_bounds_passed", "quality": quality})
                        continue
                    burned = encoding["encoding"] == BURNED_FLAGS
                    _check_column(reader, encoding, t, column, shapes_only=burned)
                    checks.append({"table": spec["path"], "column": column,
                                   "status": "intentionally_discarded_shapes_preserved" if burned else "bitwise_exact"})
        if restored is not None:
            if source is None:
                raise ValueError("--restored comparison requires --source")
            lossy_weights = [c for t in reader.manifest["tables"] if t["path"] == "." for c,s in t["columns"].items() if s["encoding"] == WEIGHT_ENCODING and not s["quality"]["original_bits_exact"]]
            checks.extend(compare_ms(source, restored,
                                     (reader.manifest["config"]["columns"] if reader.manifest["config"]["mode"] != "lossless" else []) + lossy_weights,
                                     burned_flag_columns=burned_columns(reader.manifest),
                                     keep_calib_only=reader.manifest["config"].get("keep_calib_only", False),
                                     calib_config=reader.manifest["config"]))
            for spec in reader.manifest["tables"]:
                if spec["path"] != ".":
                    continue
                with _open(Path(restored).resolve()) as exported:
                    for name, column in spec["columns"].items():
                        if column["encoding"] == BURNED_FLAGS:
                            _check_column(reader, column, exported, name)
                            checks.append({"table": ".", "column": name, "export": "all_false_shapes_and_presence_preserved"})
                        if column["encoding"] in ("quantized", WEIGHT_ENCODING):
                            values = read_column(reader, column, 0, spec["row_count"])
                            for start in range(0, spec["row_count"], 1024):
                                expected = values[start:start + 1024]
                                if not equal(expected, np.asarray(exported.getcol(name, startrow=start, nrow=len(expected)))):
                                    raise ValueError(f"Exported quantized values differ from archive decoder: {name}")
                            checks.append({"table": ".", "column": name, "export": "exactly_matches_archive_reconstruction"})
        result = {"schema": "varquanms.verification/1", "status": "passed",
                  "archive": str(archive), "integrity": "all inventory SHA-256 and HDF5 chunk checksums passed",
                  "source": str(source) if source else None, "restored": str(restored) if restored else None,
                  "column_checks": checks, "summary": archive_summary(archive, reader.manifest),
                  "flag_policy": policy_summary(reader.manifest),
                  "verification_scope": "declared archive policy; intentionally discarded flags are not checked for source-value equality",
                  "CASA_workflow": "not run", "scientific_products": "not tested",
                  "elapsed_seconds": time.perf_counter() - started}
    return result


def compare_ms(source, restored, approximate_columns=(), burned_flag_columns=(), keep_calib_only=False, calib_config=None):
    if set(burned_flag_columns) - set(FLAG_COLUMNS):
        raise ValueError("Only MAIN flag columns can be intentionally discarded")
    source, restored = Path(source).resolve(), Path(restored).resolve()
    schemas = discover_tables(source)
    output_schemas = discover_tables(restored)
    if {s["path"] for s in schemas.values()} != {s["path"] for s in output_schemas.values()}:
        raise ValueError("Exported subtable graph differs")
    result = []
    for path, schema in schemas.items():
        relative = schema["path"]
        with _open(path) as a, _open(restored / relative) as b:
            if relative == "." and keep_calib_only:
                a = CalibratedView(a, calib_config)
            expected = relocate(a.getdesc(), source)
            names = a.colnames()
            derived = relative == '.' and keep_calib_only and DATA_KEYWORD in b.getkeywords()
            if derived:
                if not equal(b.getkeyword(DATA_KEYWORD), data_marker()):
                    raise ValueError('Invalid derived DATA provenance')
                expected['DATA'] = expected['CORRECTED_DATA']
                expected.setdefault('_keywords_', {})[DATA_KEYWORD] = data_marker()
                names = names + ['DATA']
            if relative == "." and burned_flag_columns:
                expected.setdefault("_keywords_", {})[POLICY_KEYWORD] = export_marker(burned_flag_columns)
            if a.nrows() != b.nrows() or names != b.colnames() or not equal(logical_descriptor(expected), logical_descriptor(relocate(b.getdesc(), restored))) or not equal(a.info(), b.info()):
                raise ValueError(f"Exported schema/metadata differs: {relative}")
            if derived:
                for start in range(0, b.nrows(), 1024):
                    count = min(1024, b.nrows()-start)
                    if not equal(b.getcol('DATA', start, count), b.getcol('CORRECTED_DATA', start, count)):
                        raise ValueError('Derived DATA differs from decoded CORRECTED_DATA')
                result.append({'table':'.','column':'DATA','export':'derived_calibrated_copy_matches_CORRECTED_DATA_bits'})
            for column in a.colnames():
                burned = relative == "." and column in burned_flag_columns
                if relative == "." and column in approximate_columns and not burned:
                    # Archive-vs-source measured errors are checked separately.
                    continue
                def matches(original_values, exported_values):
                    if burned:
                        return (original_values.shape == exported_values.shape
                                and original_values.dtype == exported_values.dtype == np.dtype(bool)
                                and not exported_values.any())
                    return equal(original_values, exported_values)

                label = "flag reconstruction" if burned else "exact column"
                for start in range(0, a.nrows(), 1024):
                    stop = min(a.nrows(), start + 1024)
                    try:
                        original_values = np.asarray(a.getcol(column, startrow=start, nrow=stop - start))
                        exported_values = np.asarray(b.getcol(column, startrow=start, nrow=stop - start))
                    except RuntimeError:
                        for row in range(start, stop):
                            defined_a, defined_b = a.iscelldefined(column, row), b.iscelldefined(column, row)
                            if defined_a != defined_b or defined_a and not matches(np.asarray(a.getcell(column, row)), np.asarray(b.getcell(column, row))):
                                raise ValueError(f"Exported {label} differs: {relative}/{column}, row {row}")
                    else:
                        if not matches(original_values, exported_values):
                            raise ValueError(f"Exported {label} differs: {relative}/{column}, rows {start}:{stop}")
                result.append({"table": relative, "column": column,
                               "export": "all_false_shapes_and_presence_preserved" if burned else "bitwise_exact"})
    return result
