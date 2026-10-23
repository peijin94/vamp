"""Experimental ordered bundle of independently decodable single-context archives.

Runs of equal row context become partitions. All original rows and all subtables
are retained; repeated subtable bytes are counted, not claimed as deduplicated.
The original MS is not required by the bundle decoder.
"""

from collections import Counter
import json
from pathlib import Path
import tempfile
import time

from casacore.tables import table
import numpy as np

from .archive import decode, discover_tables, encode, read_column, verify
from .config import Config
from .calibrated import OMITTED, materialize_ms_v2_data
from .container import Reader, new_output
from .imaging import file_sha256
from .inspection import ROW_KEYS, _open, inspect_ms
from .metadata import dumps, logical_descriptor, pack, unpack

SCHEMA = "varquanms.ordered-partitions/1"


def context_runs(source):
    """Return contiguous [start, stop) rules using exact context-column bits."""
    with _open(source) as t:
        if not t.nrows():
            raise ValueError("Empty MAIN table")
        names = [c for c in ROW_KEYS if c in t.colnames() and c not in ("ANTENNA1", "ANTENNA2")]
        changed = np.zeros(t.nrows() - 1, dtype=bool)
        for name in names:
            a = np.asarray(t.getcol(name))
            if a.ndim != 1 or a.dtype.kind not in "biuf":
                raise ValueError("Context columns must be numeric scalars")
            bits = np.ascontiguousarray(a).view(np.uint8).reshape(len(a), -1)
            changed |= np.any(bits[1:] != bits[:-1], axis=1)
        bounds = np.r_[0, np.flatnonzero(changed) + 1, t.nrows()]
    if len(bounds) > 4097:
        raise ValueError("More than 4096 contiguous contexts; explicit interleaved maps are not implemented")
    return names, [(int(a), int(b)) for a, b in zip(bounds[:-1], bounds[1:])]


def _signature(manifest):
    """Require identical logical schemas and exact subtable contents in a bundle."""
    entries = {x["path"]: x for x in manifest["inventory"]}
    result = []
    for spec in manifest["tables"]:
        node = {"path": spec["path"], "columns": spec["column_names"],
                "descriptor": pack(logical_descriptor(unpack(spec["descriptor"]))), "info": spec["info"]}
        if spec["path"] != ".":
            node.update(rows=spec["row_count"], encodings=spec["columns"])
        result.append(node)
    def referenced(node):
        if isinstance(node, str):
            return {node} if node in entries else set()
        if isinstance(node, dict):
            return set().union(*(referenced(v) for v in node.values()))
        if isinstance(node, list):
            return set().union(*(referenced(v) for v in node))
        return set()
    paths = set().union(*(referenced(s["columns"]) for s in manifest["tables"] if s["path"] != "."))
    subtable_entries = [{k: entries[path][k] for k in ("path", "shape", "dtype", "sha256")}
                       for path in sorted(paths)]
    return dumps({"tables": result, "subtable_values": subtable_entries})


def encode_partitioned(source, output, config=None):
    """Encode contiguous contexts; publish the bundle only after all checks pass."""
    started = time.perf_counter()
    source = Path(source).resolve()
    output = new_output(output, source)
    config = config or Config()
    config.validate()
    if config.burn_flag:
        raise ValueError("Ordered partition bundles currently require exact flag preservation")
    discover_tables(source)
    names, runs = context_runs(source)
    original_bytes = sum(p.stat().st_size for p in source.rglob("*") if p.is_file())
    staged = Path(tempfile.mkdtemp(prefix=output.name + ".", suffix=".partial", dir=output.parent))
    manifest = {"schema": SCHEMA, "status": "complete", "config": config.to_dict(),
                "source_name": source.name, "source_physical_bytes": original_bytes,
                "row_count": runs[-1][1], "context_columns": names, "partitions": [],
                "mapping": "Original rows [start,stop) equal local rows [0,stop-start); no averaging or row permutation",
                "scope": "Contiguous numeric contexts; each partition must pass the existing archive preflight. All subtables are repeated exactly in every archive."}
    signature = None
    for index, (start, stop) in enumerate(runs):
        copying = time.perf_counter()
        with tempfile.TemporaryDirectory(prefix="partition-ms-", dir=output.parent) as scratch:
            part = Path(scratch) / "input.ms"
            with _open(source) as t, t.selectrows(np.arange(start, stop)) as selected:
                if config.keep_calib_only:
                    with selected.query(columns=','.join(c for c in selected.colnames() if c not in OMITTED)) as projection:
                        projection.copy(str(part), deep=True, valuecopy=True).close()
                else:
                    selected.copy(str(part), deep=True, valuecopy=True).close()
            copy_seconds = time.perf_counter() - copying
            scanning = time.perf_counter()
            inspection = inspect_ms(part, config.columns, config.io_chunk_rows)
            inspection_seconds = time.perf_counter() - scanning
            name = f"partition-{index:04d}.h5"
            summary = encode(part, staged / name, config, inspection)
            summary["archive"] = name  # Relative to the bundle after publication.
            with Reader(staged / name) as reader:
                current = _signature(reader.manifest)
            if signature is not None and current != signature:
                raise ValueError("Partition schemas/subtables differ; bundle cannot reconstruct them safely")
            signature = current
            validation = verify(staged / name, source=part)
            validation_name = f"verification-{index:04d}.json"
            (staged / validation_name).write_text(json.dumps(validation, indent=2, allow_nan=False) + "\n")
        record = {"start": start, "stop": stop, "archive": name, "bytes": (staged / name).stat().st_size,
                  "sha256": file_sha256(staged / name), "copy_seconds": copy_seconds,
                  "inspection_seconds": inspection_seconds, "summary": summary,
                  "verification": validation_name, "verification_sha256": file_sha256(staged / validation_name)}
        manifest["partitions"].append(record)
        print(json.dumps({"partition": index, "rows": stop-start, "bytes": record["bytes"],
                          "encode_seconds_excluding_inspection": summary["elapsed_seconds"]}), flush=True)
    manifest["elapsed_seconds_including_copy_inspection_verification"] = time.perf_counter() - started
    (staged / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    (staged / "manifest.sha256").write_text(file_sha256(staged / "manifest.json") + "\n")
    if output.exists():
        raise FileExistsError(output)
    staged.rename(output)
    return bundle_summary(output, manifest)


def _load(bundle):
    bundle = Path(bundle).resolve()
    if file_sha256(bundle / "manifest.json") != (bundle / "manifest.sha256").read_text().strip():
        raise ValueError("Bundle manifest checksum mismatch")
    manifest = json.loads((bundle / "manifest.json").read_text())
    if manifest.get("schema") != SCHEMA or manifest.get("status") != "complete":
        raise ValueError("Unsupported or incomplete partition bundle")
    parts = manifest["partitions"]
    if not 1 <= len(parts) <= 4096:
        raise ValueError("Invalid partition count")
    cursor, seen, signature = 0, set(), None
    for part in parts:
        start, stop = part["start"], part["stop"]
        if type(start) is not int or type(stop) is not int or start != cursor or stop <= start:
            raise ValueError("Partition row mapping has gaps, overlaps or invalid bounds")
        cursor = stop
        for key, digest_key in (("archive", "sha256"), ("verification", "verification_sha256")):
            name = part[key]
            if not isinstance(name, str) or Path(name).name != name or name in (".", "..") or name in seen:
                raise ValueError("Unsafe or duplicate bundle member")
            seen.add(name)
            path = bundle / name
            if path.is_symlink() or not path.is_file() or file_sha256(path) != part[digest_key]:
                raise ValueError("Bundle member checksum mismatch")
        with Reader(bundle / part["archive"]) as reader:
            if reader.manifest["source"]["main_row_count"] != stop-start or reader.manifest["config"] != manifest["config"]:
                raise ValueError("Partition row count/config mismatch")
            current = _signature(reader.manifest)
            if signature is not None and current != signature:
                raise ValueError("Partition schemas/subtables differ")
            signature = current
    if cursor != manifest["row_count"]:
        raise ValueError("Partition row count mismatch")
    return bundle, manifest


def decode_partitioned(bundle, output, ms_v2_data=None):
    """Restore all rows using only the bundle; integrity is checked before output."""
    started = time.perf_counter()
    output = new_output(output, bundle)
    bundle, manifest = _load(bundle)
    if ms_v2_data is None:
        ms_v2_data = manifest['config'].get('keep_calib_only', False)
    if ms_v2_data and not manifest['config'].get('keep_calib_only', False):
        raise ValueError('--ms-v2-data is only for a calibrated-only bundle')
    staged = Path(tempfile.mkdtemp(prefix=output.name + ".", suffix=".partial", dir=output.parent))
    restored = staged / "restored.ms"
    decode(bundle / manifest["partitions"][0]["archive"], restored, ms_v2_data=False)
    with table(str(restored), readonly=False, ack=False) as t:
        t.addrows(manifest["row_count"] - t.nrows())
        for part in manifest["partitions"][1:]:
            with Reader(bundle / part["archive"]) as reader:
                main = next(s for s in reader.manifest["tables"] if s["path"] == ".")
                for name, column in main["columns"].items():
                    quantized = read_column(reader, column, 0, main["row_count"]) if column["encoding"] == "quantized" else None
                    for start in range(0, main["row_count"], 1024):
                        stop = min(start+1024, main["row_count"])
                        values = quantized[start:stop] if quantized is not None else read_column(reader, column, start, stop)
                        target = part["start"] + start
                        if isinstance(values, np.ndarray) and values.size:
                            t.putcol(name, values, startrow=target, nrow=len(values))
                        else:
                            for offset, value in enumerate(values):
                                if value is not None:
                                    cell = value.item() if isinstance(value, np.ndarray) and value.ndim == 0 else value
                                    t.putcell(name, target+offset, cell)
    if ms_v2_data:
        materialize_ms_v2_data(restored)
    with _open(restored) as t:
        if t.nrows() != manifest["row_count"]:
            raise ValueError("Reassembled MS row count mismatch")
    if output.exists():
        raise FileExistsError(output)
    restored.rename(output)
    staged.rmdir()
    return {"status": "exported", "output": str(output), "row_count": manifest["row_count"],
            "elapsed_seconds": time.perf_counter()-started, "source_required": False,
            "derived_calibrated_DATA": ms_v2_data}


def bundle_summary(bundle, manifest=None):
    bundle = Path(bundle)
    manifest = manifest or json.loads((bundle / "manifest.json").read_text())
    categories = Counter()
    for part in manifest["partitions"]:
        categories.update(part["summary"]["stored_bytes_by_category"])
    total = sum(p.stat().st_size for p in bundle.iterdir() if p.is_file())
    categories["bundle_manifest_and_verification"] = total-sum(categories.values())
    return {"bundle": str(bundle), "partitions": len(manifest["partitions"]), "row_count": manifest["row_count"],
            "archive_bytes": total, "source_physical_bytes": manifest["source_physical_bytes"],
            "whole_archive_ratio": manifest["source_physical_bytes"]/total, "stored_bytes_by_category": dict(categories),
            "sum_encode_seconds_excluding_inspection": sum(p["summary"]["elapsed_seconds"] for p in manifest["partitions"]),
            "sum_inspection_seconds": sum(p["inspection_seconds"] for p in manifest["partitions"]),
            "sum_partition_copy_seconds": sum(p["copy_seconds"] for p in manifest["partitions"]),
            "elapsed_seconds_including_copy_inspection_verification": manifest["elapsed_seconds_including_copy_inspection_verification"]}
