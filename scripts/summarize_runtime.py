"""Summarize fresh-process encode timings without double-counting nested timers."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import statistics

from varquanms.cli import _write_new
from varquanms.container import Reader


def distribution(values):
    return {"count": len(values), "median": statistics.median(values),
            "minimum": min(values), "maximum": max(values),
            "mean": statistics.mean(values),
            "sample_stddev": statistics.stdev(values) if len(values) > 1 else None}


def stages(encoding):
    """Exclusive measured subregions plus a labelled residual; codec_total is nested."""
    outer = encoding["timings_seconds"]
    detail = encoding.get("visibility_block_statistics", {}).get("DATA", {}).get("timings_seconds", {})
    result = {"inspection_preflight": outer["inspection_preflight"]}
    for key in ("input_reading", "amplitude_statistics", "baseline_sorting_and_grouping",
                "memory_gathering", "quantization_and_scale_estimation",
                "packing_and_secondary_compression", "candidate_reconstruction_and_error_checks",
                "payload_writing", "quality_accounting", "parallel_block_pipeline_wall", "native_candidate_search"):
        result[key] = detail.get(key, 0.0)
    result["integrity_finalization_and_publication"] = encoding["elapsed_seconds"] - outer["before_integrity_finalization"]
    if any(value < 0 for value in result.values()):
        raise ValueError("Negative/inconsistent timing region")
    remainder = encoding["elapsed_seconds"] - sum(result.values())
    if remainder < -1e-6:
        raise ValueError("Overlapping/inconsistent timers; cannot construct exclusive stages")
    result["other_exact_columns_and_uninstrumented_overhead"] = max(0.0, remainder)
    return result


def summarize(benchmark):
    if benchmark.get("status") != "complete":
        raise ValueError("Benchmark is incomplete")
    records = []
    for entry in benchmark["experiments"]:
        if entry["verification_status"] != "passed":
            raise ValueError("Every timed archive must pass source verification")
        e = entry["encoding"]
        records.append({"name": entry["name"], "mode": entry["configuration"]["mode"],
                        "configuration": entry["configuration"], "archive_bytes": e["archive_bytes"],
                        "source_physical_file_bytes": e["source_physical_file_bytes"],
                        "selected_visibility_logical_bytes": e["selected_visibility_logical_bytes"],
                        "encode_process_seconds": entry["process_timings_seconds"]["encode_process"],
                        "encode_function_seconds": e["elapsed_seconds"],
                        "verification_process_seconds": entry["process_timings_seconds"]["verification_process"],
                        "peak_process_rss_bytes": e["peak_process_rss_bytes"],
                        "exclusive_encode_stages_seconds": stages(e),
                        "worker_count": e.get("visibility_block_statistics", {}).get("DATA", {}).get("worker_count", 1),
                        "worker_timing_sums_seconds": e.get("visibility_block_statistics", {}).get("DATA", {}).get("worker_timings_seconds", {}),
                        "source_MB_per_second": e["source_physical_file_bytes"] / 1e6 / entry["process_timings_seconds"]["encode_process"],
                        "visibility_Mcomplex_samples_per_second": e["selected_visibility_logical_bytes"] / 8e6 / entry["process_timings_seconds"]["encode_process"],
                        "verification_status": "passed"})
    adaptive = [r for r in records if r["mode"] == "adaptive"]
    if not adaptive or any(r["configuration"] != adaptive[0]["configuration"] for r in adaptive):
        raise ValueError("Adaptive repeats must have identical configurations")
    aggregate = {key: distribution([r[key] for r in adaptive]) for key in (
        "encode_process_seconds", "encode_function_seconds", "verification_process_seconds",
        "archive_bytes", "peak_process_rss_bytes", "source_MB_per_second",
        "visibility_Mcomplex_samples_per_second")}
    aggregate["exclusive_encode_stages_seconds"] = {
        key: distribution([r["exclusive_encode_stages_seconds"][key] for r in adaptive])
        for key in adaptive[0]["exclusive_encode_stages_seconds"]}
    return {"schema": "varquanms.paper-runtime/1", "status": "complete",
            "created_utc": datetime.now(timezone.utc).isoformat(), "source": benchmark["source"],
            "cache_methodology": benchmark["cache_methodology"], "methodology": benchmark["notes"],
            "records": records, "adaptive_summary": aggregate,
            "timer_scope": {
                "encode_process": "Fresh Python process: startup + inspect + encode + integrity/final publication + report serialization; no verification/export/imaging",
                "encode_function": "encode() entry through complete archive publication, including input inspection; excludes CLI startup/report serialization",
                "stages": "Exclusive measured regions plus residual. codec_total and before_integrity_finalization are nested/cumulative and are NOT summed. Stage medians need not sum to median total.",
                "writing": "payload_writing covers visibility payload writes only. quality_accounting includes quality-array writes. Exact-column/index/statistic writes and miscellaneous setup are included in the residual; integrity stage includes flushing, checksums, manifest and fsync/publication.",
                "measurement_limits": "Shared host; no cache eviction, CPU affinity, resource isolation or cold-cache claim. Packing and secondary compression are a combined timer."}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("benchmark", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--source-inventory", required=True, type=Path)
    parser.add_argument("--activity-note", default="", help="record other known work on the shared host during timing")
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("Runtime report must be new")
    benchmark = json.loads(args.benchmark.read_text())
    result = summarize(benchmark)
    environment = []
    reference_inventory = None
    repeated_archives = 0
    for entry in benchmark["experiments"]:
        with Reader(entry["encoding"]["archive"], verify=False) as reader:
            environment.append({"name": entry["name"], "environment": reader.manifest["environment"]})
            if entry["configuration"]["mode"] == "adaptive":
                inventory = {name: (item["shape"], item["dtype"], item["sha256"]) for name, item in reader.entries.items()}
                if reference_inventory is not None and inventory != reference_inventory:
                    raise ValueError("Adaptive repeat dataset contents are not identical")
                reference_inventory = inventory
                repeated_archives += 1
    result["run_environments"] = environment
    result["adaptive_repeat_dataset_identity"] = {"status": "passed", "archives_compared": repeated_archives,
        "scope": "All inventoried dataset paths, shapes, dtypes and SHA-256 values match; manifest timestamps/timings/environment are excluded"}
    result["shared_host_activity_note"] = args.activity_note
    cpu = Path("/proc/cpuinfo").read_text()
    result["hardware"] = {"cpu_model": next(line.split(":", 1)[1].strip() for line in cpu.splitlines() if line.startswith("model name")),
                          "memory_total_kib": int(Path("/proc/meminfo").read_text().splitlines()[0].split()[1])}
    inventory = json.loads(args.source_inventory.read_text())
    parent = Path(inventory["extracted_parent"])
    for entry in inventory["files"]:
        path = parent / entry["path"]
        with path.open("rb") as stream:
            digest = hashlib.sha256()
            for block in iter(lambda: stream.read(4 * 1024**2), b""):
                digest.update(block)
        if digest.hexdigest() != entry["sha256"] or path.stat().st_size != entry["bytes"]:
            raise ValueError(f"Source changed: {path}")
    source = Path(benchmark["source"])
    actual = {str(p.relative_to(parent.resolve())) for p in source.rglob("*") if p.is_file()}
    if actual != {e["path"] for e in inventory["files"]}:
        raise ValueError("Source file inventory changed")
    result["source_unchanged"] = {"sha256_matched_files": len(inventory["files"]), "reference_inventory": str(args.source_inventory)}
    result["benchmark_report"] = str(args.benchmark)
    _write_new(args.output, json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"report": str(args.output), "encode_process_seconds": result["adaptive_summary"]["encode_process_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
