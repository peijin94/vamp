"""Reproducible sequential experiments, isolated processes and explicit checks."""

import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from .config import Config


def load_matrix(path):
    matrix = json.loads(Path(path).read_text())
    if not isinstance(matrix, dict) or set(matrix) - {"schema", "base", "experiments", "export_ms", "cache_methodology", "notes"}:
        raise ValueError("Unknown benchmark matrix field")
    if matrix.get("schema") not in ("vamp.benchmark/1", "varquanms.benchmark/1") or not isinstance(matrix.get("base", {}), dict):
        raise ValueError("Unsupported benchmark schema/base configuration")
    if type(matrix.get("export_ms", False)) is not bool or not isinstance(matrix.get("cache_methodology"), str):
        raise ValueError("Benchmark requires explicit cache_methodology and boolean export_ms")
    experiments = matrix.get("experiments")
    if not isinstance(experiments, list) or not experiments:
        raise ValueError("Benchmark experiments must be a nonempty list")
    names, normalized = set(), []
    for item in experiments:
        if not isinstance(item, dict) or set(item) != {"name", "config"} or not isinstance(item["config"], dict):
            raise ValueError("Each experiment requires only name and config")
        name = item["name"]
        if not isinstance(name, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", name) or name in names:
            raise ValueError("Experiment names must be unique safe path components")
        names.add(name)
        config = Config.load(None, **{**matrix.get("base", {}), **item["config"]})
        normalized.append((name, config))
    return matrix, normalized


def run(source, matrix_path, output, selected=None):
    from .cli import _write_new
    matrix, experiments = load_matrix(matrix_path)
    if selected is not None:
        if len(set(selected)) != len(selected) or set(selected) - {name for name, _ in experiments}:
            raise ValueError("Unknown or duplicate selected benchmark experiment")
        experiments = [(name, config) for name, config in experiments if name in selected]
    source, output = Path(source).resolve(), Path(output).resolve()
    if output.exists() or output.is_relative_to(source) or not output.parent.is_dir():
        raise ValueError("Benchmark output must be a new directory outside the source")
    output.mkdir()
    records = []
    started = time.perf_counter()
    _write_new(output / "matrix.json", json.dumps(matrix, indent=2) + "\n")

    def invoke(name, stage, arguments):
        print(json.dumps({"experiment": name, "stage": stage}), flush=True)
        stamp = time.perf_counter()
        result = subprocess.run([sys.executable, "-m", "vamp", *arguments],
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _write_new(output / name / f"{stage}.stdout.log", result.stdout)
        _write_new(output / name / f"{stage}.stderr.log", result.stderr)
        if result.returncode:
            raise RuntimeError(f"Benchmark {name}/{stage} failed: {result.stderr.strip()}")
        return time.perf_counter() - stamp

    for name, config in experiments:
        directory = output / name
        directory.mkdir()
        _write_new(directory / "config.json", json.dumps(config.to_dict(), indent=2) + "\n")
        archive = directory / "archive.h5"
        timing = {}
        timing["encode_process"] = invoke(name, "encode", ["encode", str(source), str(archive), "--config", str(directory / "config.json"), "--report", str(directory / "encode.json")])
        verification = ["verify", str(archive), "--source", str(source), "--report", str(directory / "verification.json")]
        if matrix.get("export_ms", False):
            timing["decode_and_export_process"] = invoke(name, "decode", ["decode", str(archive), str(directory / "restored.ms"), "--report", str(directory / "decode.json")])
            verification += ["--restored", str(directory / "restored.ms")]
        timing["verification_process"] = invoke(name, "verify", verification)
        record = {"name": name, "configuration": config.to_dict(),
                  "encoding": json.loads((directory / "encode.json").read_text()),
                  "verification_status": json.loads((directory / "verification.json").read_text())["status"],
                  "process_timings_seconds": timing}
        records.append(record)
        _write_new(directory / "complete.json", json.dumps(record, indent=2) + "\n")
    result = {"schema": "vamp.benchmark-results/1", "status": "complete", "source": str(source),
              "cache_methodology": matrix["cache_methodology"], "notes": matrix.get("notes"),
              "experiments": records, "selected_experiments": [name for name, _ in experiments],
              "elapsed_seconds": time.perf_counter() - started,
              "Dysco": {"status": "not_run", "available_executable": shutil.which("aocompress") or shutil.which("dysco-compress"),
                        "reason": "External Dysco benchmark adapter/settings equivalence are not implemented"},
              "scientific_products": "not tested by this benchmark; numeric bounds are not scientific acceptance",
              "timing_note": "Each experiment uses fresh Python processes; RSS is that process's high-water mark. Caches are not forcibly evicted. Some inner codec stages remain combined in the recorded timer names."}
    _write_new(output / "results.json", json.dumps(result, indent=2, allow_nan=False) + "\n")
    return result
