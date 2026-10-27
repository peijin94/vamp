"""Image verified benchmark outputs as they finish; write measured FITS reports."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

from vamp.benchmark import load_matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--benchmarks", type=Path, nargs="+", required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--beam-from", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="existing imaging root; each experiment subdirectory must be new")
    parser.add_argument("--reports", type=Path, required=True)
    parser.add_argument("--experiments", nargs="+")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--report-prefix", default="wsclean-")
    parser.add_argument("--wait-seconds", type=float, default=3600)
    args = parser.parse_args()
    _, configs = load_matrix(args.matrix)
    available = {name: config.to_dict() for name, config in configs}
    names = args.experiments or list(available)
    if len(names) != len(set(names)) or set(names) - set(available):
        raise ValueError("Unknown or duplicate imaging experiments")
    if args.threads < 1 or not args.report_prefix or Path(args.report_prefix).name != args.report_prefix:
        raise ValueError("Invalid thread count or report prefix")
    if not args.output.is_dir() or not args.reports.is_dir():
        raise ValueError("Imaging root and report directory must exist")
    # Preflight existing outputs rather than silently reusing possibly unrelated data.
    for name in names:
        if (args.output / name).exists() or (args.reports / f"{args.report_prefix}{name}.json").exists():
            raise ValueError(f"Imaging output already exists: {name}")
    scripts = Path(__file__).resolve().parent
    pending = list(names)
    deadline = time.monotonic() + args.wait_seconds
    while pending:
        ready = None
        for name in pending:
            found = [root / name for root in args.benchmarks if (root / name / "complete.json").exists()]
            if len(found) > 1:
                raise ValueError(f"Ambiguous experiment identity: {name}")
            if found:
                directory = found[0]
                record = json.loads((directory / "complete.json").read_text())
                if record.get("verification_status") != "passed" or record["configuration"] != available[name]:
                    raise ValueError(f"Benchmark verification/configuration mismatch: {name}")
                verification = json.loads((directory / "verification.json").read_text())
                if verification.get("status") != "passed" or not (directory / "restored.ms").is_dir():
                    raise ValueError(f"Verified export required: {name}")
                ready = name, directory
                break
        if ready is None:
            if time.monotonic() >= deadline:
                raise RuntimeError(f"Timed out waiting for benchmark completion: {pending}")
            time.sleep(5)
            continue
        name, directory = ready
        print(json.dumps({"experiment": name, "stage": "wsclean"}), flush=True)
        subprocess.run([sys.executable, str(scripts / "run_wsclean.py"), str(directory / "restored.ms"),
                        "--output", str(args.output / name), "--beam-from", str(args.beam_from), "--threads", str(args.threads)], check=True)
        print(json.dumps({"experiment": name, "stage": "compare_images"}), flush=True)
        subprocess.run([sys.executable, str(scripts / "compare_wsclean.py"), str(args.reference), str(args.output / name),
                        "--report", str(args.reports / f"{args.report_prefix}{name}.json"),
                        "--plot", str(args.reports / f"{args.report_prefix}{name}.png")], check=True)
        pending.remove(name)
    print(json.dumps({"status": "complete", "imaged_experiments": names}), flush=True)


if __name__ == "__main__":
    main()
