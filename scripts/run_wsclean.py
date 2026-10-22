"""Run pinned-container WSClean on a private MS copy with reproducible settings."""

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from varquanms.cli import _write_new
from varquanms.imaging import file_sha256


def snapshot(path):
    result = {}
    for item in sorted(path.rglob("*")):
        if item.is_symlink():
            raise ValueError("Imaging input must not contain symlinks")
        if item.is_file():
            result[str(item.relative_to(path))] = file_sha256(item)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ms", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--container", default="docker.io/astronrd/linc:5.2.1")
    parser.add_argument("--size", type=int, default=2048)
    parser.add_argument("--scale", default="4amin")
    parser.add_argument("--niter", type=int, default=10000)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--polarizations", default="I", choices=("I", "I,Q,U,V"))
    parser.add_argument("--beam-from", type=Path, help="use this reference FITS restoring beam for every product")
    args = parser.parse_args()
    source, output = args.ms.resolve(), args.output.resolve()
    if not source.is_dir() or output.exists() or output.is_relative_to(source) or not output.parent.is_dir():
        raise ValueError("Use an existing MS and a new output directory outside it")
    if args.size < 16 or args.niter < 1 or args.threads < 1:
        raise ValueError("Invalid size, iteration count or thread count")
    info = json.loads(subprocess.check_output(["docker", "image", "inspect", args.container], text=True))[0]
    image_id = info["Id"]
    docker = ["docker", "run", "--rm", "--network", "none", "--user", f"{os.getuid()}:{os.getgid()}"]
    version = subprocess.check_output(docker + [image_id, "wsclean", "-version"], text=True)
    output.mkdir()
    started = time.perf_counter()
    before = snapshot(source)
    working = output / "working.ms"
    shutil.copytree(source, working)
    (output / "tmp").mkdir()
    options = ["-j", str(args.threads), "-abs-mem", "8", "-no-update-model-required",
               "-temp-dir", "/work/tmp", "-name", "/work/image", "-data-column", "DATA",
               "-size", str(args.size), str(args.size), "-scale", args.scale,
               "-weight", "briggs", "0", "-pol", args.polarizations,
               "-multiscale", "-mgain", "0.8", "-niter", str(args.niter),
               "-auto-mask", "5", "-auto-threshold", "1", "/work/working.ms"]
    beam = None
    if args.beam_from:
        from astropy.io import fits
        header = fits.getheader(args.beam_from)
        beam = {"major_arcsec": float(header["BMAJ"]) * 3600, "minor_arcsec": float(header["BMIN"]) * 3600,
                "position_angle_deg": float(header["BPA"])}
        options[-1:-1] = ["-beam-shape", str(beam["major_arcsec"]), str(beam["minor_arcsec"]), str(beam["position_angle_deg"])]
    command = docker + ["-e", "OPENBLAS_NUM_THREADS=1", "-v", f"{output}:/work", "-w", "/work", image_id, "wsclean", *options]
    record = {"schema": "varquanms.wsclean/1", "status": "running", "source_ms": str(source),
              "container": args.container, "container_id": image_id, "container_digests": info.get("RepoDigests", []),
              "wsclean_version": version.strip(), "command": command, "options": options,
              "parameters": {"size": args.size, "scale": args.scale, "polarizations": args.polarizations,
                             "weighting": "briggs 0", "multiscale": True, "mgain": .8, "niter": args.niter,
                             "auto_mask": 5, "auto_threshold": 1, "threads": args.threads, "memory_gb": 8, "restoring_beam": beam},
              "source_sha256": before,
              "interpretation": "Identical imaging settings are an experiment, not approved scientific-error tolerances. Auto-masks are separately estimated in each image; CLEAN path changes are part of the measured effect."}
    _write_new(output / "invocation.json", json.dumps(record, indent=2) + "\n")
    with (output / "wsclean.log").open("x") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    record.update(status="complete" if result.returncode == 0 else "failed", returncode=result.returncode,
                  elapsed_seconds=time.perf_counter() - started, source_unchanged=snapshot(source) == before,
                  fits_files=[p.name for p in sorted(output.glob("*.fits"))])
    if not record["source_unchanged"] or not record["fits_files"]:
        record["status"] = "failed"
    _write_new(output / "complete.json", json.dumps(record, indent=2) + "\n")
    print(json.dumps({k: record[k] for k in ("status", "returncode", "elapsed_seconds", "source_unchanged", "fits_files")}, indent=2))
    return 0 if record["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
