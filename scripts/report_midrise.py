"""Summarize matched unsigned-code experiments and draw reproducible figures."""
import argparse
from collections import Counter
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from casacore.tables import table

from run_wsclean import snapshot
from vamp.container import Reader
from vamp.visibility import quality

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams["svg.hashsalt"] = "varquanms-midrise-v1"


def save(fig, name):
    for ext in ("png", "svg", "pdf"):
        path = ROOT / "docs/paper/figs" / f"{name}.{ext}"
        metadata = {"Date": None} if ext == "svg" else {"CreationDate": None, "ModDate": None} if ext == "pdf" else {}
        fig.savefig(path, dpi=180, bbox_inches="tight", metadata=metadata)
        if ext == "svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    plt.close(fig)


def figures(report):
    fig, ax = plt.subplots(figsize=(10, 3.5))
    for y, denominator, label, color in ((1, 7, "Symmetric: 15 levels, exact zero", "#226dac"),
                                         (0, 7.5, "Midrise: 16 levels, no exact zero", "#b86a20")):
        levels = np.arange(-denominator, denominator+1) / denominator
        ax.hlines(y, -1, 1, colors="#bbbbbb")
        ax.scatter(levels, np.full(levels.shape, y), c=color, s=55, zorder=3)
        ax.text(-1, y+.22, label, color=color, fontsize=11)
    ax.axvline(0, color="#666666", ls=":")
    ax.annotate("ordinary 0 → +Δ/2", xy=(1/15, 0), xytext=(.28, -.45),
                arrowprops={"arrowstyle": "->"}, fontsize=10)
    ax.set(xlim=(-1.08, 1.08), ylim=(-.65, 1.55), yticks=[], xlabel="Reconstruction / common cutoff t", title="4-bit unsigned codebooks: one additional level")
    save(fig, "uint-levels")
    rows = report["codecs"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.1), constrained_layout=True)
    labels = ["Symmetric", "Midrise"]
    colors = ["#226dac", "#b86a20"]
    axes[0].bar(labels, [r["encoding"]["archive_bytes"]/1e6 for r in rows], color=colors)
    axes[0].set(ylabel="Whole archive (decimal MB)", title="Same configured error ceilings")
    axes[1].bar(labels, [r["component_mse"] for r in rows], color=colors)
    axes[1].set(ylabel="Measured component MSE", title="Achieved error is compared")
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    for i, row in enumerate(rows):
        values = row["encode_process_seconds"]
        axes[2].scatter(np.full(len(values), i), values, color=colors[i], s=35)
        axes[2].hlines(statistics.median(values), i-.2, i+.2, color=colors[i], lw=3)
    axes[2].axhline(15, color="#777777", ls="--", label="15 s target")
    axes[2].set(xticks=[0,1], xticklabels=labels, ylabel="Fresh process encode (seconds)", title="Three repeats; 24 workers", ylim=(0, max(16, max(max(r["encode_process_seconds"]) for r in rows)*1.1)))
    axes[2].legend(fontsize=9)
    save(fig, "uint-comparison")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT/"results/midrise-20261023")
    parser.add_argument("--figures-only", action="store_true")
    args = parser.parse_args()
    destination = ROOT/"reports/uint-midrise-evaluation.json"
    if args.figures_only:
        figures(json.loads(destination.read_text())); return
    benchmark = json.loads((args.root/"benchmark/results.json").read_text())
    if benchmark["status"] != "complete":
        raise ValueError("Incomplete benchmark")
    report = {"schema": "vamp.uint-midrise-study/1", "date_utc": "2026-10-23",
              "source": benchmark["source"], "methodology": benchmark["notes"],
              "cache_methodology": benchmark["cache_methodology"], "codecs": [],
              "formula": {"levels": "N=2**b; C=(N-1)/2; delta=ceil_float32(t/C)",
                          "encoder": "q=rint_even(float64(x)/float64(delta)+C), 0 <= q < N",
                          "decoder": "float32(float32(q)-float32(C))*stored_float32_delta",
                          "positive_zero_nonzero_step": "+delta/2 under nearest-even",
                          "old_to_new_step_ratio": "(N-2)/(N-1)"}}
    for name in ("uniform-symmetric-v1", "uniform-midrise-v1"):
        entries = [e for e in benchmark["experiments"] if e["configuration"]["codebook"] == name]
        if len(entries) != 3 or any(e["verification_status"] != "passed" for e in entries):
            raise ValueError("Require three verified runs per codec")
        encoding = entries[0]["encoding"]
        q = encoding["visibility_quality"]["DATA"]["by_product"]
        samples = sum(x["valid_samples"] for x in q.values())
        hist, sample_hist, accounting = Counter(), Counter(), Counter()
        for p in encoding["visibility_block_statistics"]["DATA"]["products"]:
            hist.update(p["component_bit_width_histogram"])
            sample_hist.update(p["complex_sample_bit_width_histogram"])
            accounting.update(p["logical_bytes_before_secondary_compression"])
        checksums = []
        for e in entries:
            with Reader(e["encoding"]["archive"], verify=False) as reader:
                checksums.append({n: d["sha256"] for n, d in reader.entries.items()})
                environment = reader.manifest["environment"]
        if any(c != checksums[0] for c in checksums[1:]):
            raise ValueError("Repeated payloads differ")
        report["codecs"].append({"codebook": name, "configuration": entries[0]["configuration"],
             "encoding": encoding, "environment": environment,
             "component_mse": sum(x["component_mse"]*x["valid_samples"] for x in q.values())/samples,
             "block_histogram": dict(hist), "sample_histogram": dict(sample_hist),
             "logical_byte_accounting": dict(accounting),
             "encode_process_seconds": [e["process_timings_seconds"]["encode_process"] for e in entries],
             "verification_process_seconds": [e["process_timings_seconds"]["verification_process"] for e in entries],
             "archive_bytes_repeats": [e["encoding"]["archive_bytes"] for e in entries],
             "payloads_identical_across_repeats": True, "verification": "passed"})
    source = Path(benchmark["source"])
    with table(str(source), readonly=True, lockoptions="usernoread", ack=False) as t:
        original = t.getcol("DATA"); flags = t.getcol("FLAG") | t.getcol("FLAG_ROW")[:, None, None]
        auto = t.getcol("ANTENNA1") == t.getcol("ANTENNA2")
    with table(str(args.root/"restored.ms"), readonly=True, lockoptions="usernoread", ack=False) as t:
        restored = t.getcol("DATA")
    report["midrise_by_baseline_class"] = {}
    for name, mask in (("autocorrelation", auto), ("cross_correlation", ~auto)):
        summary, _ = quality(original[mask], restored[mask], flags[mask], list(q))
        report["midrise_by_baseline_class"][name] = summary
    report["export_verification"] = json.loads((args.root/"export-verification.json").read_text())
    before = json.loads((args.root/"source-before.json").read_text()); after = snapshot(source)
    if before != after:
        raise ValueError("Source file hashes changed")
    report.update(source_files_unchanged=True, source_file_sha256=before)
    report["imaging"] = json.loads((ROOT/"reports/wsclean-midrise.json").read_text())
    run = json.loads((args.root/"wsclean/complete.json").read_text())
    report["imaging_run"] = {k: run[k] for k in ("status", "container_id", "wsclean_version", "parameters", "elapsed_seconds", "source_unchanged")}
    report["limits"] = "One integration, engineering bounds with no validated noise model; Stokes I imaging only, fixed beam, 10,000-iteration CLEAN cap. No calibrated solar peak/weak-polarization, calibration-solution, narrow-feature or temporal acceptance. Timings exclude verification, export and imaging."
    destination.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    figures(report)


if __name__ == "__main__":
    main()
