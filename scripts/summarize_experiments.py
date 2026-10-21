"""Collect measured archive points and render a small rate/distortion figure."""

import argparse
import json
from pathlib import Path

import numpy as np


def point(name, encoding, configuration=None):
    metrics = encoding.get("visibility_quality", {}).get("DATA", {}).get("by_product", {})
    count = sum(v["valid_samples"] for v in metrics.values())
    mse = sum(v["component_mse"] * v["valid_samples"] for v in metrics.values()) / count if count else 0.0
    histogram = {}
    for product in encoding.get("visibility_block_statistics", {}).get("DATA", {}).get("products", []):
        for key, value in product["component_bit_width_histogram"].items():
            histogram[key] = histogram.get(key, 0) + value
    return {"name": name, "archive_bytes": encoding["archive_bytes"],
            "whole_archive_ratio": encoding["whole_archive_ratio"],
            "mean_component_mse": mse,
            "max_complex_error": max((v["max_complex_error"] for v in metrics.values()), default=0.0),
            "block_histogram": histogram, "configuration": configuration,
            "stored_bytes_by_category": encoding["stored_bytes_by_category"],
            "by_product": metrics}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--lossless", type=Path, required=True)
    parser.add_argument("--adaptive", type=Path, required=True)
    parser.add_argument("--adaptive-archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plot", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.plot.exists():
        raise ValueError("Summary outputs must be new")
    benchmark = json.loads(args.benchmark.read_text())
    lossless, adaptive = json.loads(args.lossless.read_text()), json.loads(args.adaptive.read_text())
    # Early encode reports did not yet expose the already-stored block histogram.
    if "visibility_block_statistics" not in adaptive:
        from varquanms.container import Reader
        with Reader(args.adaptive_archive, verify=False) as reader:
            spec = reader.manifest["tables"][0]["columns"]["DATA"]
            adaptive["visibility_block_statistics"] = {"DATA": {"products": spec["products"]}}
    points = [point("Lossless with deduplication", lossless), point("Grouped original adaptive", adaptive)]
    for experiment in benchmark["experiments"]:
        points.append(point(experiment["name"], experiment["encoding"], experiment["configuration"]))
    result = {"schema": "varquanms.initial-rate-distortion/1", "points": points,
              "cache_methodology": benchmark["cache_methodology"],
              "interpretation": "Different points can have different measured errors under the same ceilings. This limited sample is not an equal-error proof, a block-size sweep, a thermal-noise guarantee or scientific acceptance.",
              "timing_comparison": "Early standalone run and later benchmark use different prototype revisions; do not compare their encode times as a codec-performance benchmark."}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    for index, entry in enumerate(points):
        ax.scatter(entry["archive_bytes"] / 1e6, entry["mean_component_mse"], s=60,
                   marker=("o", "s", "^", "D", "v", "P", "X")[index % 7], label=entry["name"])
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.set_ylim(bottom=-.02 * max(p["mean_component_mse"] for p in points))
    ax.set_xlabel("Whole archive (MB, including metadata, exceptions and mappings)")
    ax.set_ylabel("Mean per-component MSE (stored DATA units squared)")
    ax.set_title("Initial rate–distortion points: one 192-channel integration")
    ax.grid(alpha=.25)
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1, 1))
    fig.savefig(args.plot, dpi=160)
    print(json.dumps({"points": len(points), "summary": str(args.output), "plot": str(args.plot)}))


if __name__ == "__main__":
    main()
