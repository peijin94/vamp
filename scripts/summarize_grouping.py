"""Collect complete grouping experiments, numerical errors and WSClean outcomes."""

import argparse
from collections import Counter
import json
import os
from pathlib import Path

from vamp.benchmark import load_matrix
from vamp.cli import _write_new


def point(record, imaging):
    encoding = record["encoding"]
    metrics = encoding["visibility_quality"]["DATA"]["by_product"]
    count = sum(m["valid_samples"] for m in metrics.values())
    histogram, samples, logical = Counter(), Counter(), Counter()
    group_sizes = {}
    for p, product in enumerate(encoding["visibility_block_statistics"]["DATA"]["products"]):
        histogram.update(product["component_bit_width_histogram"])
        samples.update(product["complex_sample_bit_width_histogram"])
        logical.update(product["logical_bytes_before_secondary_compression"])
        group_sizes[str(p)] = product["grouping"]
    total = sum(samples.values())
    image_metrics = {}
    for name, product in imaging["products"].items():
        image_metrics[name] = {"regions": product["regions"], "peak_shift_pixels": product["peak_shift_pixels"],
                               "relative_change_at_reference_peak": product["relative_change_at_reference_peak"],
                               "aperture": product["aperture"], "clean": product["clean"]}
    return {"name": record["name"], "configuration": record["configuration"], "verification_status": record["verification_status"],
            "archive_bytes": encoding["archive_bytes"], "whole_archive_ratio": encoding["whole_archive_ratio"],
            "visibility_ratio_including_maps_statistics": encoding["visibility_dataset_ratio_including_maps_statistics"],
            "stored_bytes_by_category": encoding["stored_bytes_by_category"], "logical_bytes_before_secondary_compression": dict(logical),
            "component_mse": sum(m["component_mse"] * m["valid_samples"] for m in metrics.values()) / count if count else 0,
            "max_complex_error": max(m["max_complex_error"] for m in metrics.values()), "by_product": metrics,
            "block_histogram": dict(histogram), "sample_histogram": dict(samples),
            "sample_fraction_by_component_bits": {key: n / total for key, n in samples.items()}, "grouping_by_product": group_sizes,
            "process_timings_seconds": record["process_timings_seconds"], "peak_process_rss_bytes": encoding["peak_process_rss_bytes"],
            "imaging": image_metrics}


def dominance_vector(entry):
    vector = [entry["archive_bytes"]]
    for name in sorted(entry["by_product"]):
        value = entry["by_product"][name]
        vector.extend(value[key] for key in ("component_mse", "max_complex_error", "complex_bias_abs", "max_baseline_component_mse"))
    vector.extend(entry["imaging"][name]["regions"]["celestial_valid"]["difference_rms"]
                  for name in ("image-dirty.fits", "image-image.fits"))
    return vector


def summarize(matrix, roots, reports, repeat, image_prefix="wsclean-"):
    _, configurations = load_matrix(matrix)
    repeat = json.loads(repeat.read_text())
    points = []
    for name, config in configurations:
        candidates = [root / name / "complete.json" for root in roots if (root / name / "complete.json").exists()]
        if len(candidates) != 1:
            raise ValueError(f"Expected one completed benchmark for {name}, found {len(candidates)}")
        record = json.loads(candidates[0].read_text())
        if record["verification_status"] != "passed" or record["configuration"] != config.to_dict():
            raise ValueError(f"Verification/configuration mismatch for {name}")
        image = json.loads((reports / f"{image_prefix}{name}.json").read_text())
        if image["status"] != "measured":
            raise ValueError(f"Missing measured image comparison: {name}")
        if any(image[key] != repeat[key] for key in ("parameters", "container_id", "reference")):
            raise ValueError(f"Imaging controls differ from the repeatability reference: {name}")
        entry = point(record, image)
        entry["wsclean_run"] = json.loads((Path(image["candidate"]) / "complete.json").read_text())
        from vamp.container import Reader
        with Reader(record["encoding"]["archive"], verify=False) as reader:
            entry["archive_environment"] = reader.manifest["environment"]
        points.append(entry)
    vectors = {p["name"]: dominance_vector(p) for p in points}
    for p in points:
        vector = vectors[p["name"]]
        p["dominated_by"] = [name for name, other in vectors.items() if name != p["name"] and
                             all(x <= y for x, y in zip(other, vector)) and any(x < y for x, y in zip(other, vector))]
    return {"schema": "vamp.grouping-results/1", "status": "complete", "points": points,
            "reference_repeatability": repeat,
            "reference_wsclean_run": json.loads((Path(repeat["reference"]) / "complete.json").read_text()),
            "dominance_definition": "No worse in whole bytes, every product's MSE/maximum/bias/max-baseline-MSE, and dirty/restored valid-sky RMS; strictly better in at least one. This finite-sample relation is not a scientific guarantee or equal-error interpolation.",
            "interpretation": "Common engineering ceilings are not equal measured error. CLEAN differences include nonlinear path changes and the reported reference repeatability floor. Images are Stokes I only.",
            "science_acceptance": "not established; no calibrated noise model or approved image-error thresholds",
            "timing_policy": "Concurrent experiments on a shared host, no isolated throughput claim"}


def plot(result, path):
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/vamp-matplotlib-cache")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    for index, p in enumerate(result["points"]):
        color = f"C{index % 10}"
        marker = ("o", "s", "^", "D", "v", "P", "X", "<", ">", "h", "*", "+")[index % 12]
        x = p["archive_bytes"] / 1e6
        ys = [p["component_mse"], *[100 * p["imaging"][name]["regions"]["celestial_valid"]["relative_rms_difference"]
                                   for name in ("image-dirty.fits", "image-image.fits")]]
        for ax, y in zip(axes, ys):
            ax.scatter(x, y, label=p["name"], color=color, marker=marker, s=55)
    repeat = result["reference_repeatability"]["products"]
    for ax, name in zip(axes[1:], ("image-dirty.fits", "image-image.fits")):
        floor = 100 * repeat[name]["regions"]["celestial_valid"]["relative_rms_difference"]
        ax.axhline(floor, color="0.5", ls="--", label="Repeated lossless input")
    for ax, title, ylabel in zip(axes, ("Visibility reconstruction", "WSClean dirty image", "WSClean multiscale restored image"),
                                ("Mean per-component MSE", "Difference RMS / reference RMS (%)", "Difference RMS / reference RMS (%)")):
        ax.set_title(title)
        ax.set_xlabel("Whole archive (MB)")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=.2)
        ax.set_ylim(bottom=0)
    axes[0].ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axes[2].legend(loc="upper left", bbox_to_anchor=(1, 1), fontsize=8)
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--benchmarks", type=Path, nargs="+", required=True)
    parser.add_argument("--reports", type=Path, required=True)
    parser.add_argument("--repeat", type=Path, required=True)
    parser.add_argument("--image-prefix", default="wsclean-")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plot", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.plot.exists():
        raise ValueError("Summary output paths must be new")
    result = summarize(args.matrix, args.benchmarks, args.reports, args.repeat, args.image_prefix)
    plot(result, args.plot)
    _write_new(args.output, json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": "complete", "points": len(result["points"]), "output": str(args.output)}))


if __name__ == "__main__":
    main()
