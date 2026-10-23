"""Export trusted local cProfile evidence and plot the measured speed comparison."""

import argparse
import json
import os
from pathlib import Path
import pstats
import statistics

os.environ.setdefault("MPLCONFIGDIR", "/tmp/varquanms-performance-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from varquanms.cli import _write_new


def profile(path):
    stats = pstats.Stats(str(path))
    calls = []
    for (filename, line, function), (primitive, total, own, cumulative, _) in stats.stats.items():
        calls.append({"file": Path(filename).name, "line": line, "function": function,
                      "primitive_calls": primitive, "total_calls": total,
                      "self_seconds": own, "cumulative_seconds": cumulative})
    return {"total_seconds": stats.total_tt, "total_calls": stats.total_calls,
            "top_self_time": sorted(calls, key=lambda v: v["self_seconds"], reverse=True)[:30]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-profile", type=Path)
    parser.add_argument("--optimized-profile", type=Path)
    parser.add_argument("--runtime", type=Path, default=Path("reports/performance-runtime.json"))
    parser.add_argument("--baseline-runtime", type=Path, default=Path("reports/paper-runtime.json"))
    parser.add_argument("--identity", type=Path)
    parser.add_argument("--verification", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--plot-only", action="store_true", help="regenerate charts using only versioned runtime reports")
    parser.add_argument("--force", action="store_true", help="replace generated figures; never replace reports")
    parser.add_argument("--figure", type=Path, default=Path("docs/paper/figs/performance"))
    args = parser.parse_args()
    runtime = json.loads(args.runtime.read_text())
    baseline = json.loads(args.baseline_runtime.read_text())
    if not args.plot_only:
        if any(getattr(args, key) is None for key in ("baseline_profile", "optimized_profile", "identity", "verification", "output")):
            parser.error("Report generation requires both profiles, identity, verification and output; use --plot-only for figures")
        if args.output.exists():
            raise ValueError("Report output must be new")
        verification = json.loads(args.verification.read_text())
        if verification["portable_verification"]["status"] != "passed":
            raise ValueError("Independent verification must pass")
    old = baseline["adaptive_summary"]["encode_process_seconds"]["median"]
    new = runtime["adaptive_summary"]["encode_process_seconds"]["median"]
    cpu = {name: statistics.mean(r["worker_timing_sums_seconds"][name] for r in runtime["records"])
           for name in ("native_search_quantization_cpu", "native_search_compression_cpu")}
    result = {"all_final_repeats_under_target": runtime["adaptive_summary"]["encode_process_seconds"]["maximum"] < 15}
    if not args.plot_only:
        result.update({"schema": "varquanms.performance-profile/1", "status": "complete",
              "baseline_profile": profile(args.baseline_profile),
              "optimized_profile": profile(args.optimized_profile),
              "profile_scope": "cProfile observes the calling thread. Lock acquire time in the optimized profile is primarily waiting for worker results; it is not a CPU-cost attribution to locks. Native kernel CPU counters are separately measured with CLOCK_THREAD_CPUTIME_ID.",
              "native_mean_thread_cpu_seconds": cpu,
              "baseline_median_process_seconds": old, "optimized_median_process_seconds": new,
              "median_latency_speedup": old/new, "target_process_seconds": 15,
              "all_final_repeats_under_target": runtime["adaptive_summary"]["encode_process_seconds"]["maximum"] < 15,
              "comparison": "Original one-worker Python encoder versus optimized 24-worker encoder; both algorithm changes and parallelism contribute. Same quantization/error/flag policies, same shared host, warm/unspecified cache.",
              "dataset_identity": json.loads(args.identity.read_text()),
              "portable_verification": {"status": "passed", "native_kernels_disabled": verification["native_kernels_disabled"],
                                        "export": verification["export"], "seconds": verification["seconds"],
                                        "source_and_restored": True},
              "timing_report": str(args.runtime)})
    stages = runtime["adaptive_summary"]["exclusive_encode_stages_seconds"]
    names = {"parallel_block_pipeline_wall": "Parallel block pipeline", "inspection_preflight": "Inspection / preflight",
             "amplitude_statistics": "Amplitude statistics", "quality_accounting": "Quality accounting",
             "input_reading": "Visibility reads", "other_exact_columns_and_uninstrumented_overhead": "Exact columns / other",
             "integrity_finalization_and_publication": "Checksums / final publication", "payload_writing": "Payload writes",
             "baseline_sorting_and_grouping": "Sort / group"}
    selected = sorted(names, key=lambda key: stages[key]["mean"], reverse=True)
    plt.rcParams.update({"svg.fonttype": "none", "svg.hashsalt": "varquanms-performance-v1"})
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5), constrained_layout=True)
    bars = axes[0].bar(["Original\n1 worker", "Optimized\n24 workers"], [old, new], color=["#71808b", "#247c77"])
    for bar, value in zip(bars, [old, new]):
        axes[0].text(bar.get_x()+bar.get_width()/2, value+7, f"{value:.2f} s", ha="center")
    axes[0].axhline(15, color="#b84d52", linestyle="--", label="15 s target")
    axes[0].set(ylabel="Median complete encode-process time (s)", title=f"{old/new:.2f}× lower latency; three repeats each", ylim=(0, old*1.15))
    axes[0].legend(loc="upper right", fontsize=9)
    axes[1].barh([names[k] for k in selected], [stages[k]["mean"] for k in selected], color="#23618d")
    axes[1].invert_yaxis(); axes[1].set(xlabel="Mean exclusive function-stage wall time (s)", title="Optimized latency breakdown")
    bars = axes[2].bar(["Quantize / pack /\nprecheck", "zlib\ncompression"], list(cpu.values()), color=["#805ba5", "#cf7434"])
    for bar, value in zip(bars, cpu.values()):
        axes[2].text(bar.get_x()+bar.get_width()/2, value+1, f"{value:.2f}", ha="center")
    axes[2].set(ylabel="Mean summed native thread CPU seconds", title="Kernel CPU cost across 24 workers")
    fig.suptitle("Same 409.595 MB MS and error bounds → 106.731 MB archive\nStartup, inspection, writing and integrity included; verification/imaging excluded; warm/unspecified cache")
    args.figure.parent.mkdir(parents=True, exist_ok=True)
    for extension in (".png", ".svg", ".pdf"):
        path = args.figure.with_suffix(extension)
        if path.exists() and not args.force:
            raise ValueError(f"Figure already exists: {path}")
        metadata = {"Date": None} if extension == ".svg" else {"CreationDate": None, "ModDate": None} if extension == ".pdf" else {}
        fig.savefig(path, dpi=180, bbox_inches="tight", metadata=metadata)
        if extension == ".svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines())+"\n")
    plt.close(fig)
    if not args.plot_only:
        _write_new(args.output, json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"report": str(args.output), "speedup": old/new, "all_runs_under_15s": result["all_final_repeats_under_target"]}))


if __name__ == "__main__":
    main()
