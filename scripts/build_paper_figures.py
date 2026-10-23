"""Reproducible paper diagrams, measured plots and real FITS comparison panels.

Vector diagrams are schematic; numerical charts read versioned JSON reports.
Image panels read the measured WSClean FITS files, never synthetic sky images.
"""

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/varquanms-paper-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle
import numpy as np

from varquanms.imaging import compare_plane

ROOT = Path(__file__).resolve().parents[1]
COLORS = {"blue": "#23618d", "teal": "#247c77", "orange": "#cf7434", "purple": "#805ba5", "red": "#b84d52", "grey": "#65727e"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "svg.fonttype": "none", "svg.hashsalt": "varquanms-paper-v1"})


def read(relative):
    return json.loads((ROOT / relative).read_text())


def save(fig, output, name, force, vector=True):
    paths = [output / (name + ext) for ext in ((".png", ".svg", ".pdf") if vector else (".png",))]
    if not force and any(p.exists() for p in paths):
        raise ValueError(f"Figure exists; use --force to regenerate {name}")
    for path in paths:
        metadata = {"Date": None} if path.suffix == ".svg" else {"CreationDate": None, "ModDate": None} if path.suffix == ".pdf" else {}
        fig.savefig(path, dpi=180, bbox_inches="tight", metadata=metadata)
        if path.suffix == ".svg":
            path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    plt.close(fig)
    return paths


def box(ax, x, y, w, h, text, color="blue", size=11):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08,rounding_size=0.10",
                              fc="white", ec=COLORS[color], lw=1.8))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size, linespacing=1.5)


def arrow(ax, a, b, label=None, color="grey", bend=0):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=15,
                                color=COLORS[color], lw=1.5,
                                connectionstyle=f"arc3,rad={bend}"))
    if label:
        ax.text((a[0] + b[0]) / 2, (a[1] + b[1]) / 2 + .13, label,
                ha="center", va="bottom", fontsize=9, bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})


def grouping():
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_axes([.01, .01, .98, .98]); ax.set(xlim=(0, 14), ylim=(0, 8)); ax.axis("off")
    ax.text(7, 7.75, "Grouping: arrange baseline spectra into blocks that share a scale", ha="center", fontsize=17, weight="bold")
    labels = ["Inspect supported layout\nrow context + shape", "One correlation\nV[row, channel, p]", "Separate classes\nauto / cross; valid / no-valid",
              "Order baseline IDs\namplitude or other key", "Split into B × F blocks\nindependent scale + bits"]
    for i, text in enumerate(labels):
        x = .2 + i * 2.8
        box(ax, x, 6.15, 2.35, 1.0, text, "teal" if i in (2, 4) else "blue", size=10)
        if i < 4:
            arrow(ax, (x + 2.45, 6.65), (x + 2.7, 6.65))
    ax.text(7, 5.65, "Illustration below: one valid cross-correlation class; toy amplitudes, B = 2 and F = 4", ha="center", fontsize=11)
    amplitudes = np.asarray([8, 1, 4, 2, 16, .5, 4, 2.])
    perm = np.argsort(amplitudes, kind="stable")
    toy = np.log2(amplitudes[:, None]) + np.linspace(-.2, .2, 8)[None, :]
    for rect, values, rows, title, grouped in [([.08, .235, .28, .40], toy, np.arange(8), "Original row order", False),
                                              ([.60, .235, .28, .40], toy[perm], perm, "Stable amplitude order", True)]:
        grid = fig.add_axes(rect)
        grid.imshow(values, cmap="viridis", aspect="auto", interpolation="nearest", vmin=-1.2, vmax=4.2)
        grid.set(xticks=np.arange(8), xticklabels=np.arange(8), yticks=np.arange(8),
                 yticklabels=[f"row {r}: A={amplitudes[r]:g}" for r in rows], xlabel="Channel index (order unchanged)")
        grid.set_title(title, fontsize=13, pad=8)
        grid.set_xticks(np.arange(-.5, 8, 1), minor=True); grid.set_yticks(np.arange(-.5, 8, 1), minor=True)
        grid.grid(which="minor", color="white", alpha=.35); grid.tick_params(which="minor", length=0)
        if grouped:
            for y in range(0, 8, 2):
                for x in (0, 4):
                    grid.add_patch(Rectangle((x-.5, y-.5), 4, 2, fill=False, ec="#ef8354", lw=2.6))
    arrow(ax, (5.2, 3.65), (7.05, 3.65), color="teal")
    ax.text(6.15, 4.25, r"$A_{rp}=\mathrm{median}_{f\in\mathcal{F}_{rp}}|V_{rfp}|$", ha="center", fontsize=11)
    ax.text(6.15, 2.8, "Move whole spectra\nnever individual samples", ha="center", fontsize=10)
    ax.text(7, 1.30, "perm = [5, 1, 3, 7, 2, 6, 0, 4]     (ties retain original row index)", ha="center", fontsize=11, family="monospace")
    box(ax, .7, .20, 12.6, .65, "Store perm + group boundaries + segment edges. Restore: V̂[perm[j], f, p] = decoded[j, f, p].\nUVW, original row keys, channel/product order, flags and weights remain exact by default.", "teal", size=10)
    return fig


def selection():
    fig, ax = plt.subplots(figsize=(13, 10))
    ax.set(xlim=(0, 13), ylim=(0, 10)); ax.axis("off")
    ax.text(6.5, 9.75, "Adaptive precision: smallest serialized candidate that passes all checks", ha="center", fontsize=16, weight="bold")
    box(ax, 3.5, 8.45, 6, .8, "One complex64 block + source flags\nand configured error ceilings", size=11)
    arrow(ax, (4.7, 8.4), (2.15, 7.7)); arrow(ax, (8.0, 8.4), (8.4, 7.7))
    box(ax, .5, 6.65, 3.35, 1.0, "Always construct exact candidate\noriginal float32 components\ncompress only if smaller", "teal", 10)
    box(ax, 5.0, 6.65, 6.8, 1.0, "Try b ∈ {4, 8, 12, 16} and scale quantiles q ∈ {1, 0.995}\nstep = stored float32(cutoff / (2ᵇ⁻¹ − 1))\nadjust upward if needed to cover the cutoff", "blue", 10)
    arrow(ax, (8.4, 6.55), (8.4, 6.03))
    box(ax, 5.0, 4.95, 6.8, 1.0, "Quantize real and imaginary with the same local step\nPreserve flagged / nonfinite / signed-zero / outlier values exactly\nPack codes; append positions + exact exceptions; compress", "blue", 10)
    arrow(ax, (8.4, 4.85), (8.4, 4.30))
    box(ax, 5.0, 3.45, 6.8, .78, "Decode using stored step and actual decoder arithmetic\nMeasure complex errors on source-unflagged finite samples", "purple", 10)
    arrow(ax, (8.4, 3.35), (8.4, 2.98))
    ax.add_patch(Polygon([[8.4, 3], [11.5, 2.3], [8.4, 1.6], [5.3, 2.3]], closed=True,
                         fc="white", ec=COLORS["purple"], lw=1.8))
    ax.text(8.4, 2.3, "Finite reconstruction?\nMSE + max error + bias + per-baseline MSE pass?", ha="center", va="center", fontsize=9)
    arrow(ax, (11.55, 2.3), (12.6, 2.3), "No", "red")
    ax.text(12.35, 1.65, "Reject\ncandidate", ha="center", fontsize=10, color=COLORS["red"])
    arrow(ax, (6.1, 1.95), (4.3, 1.10), "Yes", "teal")
    arrow(ax, (.65, 6.55), (.65, 1.60), color="teal")
    ax.text(2.8, 5.8, "Exact is always eligible\n(zero reconstruction error)", ha="center", fontsize=10, color=COLORS["teal"])
    box(ax, .55, .55, 5.4, 1.0, "Choose minimum actual serialized bytes\n16-byte header + compressed body\n(codes + exceptions, or exact components)\nFixed index cost is identical across candidates", "teal", 9)
    arrow(ax, (6.05, 1.03), (6.75, 1.03))
    box(ax, 6.85, .55, 5.1, .9, "Store winner; verify whole-product bounds\nthen finalize archive + checksums", "blue", 10)
    ax.text(2.8, 4.05, "No amplitude-to-bits lookup.\nNo silent clipping.\nNo inferred thermal-noise budget.", ha="center", fontsize=10, color=COLORS["grey"], linespacing=1.8)
    ax.text(6.5, .03, "Logical decision diagram: implementation skips error evaluation for candidates already larger than the best valid encoding.", ha="center", fontsize=9)
    return fig


def rate_distortion(old, new):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
    for i, point in enumerate(old["points"]):
        ys = [point["component_mse"], *[100 * point["imaging"][name]["regions"]["celestial_valid"]["relative_rms_difference"] for name in ("image-dirty.fits", "image-image.fits")]]
        for ax, y in zip(axes, ys):
            ax.scatter(point["archive_bytes"] / 1e6, y, color="C"+str(i%10), marker=["o","s","^","D"][i//4], s=42, label=point["name"])
    ys = [new["overall_component_mse"], *[100 * new["imaging"]["relative_rms_difference"][name] for name in ("image-dirty.fits", "image-image.fits")]]
    for ax, y, label in zip(axes, ys, ("Mean component MSE", "Dirty difference / reference RMS (%)", "Restored difference / reference RMS (%)")):
        ax.scatter(new["archive_bytes"]/1e6, y, s=200, marker="*", c="black", label="128×192, new 4/8/12/16")
        ax.set(xlabel="Whole archive (MB)", ylabel=label, ylim=(0,None)); ax.grid(alpha=.2)
    axes[0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    axes[-1].legend(loc="upper left", bbox_to_anchor=(1,1), fontsize=8)
    fig.suptitle("Historical grouping sweep (2/4/6/8) plus one new-ladder regression\nSame engineering ceilings ≠ equal measured error", fontsize=13)
    return fig


def precision_bytes(new):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.7), constrained_layout=True)
    bits = ("4","8","12","16","float32")
    percentages = [new["complex_sample_percent_by_block_precision"].get(b, 0) for b in bits]
    bars = axes[0].bar(["4","8","12","16","Exact"], percentages,
                       color=[COLORS[k] for k in ("blue","orange","purple","red","grey")])
    for bar, b, percent in zip(bars,bits,percentages):
        axes[0].text(bar.get_x()+bar.get_width()/2, bar.get_height()+2, f"{percent:.2f}%\n{new['component_bit_width_block_counts'].get(b,0):,} blocks", ha="center", fontsize=9)
    axes[0].set(ylim=(0,100), ylabel="Complex samples by block precision (%)", xlabel="Bits per component", title="Local precision allocation")
    categories = new["stored_bytes_by_category"]
    small = [(k,v) for k,v in categories.items() if k != "visibility_block_payload"]
    labels = [k.replace("visibility_", "").replace("quantization_", "").replace("container_and_manifest_", "").replace("_", " ") for k,v in small]
    axes[1].barh(labels, [v/1e6 for k,v in small], color=COLORS["teal"])
    axes[1].set(xlabel="Stored MB (excluding block payload)", title=f"Archive overhead detail\nBlock payload: {categories['visibility_block_payload']/1e6:.3f} MB")
    axes[1].invert_yaxis(); axes[1].grid(axis="x",alpha=.2)
    fig.suptitle(f"Measured archive: {new['archive_bytes']/1e6:.3f} MB; exceptions are included in block payload", fontsize=13)
    return fig


def images(reference, candidate, zoom=False):
    names = ["image-image.fits"] if zoom else ["image-dirty.fits", "image-image.fits"]
    fig, axes = plt.subplots(len(names), 3, figsize=(13, 4.1 * len(names)), squeeze=False, constrained_layout=True)
    for row, name in enumerate(names):
        measured, (a,b,valid) = compare_plane(reference/name, candidate/name)
        extent = None
        if zoom:
            y,x=measured["reference_peak_yx"]; radius=32
            lo_y,hi_y=max(0,y-radius),min(a.shape[0],y+radius+1)
            lo_x,hi_x=max(0,x-radius),min(a.shape[1],x+radius+1)
            a,b,valid=[v[lo_y:hi_y,lo_x:hi_x] for v in (a,b,valid)]
            extent=(lo_x-.5,hi_x-.5,lo_y-.5,hi_y-.5)
        maximum=max(np.max(np.abs(a[valid])),np.max(np.abs(b[valid])))
        norm=SymLogNorm(linthresh=maximum*.01,vmin=-maximum,vmax=maximum)
        for ax, values, title in zip(axes[row,:2],(a,b),("Original MS", "Decoded 4/8/12/16 archive")):
            image=ax.imshow(np.where(valid,values,np.nan),origin="lower",cmap="RdBu_r",norm=norm,extent=extent)
            ax.set_title(title);fig.colorbar(image,ax=ax,shrink=.72)
        delta=b-a; maximum=np.max(np.abs(delta[valid])) or 1.
        image=axes[row,2].imshow(np.where(valid,delta,np.nan),origin="lower",cmap="RdBu_r",vmin=-maximum,vmax=maximum,extent=extent)
        axes[row,2].set_title("Decoded − original\nIndependent linear color scale");fig.colorbar(image,ax=axes[row,2],shrink=.72)
        for ax in axes[row]:ax.set(xlabel="Pixel x",ylabel="Pixel y")
        kind="Restored" if "image-image" in name else "Dirty"
        axes[row,0].set_ylabel(f"{kind}\nPixel y")
    fig.suptitle(("Fixed-reference peak neighborhood (not an identified solar source)" if zoom else "WSClean: 2048², 4 arcmin, multiscale; finite celestial footprint")+"\nFITS Jy/beam label; absolute calibration unverified",fontsize=13)
    return fig


def timing(report):
    fig,axes=plt.subplots(1,2,figsize=(13,6),constrained_layout=True)
    records=report["records"]
    bars=axes[0].bar([r["name"].replace("adaptive-repeat-","Repeat ").replace("lossless-control","Lossless") for r in records],
                    [r["encode_process_seconds"] for r in records],color=[COLORS["grey"] if r["mode"]=="lossless" else COLORS["blue"] for r in records])
    for bar,r in zip(bars,records):axes[0].text(bar.get_x()+bar.get_width()/2,bar.get_height()+2,f"{r['encode_process_seconds']:.2f} s",ha="center")
    axes[0].set(ylabel="Complete encode-process wall time (s)",title="MS → complete archive; verification excluded")
    stages=report["adaptive_summary"]["exclusive_encode_stages_seconds"]
    keys=sorted(stages,key=lambda k:stages[k]["mean"],reverse=True)
    labels={"inspection_preflight":"Inspection/preflight","input_reading":"Visibility input reads","amplitude_statistics":"Amplitude statistics",
            "baseline_sorting_and_grouping":"Sort IDs + group","memory_gathering":"Gather blocks","quantization_and_scale_estimation":"Scale + quantization",
            "packing_and_secondary_compression":"Packing + lossless compression","candidate_reconstruction_and_error_checks":"Candidate decode + error checks",
            "payload_writing":"Visibility payload writes","quality_accounting":"Global quality accounting",
            "parallel_block_pipeline_wall":"Parallel block pipeline",
            "native_candidate_search":"Native candidate search",
            "integrity_finalization_and_publication":"Integrity/final publication","other_exact_columns_and_uninstrumented_overhead":"Exact columns + other overhead"}
    axes[1].barh([labels[k] for k in keys],[stages[k]["mean"] for k in keys],color=COLORS["teal"])
    axes[1].invert_yaxis();axes[1].set(xlabel="Mean exclusive encode-function subregion (s)",title="Nested codec_total is not added again")
    fig.suptitle("Sequential fresh processes, shared host, warm/unspecified OS cache\nNo cold-cache or isolated-throughput claim",fontsize=13)
    return fig


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"docs/paper/figs")
    parser.add_argument("--timing",type=Path)
    parser.add_argument("--reference",type=Path,default=ROOT/"results/wsclean-grouping/reference-single")
    parser.add_argument("--candidate",type=Path,default=ROOT/"results/wsclean-precision-4-8-12-16")
    parser.add_argument("--force",action="store_true")
    parser.add_argument("--skip-images",action="store_true")
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    old,new=read("reports/grouping-results.json"),read("reports/precision-ladder-4-8-12-16.json")
    figures={"grouping":grouping(),"precision-selection":selection(),
             "rate-distortion":rate_distortion(old,new),"precision-and-bytes":precision_bytes(new)}
    if args.timing:figures["runtime"]=timing(json.loads(args.timing.read_text()))
    paths=[]
    for name,fig in figures.items():paths.extend(save(fig,args.output,name,args.force))
    if not args.skip_images:
        for name,zoom in (("wsclean-full",False),("wsclean-peak",True)):
            paths.extend(save(images(args.reference,args.candidate,zoom),args.output,name,args.force,vector=False))
    print(json.dumps({"figures":[str(p) for p in paths]},indent=2))


if __name__=="__main__":main()
