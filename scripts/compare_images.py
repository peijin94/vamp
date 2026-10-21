"""Compare two CASA dirty images in their common native pixel ordering."""

import argparse
import json
from pathlib import Path

from casacore.tables import table
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--plot", type=Path, required=True)
    args = parser.parse_args()
    if args.report.exists() or args.plot.exists():
        raise ValueError("Image comparison outputs must be new")
    with table(str(args.reference), readonly=True, lockoptions="usernoread", ack=False) as t:
        reference = np.asarray(t.getcell("map", 0), dtype=np.float64).squeeze()
    with table(str(args.candidate), readonly=True, lockoptions="usernoread", ack=False) as t:
        candidate = np.asarray(t.getcell("map", 0), dtype=np.float64).squeeze()
    if reference.shape != candidate.shape or reference.ndim != 2:
        raise ValueError("Image shapes differ or are not 2D")
    error = candidate - reference
    rms = float(np.sqrt(np.mean(error**2)))
    reference_rms = float(np.sqrt(np.mean(reference**2)))
    peak = float(reference.max())
    result = {"schema": "varquanms.image-comparison/1", "reference": str(args.reference),
              "candidate": str(args.candidate), "shape": list(reference.shape),
              "difference_rms": rms, "difference_max_abs": float(np.max(np.abs(error))),
              "difference_bias": float(error.mean()), "reference_rms": reference_rms,
              "relative_rms_difference": rms / reference_rms if reference_rms else None,
              "reference_peak": peak, "candidate_peak": float(candidate.max()),
              "relative_peak_change": float((candidate.max() - peak) / abs(peak)) if peak else None,
              "limitations": ["64x64 natural-weight MFS Stokes I dirty image only; no calibration, deconvolution or weak-polarization acceptance",
                              "CASA 5.4.1 warned of outdated IERS/leap-second tables and missing OVRO-LWA observatory metadata",
                              "Image RMS includes sky structure and is not a validated thermal-noise estimate",
                              "Array ordering is native to both CASA image tables; this plot has pixel coordinates, not validated celestial axes"]}
    args.report.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), constrained_layout=True)
    extent = float(max(np.max(np.abs(reference)), np.max(np.abs(candidate))))
    for axis, data, title in zip(axes[:2], (reference, candidate), ("Lossless reference", "Adaptive reconstruction")):
        image = axis.imshow(data, origin="lower", cmap="RdBu_r", vmin=-extent, vmax=extent)
        axis.set_title(title)
        fig.colorbar(image, ax=axis, shrink=.8)
    maximum = float(np.max(np.abs(error))) or 1.0
    image = axes[2].imshow(error, origin="lower", cmap="RdBu_r", vmin=-maximum, vmax=maximum)
    axes[2].set_title("Difference (independent scale)")
    fig.colorbar(image, ax=axes[2], shrink=.8)
    for axis in axes:
        axis.set_xlabel("Pixel")
        axis.set_ylabel("Pixel")
    fig.savefig(args.plot, dpi=160)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
