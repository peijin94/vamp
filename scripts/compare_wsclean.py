"""Compare all matched WSClean products and optionally plot restored images."""

import argparse
import json
from pathlib import Path

from varquanms.cli import _write_new
from varquanms.imaging import compare_run, plot_comparison


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()
    if args.report.exists() or args.plot and args.plot.exists():
        raise ValueError("Comparison outputs must be new")
    result = compare_run(args.reference, args.candidate)
    if args.plot:
        names = [n for n in result["products"] if n.endswith("-image.fits")]
        name = "image-I-image.fits" if "image-I-image.fits" in names else names[0]
        plot_comparison(args.reference / name, args.candidate / name, args.plot)
    _write_new(args.report, json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": result["status"], "report": str(args.report),
                      "products": {name: entry["regions"]["celestial_valid"]["relative_rms_difference"] for name, entry in result["products"].items()}}, indent=2))


if __name__ == "__main__":
    main()
