#!/usr/bin/env python
"""CASA 5/6 table reopening and optional small dirty-image workflow (Python 2/3)."""
from __future__ import print_function

import argparse
import hashlib
import json
import os
import shutil
import sys
import time

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ms")
    parser.add_argument("--report", required=True)
    parser.add_argument("--image", help="new image path for a 64x64 natural-weight Stokes I dirty image")
    parser.add_argument("--data-column", default="DATA", choices=("DATA", "CORRECTED_DATA"))
    args = parser.parse_args()
    if os.path.exists(args.report) or args.image and os.path.exists(args.image):
        raise ValueError("CASA check outputs must be new")
    result = {"schema": "vamp.casa-check/1", "ms": os.path.abspath(args.ms),
              "data_column": args.data_column,
              "python_executable": sys.executable, "table_reopen": "not run",
              "imaging": "not run", "science_acceptance": "not established by this smoke workflow"}
    started = time.time()
    try:
        working = os.path.abspath((args.image or args.report) + ".work.ms")
        if os.path.exists(working):
            raise ValueError("CASA working copy must be new")
        # CASA imager appends HISTORY even with usescratch=False. Work on a copy.
        shutil.copytree(os.path.abspath(args.ms), working)
        result["working_ms"] = working
        try:
            from casatools import table, imager, image
            Table, Imager, Image = table, imager, image
        except ImportError:
            from casac import casac
            Table, Imager, Image = casac.table, casac.imager, casac.image
        tb = Table()
        try:
            tb.open(working, nomodify=True)
            result.update(table_reopen="passed", row_count=int(tb.nrows()), column_names=tb.colnames(),
                          casa_data_cell_shape=list(tb.getcell(args.data_column, 0).shape))
        finally:
            tb.close()
        if args.image:
            if not any(c in result['column_names'] for c in ('DATA', 'FLOAT_DATA')):
                raise ValueError('CASA imager requires DATA or FLOAT_DATA; use the calibrated DATA compatibility export')
            stamp = time.time()
            im = Imager()
            try:
                im.open(working, usescratch=False)
                im.selectvis(field="0", spw="0")
                im.defineimage(nx=64, ny=64, cellx="600arcsec", celly="600arcsec", stokes="I", mode="mfs")
                im.weight(type="natural")
                im.makeimage(type="corrected" if args.data_column == "CORRECTED_DATA" else "observed", image=os.path.abspath(args.image))
            finally:
                im.close()
            ia = Image()
            try:
                ia.open(os.path.abspath(args.image))
                pixels = np.asarray(ia.getchunk())
                if not np.all(np.isfinite(pixels)):
                    raise ValueError("Nonfinite image pixels")
                result.update(imaging="passed", image=os.path.abspath(args.image), image_shape=list(pixels.shape),
                              image_sha256=hashlib.sha256(pixels.astype("<f4").tobytes()).hexdigest(),
                              image_min=float(pixels.min()), image_max=float(pixels.max()),
                              image_rms=float(np.sqrt(np.mean(pixels.astype(np.float64)**2))),
                              imaging_seconds=time.time() - stamp,
                              imaging_parameters={"nx": 64, "ny": 64, "cell": "600arcsec", "stokes": "I", "mode": "mfs", "weighting": "natural", "type": args.data_column + "/dirty"})
            finally:
                ia.close()
    except Exception as error:
        result["error"] = type(error).__name__ + ": " + str(error)
        if result["table_reopen"] == "not run":
            result["table_reopen"] = "failed"
        elif args.image:
            result["imaging"] = "failed"
    result["elapsed_seconds"] = time.time() - started
    encoded = (json.dumps(result, indent=2, separators=(",", ": "), allow_nan=False) + "\n").encode("utf-8")
    fd = os.open(args.report, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(fd, "wb") as stream:
        stream.write(encoded)
    print(encoded.decode("utf-8"))
    return 1 if "error" in result else 0


if __name__ == "__main__":
    sys.exit(main())
