"""Measured WSClean FITS differences with matched coordinates and image units.

Astropy is an optional validation dependency; it is not needed by the decoder.
"""

import hashlib
import json
import os
from pathlib import Path
import re

import numpy as np


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def _grid(header):
    prefixes = r"(?:NAXIS\d*|CTYPE\d+|CUNIT\d+|CRPIX\d+|CRVAL\d+|CDELT\d+|CROTA\d+|(?:PC|CD|PV|PS)\d+_\d+)"
    exact = {"BUNIT", "BMAJ", "BMIN", "BPA", "EQUINOX", "LONPOLE", "LATPOLE", "RADESYS", "SPECSYS", "RESTFRQ"}
    return {key: header[key] for key in header if key in exact or re.fullmatch(prefixes, key)}


def read_fits(path):
    from astropy.io import fits
    with fits.open(path, memmap=False) as hdus:
        header = hdus[0].header.copy()
        pixels = np.asarray(hdus[0].data, dtype=np.float64).squeeze()
    if pixels.ndim != 2:
        raise ValueError("Expected a single two-dimensional FITS image plane")
    if not header.get("BUNIT"):
        raise ValueError("Image units must be explicit")
    return pixels, header


def sky_mask(shape, header):
    from astropy.wcs import WCS
    y, x = np.indices(shape, dtype=np.float64)
    longitude, latitude = WCS(header).celestial.pixel_to_world_values(x, y)
    return np.isfinite(longitude) & np.isfinite(latitude)


def difference_metrics(reference, candidate, mask):
    count = int(mask.sum())
    if not count:
        return {"pixels": 0}
    a, b = reference[mask], candidate[mask]
    error = b - a
    rms = float(np.sqrt(np.mean(error**2)))
    ref_rms = float(np.sqrt(np.mean(a**2)))
    absolute = np.abs(error)
    return {"pixels": count, "difference_rms": rms, "difference_bias": float(error.mean()),
            "difference_max_abs": float(absolute.max()), "difference_p99_abs": float(np.quantile(absolute, .99)),
            "difference_p999_abs": float(np.quantile(absolute, .999)),
            "reference_rms": ref_rms, "relative_rms_difference": rms / ref_rms if ref_rms else None,
            "reference_abs_peak": float(np.max(np.abs(a))),
            "rms_normalized_to_reference_peak": rms / float(np.max(np.abs(a))) if np.any(a) else None}


def compare_plane(reference_path, candidate_path, aperture_radius=4, comparison_mask=None):
    reference, header = read_fits(reference_path)
    candidate, candidate_header = read_fits(candidate_path)
    if reference.shape != candidate.shape or _grid(header) != _grid(candidate_header):
        raise ValueError("FITS image grid, coordinates, beam or units differ")
    if comparison_mask is not None and (comparison_mask.shape != reference.shape or comparison_mask.dtype != np.bool_):
        raise ValueError("Comparison mask must be boolean and match the image shape")
    finite_domain = comparison_mask if comparison_mask is not None else np.ones(reference.shape, dtype=bool)
    if not np.array_equal(np.isfinite(reference[finite_domain]), np.isfinite(candidate[finite_domain])):
        raise ValueError("FITS finite pixel masks differ")
    valid = sky_mask(reference.shape, header) & np.isfinite(reference)
    region_name = "celestial_valid"
    if comparison_mask is not None:
        valid &= comparison_mask
        region_name = "above_elevation_cut"
    if not np.any(valid):
        raise ValueError("No finite pixels with valid celestial coordinates")
    peak = np.unravel_index(np.argmax(np.where(valid, np.abs(reference), -np.inf)), reference.shape)
    candidate_peak = np.unravel_index(np.argmax(np.where(valid, np.abs(candidate), -np.inf)), reference.shape)
    y, x = np.ogrid[:reference.shape[0], :reference.shape[1]]
    aperture = valid & ((y - peak[0])**2 + (x - peak[1])**2 <= aperture_radius**2)
    bright = valid & (np.abs(reference) >= .1 * abs(reference[peak]))
    faint = valid & (np.abs(reference) < .01 * abs(reference[peak]))
    intermediate = valid & ~bright & ~faint
    reference_sum, candidate_sum = float(reference[aperture].sum()), float(candidate[aperture].sum())
    peak_value = float(reference[peak])
    result = {"reference": str(reference_path), "candidate": str(candidate_path), "shape": list(reference.shape),
              "grid_and_beam": _grid(header), "finite_mask_matches": True,
              "finite_mask_comparison_domain": "provided comparison mask" if comparison_mask is not None else "whole image",
              "excluded_pixels": int(reference.size - valid.sum()),
              "primary_region": region_name,
              "regions": {region_name: difference_metrics(reference, candidate, valid),
                          "bright_reference": difference_metrics(reference, candidate, bright),
                          "intermediate_reference": difference_metrics(reference, candidate, intermediate),
                          "faint_reference": difference_metrics(reference, candidate, faint),
                          "reference_peak_aperture": difference_metrics(reference, candidate, aperture)},
              "brightness_regions": {"bright_reference": "abs(reference) >= 10% of masked reference absolute peak",
                                     "intermediate_reference": "1% <= abs(reference)/masked reference absolute peak < 10%",
                                     "faint_reference": "abs(reference) < 1% of masked reference absolute peak; not a noise estimate"},
              "reference_peak_yx": [int(v) for v in peak], "candidate_peak_yx": [int(v) for v in candidate_peak],
              "peak_shift_pixels": float(np.hypot(candidate_peak[0] - peak[0], candidate_peak[1] - peak[1])),
              "reference_peak_value": peak_value, "candidate_at_reference_peak": float(candidate[peak]),
              "relative_change_at_reference_peak": float((candidate[peak] - peak_value) / abs(peak_value)) if peak_value else None,
              "aperture": {"radius_pixels": aperture_radius, "reference_sum_pixels": reference_sum,
                           "candidate_sum_pixels": candidate_sum,
                           "relative_sum_change": (candidate_sum - reference_sum) / abs(reference_sum) if reference_sum else None,
                           "interpretation": "Fixed pixel aperture around brightest absolute reference pixel; not identified as Sun and not a calibrated integrated-flux measurement"},
              "clean": {"reference_iterations": header.get("WSCMINOR"), "candidate_iterations": candidate_header.get("WSCMINOR"),
                        "reference_iteration_limit": header.get("WSCNITER"), "candidate_iteration_limit": candidate_header.get("WSCNITER")}}
    for label, path in (("reference", reference_path), ("candidate", candidate_path)):
        result[label + "_file_sha256"] = file_sha256(path)
    return result, (reference, candidate, valid)


def compare_run(reference_dir, candidate_dir, mask_path=None):
    reference_dir, candidate_dir = Path(reference_dir), Path(candidate_dir)
    runs = [json.loads((path / "complete.json").read_text()) for path in (reference_dir, candidate_dir)]
    if any(r.get("status") != "complete" or r.get("source_unchanged") is not True for r in runs):
        raise ValueError("Imaging runs must be complete with unchanged source checksums")
    if runs[0]["parameters"] != runs[1]["parameters"] or runs[0]["container_id"] != runs[1]["container_id"]:
        raise ValueError("Imaging settings or container versions differ")
    files = sorted(p.name for p in reference_dir.glob("image-*.fits"))
    if files != sorted(p.name for p in candidate_dir.glob("image-*.fits")):
        raise ValueError("Imaging output file sets differ")
    if not files or not any(name.endswith("-image.fits") for name in files):
        raise ValueError("No restored FITS images found")
    if 'image-horizon-mask.fits' in files:
        from astropy.io import fits
        if not np.array_equal(fits.getdata(reference_dir/'image-horizon-mask.fits'),
                              fits.getdata(candidate_dir/'image-horizon-mask.fits')):
            raise ValueError('WSClean horizon masks differ')
        files.remove('image-horizon-mask.fits')
    products = {}
    mask, mask_info = None, None
    if mask_path is not None:
        from astropy.io import fits
        with fits.open(mask_path) as hdus:
            raw_mask = np.asarray(hdus[0].data).squeeze()
            if not np.all((raw_mask == 0) | (raw_mask == 1)):
                raise ValueError('Comparison mask must contain only zero and one')
            mask = raw_mask.astype(bool)
            mask_header = hdus[0].header.copy()
        # Compare celestial grids; the binary mask has no restoring beam or flux unit.
        from astropy.wcs import WCS
        image_header = read_fits(reference_dir / files[0])[1]
        if WCS(mask_header).celestial.to_header() != WCS(image_header).celestial.to_header():
            raise ValueError("Comparison mask celestial grid differs")
        mask_info = {"path": str(mask_path), "sha256": file_sha256(mask_path),
                     "minimum_elevation_deg": mask_header.get("ELCUT"),
                     "time_policy": mask_header.get("ELPOLICY"), "pixels": int(mask.sum())}
    for name in files:
        measured, _ = compare_plane(reference_dir / name, candidate_dir / name, comparison_mask=mask)
        products[name] = measured
    result = {"schema": "varquanms.wsclean-comparison/1", "status": "measured",
              "reference": str(reference_dir), "candidate": str(candidate_dir),
              "parameters": runs[0]["parameters"], "container_id": runs[0]["container_id"], "products": products,
              "comparison_mask": mask_info,
              "limitations": ["No scientific acceptance thresholds have been approved; measured does not mean scientifically accepted",
                              "Dirty, restored, residual and model images include sky structure; image RMS is not thermal noise",
                              "Automatic masks and nonlinear CLEAN paths may differ; compare the repeated lossless reference",
                              ("All metrics use the supplied fixed elevation mask; no primary-beam correction is implied" if mask is not None else
                               "Celestial WCS validity excludes SIN projection corners; this is not a primary-beam or sky mask"),
                              "No calibration-solution, solar-identification, absolute-flux or temporal validation is implied"]}
    return result


def plot_comparison(reference_path, candidate_path, output, comparison_mask=None):
    # Reuse a writable task-specific font cache on restricted execution hosts.
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/varquanms-matplotlib-cache")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    measured, (reference, candidate, valid) = compare_plane(reference_path, candidate_path, comparison_mask=comparison_mask)
    delta = candidate - reference
    extent = max(float(np.max(np.abs(reference[valid]))), float(np.max(np.abs(candidate[valid])))) or 1.0
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.4), constrained_layout=True)
    for ax, data, title in zip(axes[:2], (reference, candidate), ("Lossless reference", "Archive reconstruction")):
        im = ax.imshow(np.where(valid, data, np.nan), origin="lower", cmap="RdBu_r",
                       norm=SymLogNorm(linthresh=extent * .01, vmin=-extent, vmax=extent))
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=.75)
    limit = float(np.max(np.abs(delta[valid]))) or 1.0
    im = axes[2].imshow(np.where(valid, delta, np.nan), origin="lower", cmap="RdBu_r", vmin=-limit, vmax=limit)
    axes[2].set_title("Difference (independent linear scale)")
    fig.colorbar(im, ax=axes[2], shrink=.75)
    for ax in axes:
        ax.set_xlabel("Pixel x")
        ax.set_ylabel("Pixel y")
    fig.suptitle(f"WSClean multiscale CLEAN • {reference.shape[1]} × {reference.shape[0]} • {measured['grid_and_beam']['BUNIT']}")
    fig.savefig(output, dpi=160)
    plt.close(fig)
