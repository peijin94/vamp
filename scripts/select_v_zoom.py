"""Select a reference-only V zoom well away from the horizon for the paper."""
import json
from pathlib import Path

import astropy
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation
from astropy.io import fits
from astropy.time import Time
from astropy.utils import iers
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
import numpy as np

from varquanms.imaging import file_sha256, read_fits


OUTPUT = Path('reports/calibrated-v-zoom.json')


def select_zoom(study):
    reference_path = Path(study['reference_directory']) / 'image-image.fits'
    reference, header = read_fits(reference_path)
    reference_sha = file_sha256(reference_path)
    if reference_sha != study['reference_file_sha256']:
        raise ValueError('V reference changed since the measured study')
    mask_path = Path(study['comparison_mask']['path'])
    provenance_path = mask_path.with_suffix('.json')
    provenance = json.loads(provenance_path.read_text())
    if file_sha256(mask_path) != provenance['mask_sha256'] or provenance['mask_sha256'] != study['comparison_mask']['sha256']:
        raise ValueError('Elevation-mask provenance does not match the study')
    orientation_path = Path(provenance['iers_file'])
    if file_sha256(orientation_path) != provenance['iers_sha256']:
        raise ValueError('Earth-orientation table differs from the original mask')
    mask = fits.getdata(mask_path).astype(bool)
    if mask.shape != reference.shape or not np.all(np.isfinite(reference[mask])):
        raise ValueError('Invalid V reference comparison domain')
    epochs = Time(np.array(provenance['times_utc_mjd_seconds']) / 86400., format='mjd', scale='utc')
    location = EarthLocation.from_geocentric(*provenance['antenna_mean_itrf_metres'], unit=u.m)
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    orientation = iers.IERS_A.open(str(orientation_path))
    _, status = orientation.ut1_utc(epochs, return_status=True)
    if np.any(status < 0):
        raise ValueError('Earth-orientation table does not cover the observation')
    wcs = WCS(header).celestial
    minimum = np.full(reference.shape, -np.inf)
    # Use the original mask's epochs, location and no-refraction convention.
    # This is an image-display calculation; no MS or UVW values are changed.
    with iers.earth_orientation_table.set(orientation):
        for start in range(0, reference.shape[0], 64):
            rows = slice(start, start + 64)
            y, x = np.nonzero(mask[rows])
            if not len(y):
                continue
            y += start
            sky = pixel_to_skycoord(x, y, wcs)
            values = np.full(len(y), np.inf)
            for epoch in epochs:
                altitude = sky.transform_to(AltAz(obstime=epoch, location=location,
                                                 pressure=provenance['pressure_hpa'] * u.hPa)).alt.deg
                values = np.minimum(values, altitude)
            minimum[y, x] = values
    eligible = mask & (minimum > 60.)
    if not eligible.any():
        raise ValueError('No V reference pixels are above 60 degrees throughout the observation')
    peak = np.unravel_index(np.argmax(np.where(eligible, np.abs(reference), -np.inf)), reference.shape)
    y, x = map(int, peak)
    bounds = dict(y_start=y-64, y_stop=y+64, x_start=x-64, x_stop=x+64)
    if y < 64 or x < 64 or y+64 > reference.shape[0] or x+64 > reference.shape[1]:
        raise ValueError('High-elevation peak cannot support a complete 128-pixel crop')
    crop = np.s_[y-64:y+64, x-64:x+64]
    lowest = float(minimum[crop].min())
    if lowest <= 45.:
        raise ValueError('The entire V zoom must stay above 45 degrees at every epoch')
    return dict(schema='varquanms.paper-v-zoom/1', reference_fits=str(reference_path),
                reference_fits_sha256=reference_sha, comparison_mask_sha256=provenance['mask_sha256'],
                mask_provenance_path=str(provenance_path), mask_provenance_sha256=file_sha256(provenance_path),
                selection='Largest absolute reference V pixel above 60 degrees at all integration TIME values; first row-major tie',
                center_yx=[y, x], bounds=bounds, pixel_indexing='zero-based; stop indices exclusive',
                center_reference_jy_per_beam=float(reference[peak]),
                center_minimum_elevation_deg=float(minimum[peak]),
                entire_crop_minimum_elevation_deg=lowest,
                required_entire_crop_minimum_elevation_deg=45.,
                eligible_pixels=int(eligible.sum()), astropy_version=astropy.__version__,
                times_utc_mjd_seconds=provenance['times_utc_mjd_seconds'],
                antenna_mean_itrf_metres=provenance['antenna_mean_itrf_metres'],
                pressure_hpa=provenance['pressure_hpa'], iers_file=str(orientation_path),
                iers_sha256=provenance['iers_sha256'],
                note='Display selection only; all reported errors retain the full above-10-degree mask and original residual RMS.')


if __name__ == '__main__':
    record = select_zoom(json.loads(Path('reports/calibrated-epsd-v.json').read_text()))
    OUTPUT.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({k: record[k] for k in ('center_yx', 'bounds', 'center_reference_jy_per_beam',
                                           'center_minimum_elevation_deg', 'entire_crop_minimum_elevation_deg')}, indent=2))
