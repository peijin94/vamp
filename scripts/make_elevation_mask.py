"""Make a fixed geometric elevation mask from original MS times and FITS WCS."""
import argparse
import json
from pathlib import Path
import time

import astropy
from astropy import units as u
from astropy.coordinates import AltAz, EarthLocation
from astropy.io import fits
from astropy.time import Time
from astropy.utils import iers
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord
from casacore.tables import table
import numpy as np

from vamp.imaging import file_sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ms", type=Path)
    parser.add_argument("image", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--minimum-elevation", type=float, default=10.)
    parser.add_argument("--iers-file", type=Path, help="Local IERS finals2000A table; retained with the experiment")
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix('.json').exists():
        raise ValueError("Mask outputs must be new")
    if not 0 <= args.minimum_elevation < 90:
        raise ValueError("Elevation must be in [0,90)")
    started = time.perf_counter()
    with table(str(args.ms), ack=False, readonly=True, lockoptions='usernoread') as t:
        keywords = t.getcolkeywords('TIME')
        if keywords.get('MEASINFO', {}).get('Ref') != 'UTC' or keywords.get('QuantumUnits') != ['s']:
            raise ValueError("Mask requires UTC TIME in seconds")
        times = np.unique(t.getcol('TIME'))
    with table(str(args.ms / 'ANTENNA'), ack=False, readonly=True, lockoptions='usernoread') as t:
        keywords = t.getcolkeywords('POSITION')
        if keywords.get('MEASINFO', {}).get('Ref') != 'ITRF' or keywords.get('QuantumUnits') != ['m'] * 3:
            raise ValueError("Mask requires antenna ITRF positions in metres")
        position = t.getcol('POSITION').mean(axis=0)
    location = EarthLocation.from_geocentric(*position, unit=u.m)
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None
    epochs = Time(times / 86400., format='mjd', scale='utc')
    orientation = iers.IERS_A.open(str(args.iers_file)) if args.iers_file else iers.IERS_Auto.open()
    iers.earth_orientation_table.set(orientation)
    _, status = orientation.ut1_utc(epochs, return_status=True)
    if np.any(status < 0):
        raise ValueError("Bundled Earth-orientation table does not cover MS times")
    header = fits.getheader(args.image)
    shape = (header['NAXIS2'], header['NAXIS1'])
    wcs = WCS(header).celestial
    mask = np.zeros(shape, dtype=bool)
    for start in range(0, shape[0], 64):
        y, x = np.mgrid[start:min(start + 64, shape[0]), :shape[1]]
        coords = pixel_to_skycoord(x, y, wcs)
        valid = np.isfinite(coords.spherical.lon.deg) & np.isfinite(coords.spherical.lat.deg)
        selected = coords[valid]
        above = np.ones(len(selected), dtype=bool)
        for epoch in epochs:
            altitude = selected.transform_to(AltAz(obstime=epoch, location=location, pressure=0*u.hPa)).alt.deg
            above &= altitude > args.minimum_elevation
        valid[valid] = above
        mask[start:start + len(valid)] = valid
    mask_header = wcs.to_header()
    mask_header['ELCUT'] = args.minimum_elevation
    mask_header['ELPOLICY'] = 'all integrations; strict > cut'
    fits.writeto(args.output, mask.astype(np.uint8), mask_header)
    record = dict(schema='vamp.elevation-mask/1', source_ms=str(args.ms), reference_fits=str(args.image),
                  reference_fits_sha256=file_sha256(args.image), mask_sha256=file_sha256(args.output),
                  minimum_elevation_deg=args.minimum_elevation, time_policy='intersection over all integration TIME values',
                  times_utc_mjd_seconds=times.tolist(), antenna_mean_itrf_metres=position.tolist(),
                  pressure_hpa=0., astropy_version=astropy.__version__, iers_class=type(orientation).__name__,
                  iers_file=str(args.iers_file) if args.iers_file else None,
                  iers_sha256=file_sha256(args.iers_file) if args.iers_file else None,
                  iers_mjd_range=[float(orientation['MJD'][0].value), float(orientation['MJD'][-1].value)],
                  iers_ut1_status=np.asarray(status).tolist(), shape=list(shape), pixels=int(mask.sum()),
                  elapsed_seconds=time.perf_counter()-started,
                  note='Coordinate calculation is used only for the image comparison mask; original UVW remains untouched.')
    args.output.with_suffix('.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
