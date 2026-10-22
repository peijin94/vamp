import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from varquanms.imaging import compare_plane, compare_run, difference_metrics


@unittest.skipUnless(importlib.util.find_spec("astropy"), "Install the optional validation extra")
class ImagingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def write(self, name, data, **changes):
        from astropy.io import fits
        header = fits.Header({"CTYPE1": "RA---SIN", "CTYPE2": "DEC--SIN", "CRVAL1": 150., "CRVAL2": 37.,
                              "CRPIX1": 5., "CRPIX2": 5., "CDELT1": -.1, "CDELT2": .1,
                              "CUNIT1": "deg", "CUNIT2": "deg", "BUNIT": "JY/BEAM",
                              "BMAJ": .2, "BMIN": .1, "BPA": 0.})
        header.update(changes)
        path = self.root / name
        fits.writeto(path, np.asarray(data, dtype=np.float32), header)
        return path

    def test_known_difference_and_fixed_peak_aperture(self):
        a = np.ones((10, 10), dtype=np.float32)
        a[4, 4] = 4
        b = a + .25
        report, _ = compare_plane(self.write("a.fits", a), self.write("b.fits", b))
        metrics = report["regions"]["celestial_valid"]
        self.assertEqual(metrics["pixels"], 100)
        self.assertEqual(metrics["difference_rms"], .25)
        self.assertEqual(metrics["difference_bias"], .25)
        self.assertEqual(report["relative_change_at_reference_peak"], .0625)
        self.assertEqual(report["peak_shift_pixels"], 0)

    def test_grid_beam_units_and_finite_mask_mismatch_rejected(self):
        a = np.ones((10, 10))
        reference = self.write("reference.fits", a)
        for i, changes in enumerate(({"CDELT1": -.2}, {"BMAJ": .3}, {"BUNIT": "JY/PIXEL"})):
            with self.assertRaisesRegex(ValueError, "grid, coordinates, beam or units"):
                compare_plane(reference, self.write(f"bad{i}.fits", a, **changes))
        a[0, 0] = np.nan
        with self.assertRaisesRegex(ValueError, "finite pixel masks"):
            compare_plane(reference, self.write("nonfinite.fits", a))

    def test_sin_projection_corners_are_excluded(self):
        a = np.ones((10, 10))
        reference = self.write("wide.fits", a, CDELT1=-20., CDELT2=20.)
        report, _ = compare_plane(reference, reference)
        self.assertGreater(report["excluded_pixels"], 0)
        self.assertLess(report["regions"]["celestial_valid"]["pixels"], 100)
        self.assertEqual(report["regions"]["celestial_valid"]["difference_rms"], 0)

    def test_run_settings_cannot_be_silently_mixed(self):
        a, b = self.root / "a", self.root / "b"
        a.mkdir()
        b.mkdir()
        record = {"status": "complete", "source_unchanged": True, "parameters": {"size": 2048}, "container_id": "same"}
        (a / "complete.json").write_text(json.dumps(record))
        record["parameters"] = {"size": 1024}
        (b / "complete.json").write_text(json.dumps(record))
        with self.assertRaisesRegex(ValueError, "settings"):
            compare_run(a, b)

    def test_empty_and_zero_reference_regions_have_no_invented_ratios(self):
        a = np.zeros((2, 2))
        self.assertEqual(difference_metrics(a, a, np.zeros_like(a, dtype=bool)), {"pixels": 0})
        result = difference_metrics(a, np.ones_like(a), np.ones_like(a, dtype=bool))
        self.assertIsNone(result["relative_rms_difference"])
        self.assertEqual(result["difference_rms"], 1)


if __name__ == "__main__":
    unittest.main()
