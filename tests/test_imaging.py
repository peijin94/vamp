import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from vamp.imaging import compare_plane, compare_run, difference_metrics, image_parameters_match


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

    def test_elevation_mask_excludes_errors_and_peak_outside_region(self):
        a = np.ones((10, 10), dtype=np.float32)
        a[4, 4] = 2
        b = a.copy()
        b[0, 0] = 1000
        mask = np.ones_like(a, dtype=bool)
        mask[0, 0] = False
        report, _ = compare_plane(self.write('a.fits', a), self.write('b.fits', b), comparison_mask=mask)
        self.assertEqual(report['primary_region'], 'above_elevation_cut')
        self.assertNotIn('celestial_valid', report['regions'])
        self.assertEqual(report['regions']['above_elevation_cut']['difference_rms'], 0)
        self.assertEqual(report['regions']['above_elevation_cut']['pixels'], 99)
        self.assertEqual(report['peak_shift_pixels'], 0)
        b[0, 0] = np.nan
        report, _ = compare_plane(self.root/'a.fits', self.write('nanoutside.fits', b), comparison_mask=mask)
        self.assertEqual(report['regions']['above_elevation_cut']['difference_rms'], 0)
        b[4, 4] = np.nan
        with self.assertRaisesRegex(ValueError, 'finite pixel masks'):
            compare_plane(self.root/'a.fits', self.write('naninside.fits', b), comparison_mask=mask)
        with self.assertRaisesRegex(ValueError, 'mask must be boolean'):
            compare_plane(self.root/'a.fits', self.root/'b.fits', comparison_mask=mask.astype(int))

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

    def test_calibrated_dysco_column_alias_requires_opt_in_and_matching_settings(self):
        reference = dict(data_column='CORRECTED_DATA', polarizations='V', threads=8)
        candidate = dict(reference, data_column='DATA')
        self.assertFalse(image_parameters_match(reference, candidate))
        self.assertTrue(image_parameters_match(reference, candidate, equivalent_calibrated_columns=True))
        self.assertFalse(image_parameters_match(candidate, reference, equivalent_calibrated_columns=True))
        self.assertFalse(image_parameters_match(reference, dict(candidate, threads=4),
                                                equivalent_calibrated_columns=True))
        self.assertFalse(image_parameters_match(reference, dict(candidate, data_column='MODEL_DATA'),
                                                equivalent_calibrated_columns=True))

    def test_run_mask_grid_and_binary_values_checked(self):
        import shutil
        a, b = self.root / 'a', self.root / 'b'
        a.mkdir(); b.mkdir()
        record = dict(status='complete', source_unchanged=True, parameters={}, container_id='same')
        image = self.write('image.fits', np.ones((10, 10)))
        for directory in (a, b):
            (directory/'complete.json').write_text(json.dumps(record))
            shutil.copyfile(image, directory/'image-image.fits')
        mask = self.write('mask.fits', np.ones((10,10)), ELCUT=10.)
        report = compare_run(a, b, mask)
        self.assertEqual(report['comparison_mask']['pixels'], 100)
        self.assertNotIn('data_columns', compare_run(a, b, mask, equivalent_calibrated_columns=True))
        invalid = self.write('badmask.fits', np.full((10,10), 2))
        with self.assertRaisesRegex(ValueError, 'only zero and one'):
            compare_run(a, b, invalid)
        shifted = self.write('shiftedmask.fits', np.ones((10,10)), CRVAL1=180.)
        with self.assertRaisesRegex(ValueError, 'mask celestial grid'):
            compare_run(a, b, shifted)

    def test_empty_and_zero_reference_regions_have_no_invented_ratios(self):
        a = np.zeros((2, 2))
        self.assertEqual(difference_metrics(a, a, np.zeros_like(a, dtype=bool)), {"pixels": 0})
        result = difference_metrics(a, np.ones_like(a), np.ones_like(a, dtype=bool))
        self.assertIsNone(result["relative_rms_difference"])
        self.assertEqual(result["difference_rms"], 1)

    def test_stokes_v_residual_normalization_and_zero_denominator(self):
        a,b=self.root/'reference-v',self.root/'candidate-v'
        a.mkdir();b.mkdir()
        record=dict(status='complete',source_unchanged=True,
                    parameters=dict(polarizations='V'),container_id='same')
        for directory in (a,b):
            (directory/'complete.json').write_text(json.dumps(record))
        stokes=dict(CTYPE3='STOKES',CRVAL3=4.,CRPIX3=1.,CDELT3=1.)
        reference=np.ones((1,10,10))
        candidate=reference+.25;candidate[0,0,0]=1000
        self.write('reference-v/image-V-image.fits',reference,**stokes)
        self.write('candidate-v/image-V-image.fits',candidate,**stokes)
        for directory in ('reference-v','candidate-v'):
            self.write(directory+'/image-V-residual.fits',np.full((1,10,10),.5),**stokes)
        mask=np.ones((10,10));mask[0,0]=0
        maskpath=self.write('mask.fits',mask,ELCUT=10.)
        report=compare_run(a,b,maskpath)
        ratio=report['residual_normalized']['image-V-image.fits']
        self.assertEqual(ratio['difference_rms'],.25)
        self.assertEqual(ratio['reference_residual_rms'],.5)
        self.assertEqual(ratio['difference_over_reference_residual_rms'],.5)
        (a/'image-V-residual.fits').unlink()
        self.write('reference-v/image-V-residual.fits',np.zeros((1,10,10)),**stokes)
        ratio=compare_run(a,b,maskpath)['residual_normalized']['image-V-image.fits']
        self.assertIsNone(ratio['difference_over_reference_residual_rms'])
        self.assertEqual(ratio['denominator_status'],'zero; ratio undefined')

    def test_stokes_identity_must_match_in_metadata_and_fits(self):
        a=np.ones((1,10,10))
        reference=self.write('v.fits',a,CTYPE3='STOKES',CRVAL3=4.)
        candidate=self.write('i.fits',a,CTYPE3='STOKES',CRVAL3=1.)
        with self.assertRaisesRegex(ValueError,'grid, coordinates, beam or units'):
            compare_plane(reference,candidate)
        for directory,pol in (('ref','V'),('candidate','I')):
            root=self.root/directory;root.mkdir()
            (root/'complete.json').write_text(json.dumps(dict(status='complete',source_unchanged=True,
                parameters=dict(polarizations=pol),container_id='same')))
        with self.assertRaisesRegex(ValueError,'settings'):
            compare_run(self.root/'ref',self.root/'candidate')


if __name__ == "__main__":
    unittest.main()
