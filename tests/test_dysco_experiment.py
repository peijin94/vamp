"""Guard the supported layout and byte accounting of the Dysco benchmark."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np


@unittest.skipUnless(importlib.util.find_spec('casacore'), 'python-casacore is required')
class DyscoExperimentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location(
            'dysco_encode', Path(__file__).resolve().parents[1] / 'scripts' / 'dysco_encode.py')
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
        import run_dysco_benchmark
        cls.case_name = staticmethod(run_dysco_benchmark.case_name)
        from verify_dysco import flagged_value_counts
        cls.flagged_value_counts = staticmethod(flagged_value_counts)
        from report_dysco_comparison import image_metrics
        cls.image_metrics = staticmethod(image_metrics)
        import prepare_dysco_input_identity
        cls.identity = prepare_dysco_input_identity

    def source(self, **changes):
        columns = dict(TIME=np.array([1., 1., 2., 2.]),
                       ANTENNA1=np.array([0, 1, 0, 1]),
                       ANTENNA2=np.array([1, 2, 1, 2]),
                       DATA_DESC_ID=np.array([0, 0, 0, 0]),
                       FIELD_ID=np.array([0, 0, 0, 0]))
        columns.update(changes)

        class Source:
            def getcol(self, name, start=0, count=None):
                array = columns[name]
                return array[start:] if count is None else array[start:start + count]
        return Source()

    def test_time_blocks_and_order_are_required(self):
        self.assertEqual(self.module.regular_layout(self.source())['rows_per_time_block'], 2)
        with self.assertRaisesRegex(ValueError, 'at least two'):
            self.module.regular_layout(self.source(TIME=np.ones(4)))
        with self.assertRaisesRegex(ValueError, 'ordered'):
            self.module.regular_layout(self.source(TIME=np.array([2., 2., 1., 1.])))
        with self.assertRaisesRegex(ValueError, 'ordering differs in ANTENNA2'):
            self.module.regular_layout(self.source(ANTENNA2=np.array([1, 2, 2, 1])))
        with self.assertRaisesRegex(ValueError, 'ordering differs in FIELD_ID'):
            self.module.regular_layout(self.source(FIELD_ID=np.array([0, 0, 0, 1])))

    def test_data_manager_files_are_counted_without_other_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            (directory / 'table.f6').write_bytes(b'abc')
            (directory / 'table.f6_1').write_bytes(b'de')
            (directory / 'table.f60').write_bytes(b'fghi')
            (directory / 'SUB').mkdir()
            (directory / 'SUB' / 'table.f6').write_bytes(b'j')
            result = self.module.inventory(directory, 6)
            self.assertEqual(result['dysco_data_manager_bytes'], 5)
            self.assertEqual(result['total_bytes'], 10)
            self.assertEqual(set(result['dysco_data_manager_files']), {'table.f6', 'table.f6_1'})

    def test_truncation_variants_do_not_overwrite_default_runs(self):
        self.assertEqual(self.case_name('RF', 3, 2.5), 'dysco-rf-03')
        self.assertEqual(self.case_name('RF', 3, 1.5), 'dysco-rf-03-t15')
        self.assertEqual(self.case_name('AF', 8, 3.5), 'dysco-af-08-t35')
        self.assertEqual(self.case_name('RF', 3, 1.5, 'nan'), 'dysco-rf-03-t15-nan')

    def test_flagged_value_marker_policies_are_checked_separately(self):
        zero = np.array([complex(0., 0.), complex(-0., 0.), complex(1., 0.)], dtype='c8')
        self.assertEqual(self.flagged_value_counts(zero, 'positive complex zero'), (1, 2))
        nan = np.array([complex(np.nan, np.nan), complex(np.nan, 0.), complex(0., 0.)], dtype='c8')
        self.assertEqual(self.flagged_value_counts(nan, 'complex NaN'), (2, 0))
        with self.assertRaisesRegex(ValueError, 'Unknown'):
            self.flagged_value_counts(zero, 'ignore')

    def test_unconverged_clean_is_measured_but_not_reported_as_converged(self):
        comparison = {'status':'measured', 'final_threshold_reached':False,
                      'residual_normalized':{'image-image.fits':{'denominator_status':'positive',
                          'difference_rms':.5, 'reference_residual_rms':2.,
                          'difference_over_reference_residual_rms':.25}},
                      'products':{'image-dirty.fits':{'regions':{'above_elevation_cut':{'difference_rms':.1}}}}}
        result = self.image_metrics(comparison)
        self.assertEqual(result['difference_over_reference_residual_percent'],25.)
        self.assertFalse(result['clean_final_threshold_reached'])

    def test_input_identity_rejects_lossy_or_changed_prepared_ms(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); original=root/'original'; prepared=root/'prepared'
            original.mkdir(); prepared.mkdir(); archive=root/'archive.h5'; archive.write_bytes(b'archive')
            verification=root/'verification.json'; report=root/'identity.json'
            proof={'status':'passed','keep_calib_only':True,'source_files_unchanged':True,
                   'source':str(original),'restored':str(prepared),'source_file_sha256':{'x':'source'},
                   'export':{'derived_calibrated_DATA':True},
                   'quantized_columns_match_archive_decoder_bitwise':['CORRECTED_DATA'],
                   'column_checks':[{'table':'.','column':'CORRECTED_DATA','export':'bitwise_exact'},
                                    {'table':'.','column':'DATA','export':'derived_calibrated_copy_matches_CORRECTED_DATA_bits'}]}
            verification.write_text(json.dumps(proof))
            argv=['prepare',str(original),str(prepared),'--archive',str(archive),
                  '--verification',str(verification),'--report',str(report)]
            with patch.object(sys,'argv',argv):
                with self.assertRaisesRegex(ValueError,'lossless calibrated-only'):
                    self.identity.main()
            proof['quantized_columns_match_archive_decoder_bitwise']=[]
            verification.write_text(json.dumps(proof))
            manifest={'config':{'mode':'lossless','keep_calib_only':True}}
            hashes=[{'x':'source'},{'x':'prepared'},{'x':'source'},{'x':'changed'}]
            with patch.object(sys,'argv',argv), patch.object(self.identity,'snapshot',side_effect=hashes), \
                    patch.object(self.identity,'_load',return_value=(archive,manifest)), \
                    patch.object(self.identity,'compare_ms',return_value=[]):
                with self.assertRaisesRegex(ValueError,'changed during'):
                    self.identity.main()
            self.assertFalse(report.exists())


if __name__ == '__main__':
    unittest.main()
