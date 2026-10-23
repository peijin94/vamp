from pathlib import Path
import tempfile
import unittest

from casacore.tables import table, makearrcoldesc, maketabdesc, tabledefinehypercolumn
import numpy as np

from varquanms.archive import encode, decode, verify, compare_ms
from varquanms.calibrated import KEYWORD, descriptor
from varquanms.cli import main
from varquanms.config import Config
from varquanms.container import Reader
from varquanms.partitioned import encode_partitioned, decode_partitioned
from test_inspection import make_ms, snapshot


def calibrated_fixture(path):
    make_ms(path)
    with table(str(path), readonly=False, ack=False) as t:
        desc = maketabdesc([makearrcoldesc(c, 0j, ndim=2, valuetype='complex')
                           for c in ('CORRECTED_DATA', 'MODEL_DATA')])
        tabledefinehypercolumn(desc, 'TiledModel', 3, ['MODEL_DATA'])
        t.addcols(desc)
        values = np.arange(18, dtype=np.float32).reshape(3, 3, 2).astype(np.complex64)
        values.imag = values.real * .125
        values[0, 0, 0] = complex(np.nan, np.inf)
        values[1, 1, 0] = complex(-0., -0.)
        t.putcol('CORRECTED_DATA', values)
        t.putcol('MODEL_DATA', values * 2)
        flags = t.getcol('FLAG'); flags[0, 0, 0] = True; t.putcol('FLAG', flags)
        t.putcell('FLAG_ROW', 2, True)
        weights = t.getcol('WEIGHT'); weights[1] = np.float32(.123456789); t.putcol('WEIGHT', weights)
    with table(str(path/'OBSERVATION'), readonly=False, ack=False) as t:
        t.addrows(1)
        t.putcell('LOG', 0, [])
        t.putcell('SCHEDULE', 0, [])
    return path


class CalibratedTests(unittest.TestCase):
    def test_both_hypercolumn_representations_drop_removed_column_references(self):
        hyper = {'ndim':3,'data':['DATA'],'coord':[],'id':[]}
        d = descriptor({'DATA':{},'_private_keywords_':{'Hypercolumn_Data':hyper},
                        '_define_hypercolumn_':{'Data':{'HCndim':3,'HCdatanames':['DATA']}}})
        self.assertEqual(d['_private_keywords_'], {})
        self.assertEqual(d['_define_hypercolumn_'], {})

    def test_profile_lossless_and_adaptive_preserve_flags_weights_and_unflagged_specials(self):
        for mode in ('lossless', 'adaptive'):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as d:
                root = Path(d); source = calibrated_fixture(root/'source.ms')
                before = snapshot(source)
                cfg = Config(mode=mode, columns=['CORRECTED_DATA'], keep_calib_only=True,
                             max_component_mse=.01, max_complex_error=.5, codebook='power-midrise-v1')
                result = encode(source, root/'archive.h5', cfg)
                self.assertFalse(result['exact_ms_round_trip'])
                self.assertEqual(result['omitted_visibility_columns'], ['DATA', 'MODEL_DATA'])
                decode(root/'archive.h5', root/'restored.ms', ms_v2_data=False)
                self.assertEqual(verify(root/'archive.h5', source, root/'restored.ms')['status'], 'passed')
                self.assertEqual(snapshot(source), before)
                with table(str(root/'restored.ms'), ack=False, readonly=True) as t:
                    self.assertNotIn('DATA', t.colnames()); self.assertNotIn('MODEL_DATA', t.colnames())
                    self.assertIn(KEYWORD, t.getkeywords())
                    self.assertNotIn('Hypercolumn_TiledModel', t.getdesc().get('_private_keywords_', {}))
                    values = t.getcol('CORRECTED_DATA')
                    flags = t.getcol('FLAG') | t.getcol('FLAG_ROW')[:, None, None]
                    self.assertTrue(np.all(values[flags].view(np.uint32) == 0))
                    self.assertTrue(np.signbit(values[1, 1, 0].real))
                with Reader(root/'archive.h5') as reader:
                    spec = reader.manifest['tables'][0]
                    self.assertNotIn('DATA', spec['columns'])
                    self.assertNotIn('MODEL_DATA', spec['columns'])
                decode(root/'archive.h5', root/'compatible.ms')
                self.assertEqual(verify(root/'archive.h5', source, root/'compatible.ms')['status'], 'passed')
                with table(str(root/'compatible.ms'), ack=False, readonly=True) as t:
                    self.assertEqual(t.getcol('DATA').tobytes(), t.getcol('CORRECTED_DATA').tobytes())

    def test_partition_bundle_and_cli_option(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); source = calibrated_fixture(root/'source.ms')
            with table(str(source), readonly=False, ack=False) as t:
                t.putcol('TIME', np.array([1., 1., 2.]))
            cfg = Config(columns=['CORRECTED_DATA'], keep_calib_only=True)
            encode_partitioned(source, root/'archive.bundle', cfg)
            decode_partitioned(root/'archive.bundle', root/'restored.ms')
            self.assertTrue(compare_ms(source, root/'restored.ms', keep_calib_only=True))
            decode_partitioned(root/'archive.bundle', root/'compatible.ms', ms_v2_data=True)
            checks = compare_ms(source, root/'compatible.ms', keep_calib_only=True)
            self.assertTrue(any(c.get('export') == 'derived_calibrated_copy_matches_CORRECTED_DATA_bits' for c in checks))
            self.assertEqual(main(['encode', str(source), str(root/'cli.h5'), '--keep-calib-only']), 0)
            with Reader(root/'cli.h5') as reader:
                self.assertTrue(reader.manifest['config']['keep_calib_only'])

    def test_conflicting_flags_and_missing_calibrated_column_fail_before_output(self):
        with self.assertRaisesRegex(ValueError, 'preserved flags'):
            Config(columns=['CORRECTED_DATA'], keep_calib_only=True, burn_flag=True).validate()
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); source = make_ms(root/'source.ms')
            with self.assertRaisesRegex(ValueError, 'CORRECTED_DATA'):
                encode(source, root/'bad.h5', Config(columns=['CORRECTED_DATA'], keep_calib_only=True))
            self.assertFalse((root/'bad.h5').exists())


if __name__ == '__main__': unittest.main()
