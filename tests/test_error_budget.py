"""Relative thresholds must constrain stored reconstructions and verification."""
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import importlib.util
import shutil

import numpy as np
from casacore.tables import table, makearrcoldesc, maketabdesc

from vamp import codec
from vamp.archive import encode, decode, verify
from vamp.config import Config
from vamp.error_budget import component_signal_power, resolve_bounds
from test_inspection import make_ms, snapshot


class ErrorBudgetTests(unittest.TestCase):
    def config(self, **changes):
        cfg = replace(Config(mode='adaptive', codebook='power-midrise-v1',
            relative_bounds=dict(reference='block-rms-v1', max_component_rmse=.01,
                                 max_complex_error=.1, max_complex_bias=.01,
                                 max_baseline_component_rmse=.02)), **changes)
        cfg.validate()
        return cfg

    def test_analytic_reference_protection_zero_and_hybrid(self):
        cfg = self.config()
        a = np.array([[3+4j, 0j, complex(float('nan'), 1), 1e30j]], dtype='c8')
        f = np.array([[False, False, False, True]])
        self.assertEqual(component_signal_power(a, f), 6.25)
        b = resolve_bounds(a, f, cfg.bounds(), cfg.relative_bounds)
        self.assertAlmostEqual(b['max_component_mse'], .000625)
        self.assertAlmostEqual(b['max_complex_error'], .25)
        self.assertAlmostEqual(b['max_baseline_component_mse'], .0025)
        hybrid = replace(cfg, max_component_mse=.0001, max_complex_error=.5)
        h = resolve_bounds(a, f, hybrid.bounds(), hybrid.relative_bounds)
        self.assertEqual(h['max_component_mse'], .0001)
        self.assertAlmostEqual(h['max_complex_error'], .25)
        for values, flags in [(np.zeros_like(a), f), (a, np.ones_like(f))]:
            z = resolve_bounds(values, flags, cfg.bounds(), cfg.relative_bounds)
            self.assertTrue(all(v is None or v == 0 for v in z.values()))

    def test_scaling_native_portable_and_fallback(self):
        rng = np.random.default_rng(39)
        a = (rng.normal(size=(8, 32)) + 1j*rng.normal(size=(8, 32))).astype('c8')
        a[0,0] = 0j
        f = np.zeros_like(a, dtype=bool)
        cfg = self.config(compression='none', scale_quantiles=[1.])
        bits = []
        for factor in (2.**-20, 1., 2.**20):
            x = a*factor
            native = codec.encode_block(x, f, cfg, cfg.bounds())
            with patch.object(codec, '_native', None):
                portable = codec.encode_block(x, f, cfg, cfg.bounds())
            self.assertEqual(native.blob, portable.blob)
            limits = resolve_bounds(x, f, cfg.bounds(), cfg.relative_bounds)
            self.assertTrue(codec.within_bounds(native.metrics, limits))
            bits.append(native.bits)
        self.assertEqual(len(set(bits)), 1)
        cfg = replace(cfg, relative_bounds=dict(reference='block-rms-v1', max_component_rmse=0., max_complex_error=0.))
        exact = codec.encode_block(a, f, cfg, cfg.bounds())
        self.assertEqual(exact.bits, 0)
        self.assertEqual(exact.decoded.tobytes(), a.tobytes())

    def test_strict_configuration(self):
        for spec in ({}, {'reference':'unknown'}, {'reference':'block-rms-v1'},
                     {'reference':'block-rms-v1','max_component_rmse':True},
                     {'reference':'block-rms-v1','max_component_rmse':-1},
                     {'reference':'block-rms-v1','max_component_rmse':float('nan')},
                     {'reference':'block-rms-v1','max_component_rmse':.01},
                     {'reference':'block-rms-v1','max_complex_error':.1}):
            with self.assertRaises(ValueError):
                self.config(relative_bounds=spec)
        cfg = self.config(relative_bounds=dict(reference='block-rms-v1',max_component_rmse=1e308,max_complex_error=1.))
        a = np.full((1,1), 1e30j, dtype='c8')
        with self.assertRaisesRegex(ValueError, 'overflows'):
            codec.encode_block(a, np.zeros_like(a, bool), cfg, cfg.bounds())

    def test_archive_verifies_relative_blocks_and_keeps_exact_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); source = root/'input.ms'
            make_ms(source)
            before = snapshot(source)
            for codebook in ('uniform-symmetric-v1','power-midrise-v1'):
                cfg = self.config(codebook=codebook, baseline_group=2, channel_segment=2)
                path = root/(codebook+'.h5')
                encode(source, path, cfg)
                self.assertEqual(verify(path, source)['status'], 'passed')
                decode(path, root/(codebook+'.ms'))
            self.assertEqual(before, snapshot(source))

    def test_visibility_diagnostics_separate_autos_and_respect_product_order(self):
        spec = importlib.util.spec_from_file_location('analyze_visibility_errors',
            Path(__file__).resolve().parents[1]/'scripts/analyze_visibility_errors.py')
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = make_ms(root/'source.ms', pairs=((0,0),(0,1)), channels=3,
                             correlations=(12,10,9,11))  # YY,XY,XX,YX
            with table(str(source),readonly=False,ack=False) as t:
                t.addcols(maketabdesc(makearrcoldesc('CORRECTED_DATA',0j,ndim=2,valuetype='complex')))
                a=np.full((2,3,4),3+4j,dtype='c8');a[0]*=1000
                t.putcol('CORRECTED_DATA',a)
                flags=np.zeros(a.shape,bool);flags[1,0,2]=True;t.putcol('FLAG',flags)
            decoded=root/'decoded.ms';shutil.copytree(source,decoded)
            with table(str(decoded),readonly=False,ack=False) as t:
                b=a.copy();b[1,:,0]+=2;b[1,:,2]-=2;b[0]+=100
                t.putcol('CORRECTED_DATA',b)
            result=module.analyze(source,decoded)
            self.assertEqual(result['products']['cross/XX']['samples'],2)
            self.assertEqual(result['products']['cross/XX']['component_mse'],2.)
            self.assertEqual(result['products']['auto/XX']['component_mse'],5000.)
            self.assertEqual(result['products']['cross/I-paired']['samples'],2)
            self.assertEqual(result['products']['cross/I-paired']['component_mse'],0.)
            self.assertEqual(sum(x['samples'] for x in result['by_baseline_length']['cross/XX']),2)


if __name__ == '__main__': unittest.main()
