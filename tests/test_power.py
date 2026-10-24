"""Power-codec conformance, nonlinear selection, protection and archival checks."""
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from varquanms import codec, power
from varquanms.archive import encode, decode, verify
from varquanms.config import Config, SUPPORTED_BITS
from varquanms.container import Reader
from test_inspection import make_ms, snapshot
from test_reference_decoder import reference


class PowerTests(unittest.TestCase):
    def config(self, **changes):
        cfg = replace(Config(mode="adaptive", codebook="power-midrise-v1", compression="none",
                             max_component_mse=1., max_complex_error=10., scale_quantiles=[1.]), **changes)
        cfg.validate(); return cfg

    def test_literal_levels_and_gamma_direction(self):
        # b=4, t=225, medians=(10,-20), gamma=1/2, codes 0,7,8,15.
        # u=(-1,-1/15,+1/15,+1), residual=(-225,-1,+1,+225).
        blob = bytes.fromhex("04000200020000000000000000006143" "000020410000a0c10000003f" "70f8")
        expected = np.asarray([[-215-21j,11+205j]], dtype="<c8")
        self.assertEqual(codec.decode_block(blob, (1,2)).tobytes(), expected.tobytes())
        self.assertEqual(reference.decode(blob, 1, 2), expected.tobytes())
        constant = codec.HEADER.pack(4,0,2,1,0,0.)+power.PARAMETERS.pack(-0.,7.,.5)+b'\x88'
        expected_bits = bytes.fromhex('000000800000e040')
        self.assertEqual(codec.decode_block(constant,(1,1)).tobytes(),expected_bits)
        self.assertEqual(reference.decode(constant,1,1),expected_bits)
        for bits in SUPPORTED_BITS:
            codes = np.arange(2**bits, dtype=np.uint16)
            old_gap = None
            for gamma in (1., .5, .25):
                levels = power.reconstruct(codes, bits, 1., 0., gamma)
                self.assertTrue(np.all(levels[1:] > levels[:-1]))
                self.assertEqual(levels[0], -1); self.assertEqual(levels[-1], 1)
                gap = levels[2**(bits-1)]-levels[2**(bits-1)-1]
                if old_gap is not None: self.assertLess(gap, old_gap)
                old_gap = gap
                blob = codec.HEADER.pack(bits,0,2,len(codes)//2,0,1.)+power.PARAMETERS.pack(0,0,gamma)+codec.pack_codes(codes,bits)
                self.assertEqual(reference.decode(blob,1,len(codes)//2),codec.decode_block(blob,(1,len(codes)//2)).tobytes())

    def test_nearest_is_in_visibility_space_and_center_storage_is_simulated(self):
        rng = np.random.default_rng(12)
        a = (10+rng.uniform(-1,1,(8,31)) + 1j*(-20+rng.uniform(-1,1,(8,31)))).astype('c8')
        flags = np.zeros(a.shape,bool)
        for gamma in (1., .5, .25):
            cfg = self.config(bits=[4],power_gammas=[gamma])
            b = codec.encode_block(a, flags,cfg,cfg.bounds())
            self.assertEqual(b.bits,4)
            par = b.parameters
            components = np.stack([a.real,a.imag],-1).astype('f8')
            center = np.median(components,axis=(0,1)).astype('f4')
            np.testing.assert_array_equal(center,[par['center_real'],par['center_imag']])
            for c in range(2):
                levels = power.reconstruct(np.arange(16),4,par['radius'],center[c],gamma)
                best = np.min(np.abs(components[...,c,None]-levels),axis=-1)
                actual = b.decoded.real if c==0 else b.decoded.imag
                np.testing.assert_array_equal(np.abs(actual.astype('f8')-components[...,c]),best)
        cfg = self.config(power_gammas=[1.,.5,.25])
        dense = ((rng.standard_normal((16,64))**3)*.001 + 3 + 1j*5).astype('c8')
        cfg = replace(cfg,bits=[4],max_component_mse=.0001,max_complex_error=.1)
        b = codec.encode_block(dense,np.zeros(dense.shape,bool),cfg,cfg.bounds())
        self.assertTrue(codec.within_bounds(b.metrics,cfg.bounds()))

    def test_exact_constants_protection_exceptions_and_lossless_fallback(self):
        cfg = self.config(bits=[4],scale_quantiles=[.95])
        rng=np.random.default_rng(57)
        a=(10+rng.normal(size=(8,31))+1j*rng.normal(size=(8,31))).astype('c8')
        a.view('u4').flat[:6]=[0x80000000,0,0x7fc12345,0xff800000,0x7f800000,0]
        a[-1,-1]=1e30+1e30j
        flags=np.zeros(a.shape,bool);flags[0,3]=True
        b=codec.encode_block(a,flags,cfg,cfg.bounds())
        self.assertEqual(b.bits,4)
        self.assertEqual(b.decoded[0,:4].tobytes(),a[0,:4].tobytes())
        self.assertEqual(b.decoded[-1,-1].tobytes(),a[-1,-1].tobytes())
        for value in (0j,3+7j,complex(np.finfo('f4').max,0)):
            a=np.full((3,17),value,dtype='c8')
            for flagged in (False,True):
                b=codec.encode_block(a,np.full(a.shape,flagged),cfg,cfg.bounds())
                self.assertEqual(b.decoded.tobytes(),a.tobytes())
        a=(rng.normal(size=(8,31))+1j*rng.normal(size=(8,31))).astype('c8')
        strict=replace(cfg,max_component_mse=0,max_complex_error=0)
        b=codec.encode_block(a,np.zeros(a.shape,bool),strict,strict.bounds())
        self.assertEqual(b.bits,0);self.assertEqual(b.decoded.tobytes(),a.tobytes())
        a[0,:2]=[-np.finfo('f4').max,np.finfo('f4').max]
        with np.errstate(over='ignore',invalid='ignore'):
            b=codec.encode_block(a,np.zeros(a.shape,bool),cfg,cfg.bounds())
        self.assertTrue(codec.within_bounds(b.metrics,cfg.bounds()))

    def test_corruption_unknown_config_and_truncation(self):
        for options in ({'power_gammas':[.75]}, {'power_gammas':[True]}, {'power_gammas':[]},
                        {'power_center':'mean'}, {'rounding':'stochastic-pcg64-v1'},
                        {'power_gamma_policy':'unknown'}, {'power_ratio_reference':0},
                        {'power_ratio_reference':float('nan')}, {'power_ratio_reference':True}):
            with self.assertRaises(ValueError):self.config(**options)
        blob=codec.HEADER.pack(4,0,2,2,0,1.)+power.PARAMETERS.pack(0,0,.5)+bytes.fromhex('70f8')
        for bad in (blob[:-1],blob+b'\0',blob[:25],blob[:24]+power.PARAMETERS.pack(.75,0,0)[:4]+blob[28:]):
            with self.assertRaises(ValueError):codec.decode_block(bad,(1,2))
        with self.assertRaisesRegex(ValueError,'codebook'):
            codec.decode_block(blob,(1,2),'uniform-midrise-v1')
        with self.assertRaises(ValueError):
            codec.decode_block(codec.HEADER.pack(0,0,2,1,0,0.)+bytes(20),(1,1))

    def test_ratio_gamma_statistics_boundaries_and_protection(self):
        residuals = np.array([[-1., -2.], [1., 2.]])
        s = power.maximum_std_ratio(residuals)
        self.assertAlmostEqual(s['component_std'], np.sqrt(2.5))
        self.assertAlmostEqual(s['max_std_ratio'], 2/np.sqrt(2.5))
        for factor in (1e-30, 1e30):
            self.assertAlmostEqual(power.maximum_std_ratio(residuals*factor)['max_std_ratio'], s['max_std_ratio'])
        for ratio, gamma in [(None,1.),(4.,1.),(np.nextafter(4.,np.inf),.5),(8.,.5),(9.,.25)]:
            self.assertEqual(power.gamma_from_ratio(ratio,4.,power.GAMMAS),gamma)
        self.assertEqual(power.gamma_from_ratio(20.,4.,[1.,.5]),.5)
        self.assertIsNone(power.maximum_std_ratio(np.empty((0,2)))['max_std_ratio'])
        self.assertIsNone(power.maximum_std_ratio(np.ones((4,2)))['max_std_ratio'])
        # Excluding protected values keeps even enormous flagged samples from
        # changing the selector; inlier mode also excludes range exceptions.
        a = np.tile(np.linspace(-1,1,64,dtype='f4'),(8,1)).astype('c8')
        a.imag = a.real
        flags = np.zeros(a.shape,bool)
        a[0,0] = 1e20+1e20j
        cfg = self.config(power_gammas=[1.,.5,.25],power_gamma_policy='max-std-v1')
        self.assertEqual(set(power._parameters(a,flags,cfg)[4][:,3]),{.25})
        flags[0,0] = True
        self.assertEqual(set(power._parameters(a,flags,cfg)[4][:,3]),{1.})
        flags[0,0] = False
        trimmed = replace(cfg,power_gamma_policy='max-std-inlier-v1',scale_quantiles=[.99])
        self.assertEqual(set(power._parameters(a,flags,trimmed)[4][:,3]),{1.})
        all_flagged = np.ones_like(flags)
        self.assertEqual(set(power._parameters(a,all_flagged,trimmed)[4][:,3]),{1.})
        with self.assertRaisesRegex(ValueError,'requires power-midrise'):
            self.config(codebook='uniform-symmetric-v1',power_gamma_policy='max-std-v1')

    def test_ratio_policy_preserves_error_checks_and_decoder_equivalence(self):
        rng = np.random.default_rng(91)
        a = (rng.normal(size=(8,31))+1j*rng.normal(size=(8,31))).astype('c8')
        a[-1,-1] = 1000+1000j
        flags = np.zeros(a.shape,bool)
        for policy in ('max-std-v1','max-std-inlier-v1'):
            cfg = self.config(power_gamma_policy=policy,power_gammas=[1.,.5,.25],
                              scale_quantiles=[1.,.99],compression='gzip')
            selected = codec.encode_block(a,flags,cfg,cfg.bounds())
            with patch.object(codec,'_native',None):
                portable = codec.encode_block(a,flags,cfg,cfg.bounds())
            self.assertEqual(selected.blob,portable.blob)
            self.assertEqual(reference.decode(selected.blob,*a.shape),selected.decoded.tobytes())
            self.assertTrue(codec.within_bounds(selected.metrics,cfg.bounds()))
            searched = replace(cfg,power_gamma_policy='search')
            full = codec.encode_block(a,flags,searched,searched.bounds())
            self.assertLessEqual(len(full.blob),len(selected.blob))
            exact = replace(cfg,max_component_mse=0,max_complex_error=0)
            escaped = codec.encode_block(a,flags,exact,exact.bounds())
            self.assertEqual(escaped.bits,0)
            self.assertEqual(escaped.decoded.tobytes(),a.tobytes())

    @unittest.skipUnless(getattr(codec._native,'POWER_V1',0),'native power extension unavailable')
    def test_native_portable_equivalence_and_near_thresholds(self):
        rng=np.random.default_rng(9)
        for factor in (0.,1e-38,.001,1.,1e25):
            a=((rng.normal(size=(5,31))+1j*rng.normal(size=(5,31)))*factor).astype('c8')
            a.view('u4').flat[:6]=[0x80000000,0,0x7fc12345,0xff800000,0x7f800000,0]
            flags=rng.random(a.shape)<.1
            for bits in SUPPORTED_BITS:
                cfg=self.config(bits=[bits],compression='gzip',scale_quantiles=[1.,.95],power_gammas=[1.,.5,.25],
                                max_component_mse=max(1e-80,factor**2*.1),max_complex_error=max(1e-40,factor))
                n=codec.encode_block(a,flags,cfg,cfg.bounds())
                with patch.object(codec,'_native',None):p=codec.encode_block(a,flags,cfg,cfg.bounds())
                self.assertEqual(n.blob,p.blob);self.assertEqual(n.metrics,p.metrics)
        a=(rng.normal(size=(8,32))+1j*rng.normal(size=(8,32))).astype('c8');flags=np.zeros(a.shape,bool)
        cfg=self.config(bits=[8]);base=codec.encode_block(a,flags,cfg,cfg.bounds())
        for bound,metric in (('max_component_mse','component_mse'),('max_complex_bias','complex_bias_abs'),
                             ('max_complex_error','max_complex_error'),('max_baseline_component_mse','max_baseline_component_mse'),('max_p99_error','p99_error')):
            for direction in (0.,np.inf):
                c=replace(cfg,**{bound:float(np.nextafter(base.metrics[metric],direction))})
                n=codec.encode_block(a,flags,c,c.bounds())
                with patch.object(codec,'_native',None):p=codec.encode_block(a,flags,c,c.bounds())
                self.assertEqual(n.blob,p.blob);self.assertTrue(codec.within_bounds(n.metrics,c.bounds()))

    def test_archive_export_and_independent_decoder(self):
        for policy in ('search','max-std-v1','max-std-inlier-v1'):
            with self.subTest(policy=policy), tempfile.TemporaryDirectory() as d:
                root=Path(d);source=make_ms(root/'source.ms',channels=17);before=snapshot(source)
                cfg=self.config(baseline_group=2,channel_segment=8,workers=2,power_gamma_policy=policy)
                encode(source,root/'archive.h5',cfg);decode(root/'archive.h5',root/'restored.ms')
                with patch.object(codec,'_native',None):
                    self.assertEqual(verify(root/'archive.h5',source,root/'restored.ms')['status'],'passed')
                with Reader(root/'archive.h5') as reader:
                    self.assertEqual(reader.manifest['config']['power_gamma_policy'],policy)
                self.assertEqual(before,snapshot(source))
        a=np.full((2,17),3+5j,dtype='c8');cfg=self.config()
        with patch.object(codec,'_native',object()):b=codec.encode_block(a,np.zeros(a.shape,bool),cfg,cfg.bounds())
        self.assertEqual(reference.decode(b.blob,*a.shape),a.tobytes())
