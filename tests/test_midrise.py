"""Full-range mapping conformance, zero bias, protection and archival tests."""
from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from varquanms import codec
from varquanms.archive import encode, decode, verify
from varquanms.config import Config, SUPPORTED_BITS
from test_inspection import make_ms, snapshot
from test_reference_decoder import reference


class MidriseTests(unittest.TestCase):
    def config(self, **changes):
        cfg = replace(Config(mode="adaptive", codebook="uniform-midrise-v1",
                             max_component_mse=1., max_complex_error=5.,
                             compression="none", scale_quantiles=[1.]), **changes)
        cfg.validate()
        return cfg

    def test_all_levels_literal_fixture_and_header_validation(self):
        # 4-bit q=[0,7,8,15], delta=1 -> [-7.5,-0.5,+0.5,+7.5].
        blob = bytes.fromhex("0400010002000000000000000000803f70f8")
        expected = bytes.fromhex("0000f0c0000000bf0000003f0000f040")
        self.assertEqual(codec.decode_block(blob, (1, 2)).tobytes(), expected)
        self.assertEqual(reference.decode(blob, 1, 2), expected)
        with self.assertRaisesRegex(ValueError, "codebook"):
            codec.decode_block(blob, (1, 2), "uniform-symmetric-v1")
        for mapping in (2, 256, 65535):
            bad = bytearray(blob); bad[2:4] = mapping.to_bytes(2, "little")
            with self.assertRaises(ValueError):
                codec.decode_block(bad, (1, 2))
            with self.assertRaises(ValueError):
                reference.decode(bad, 1, 2)
        for bad in (blob[:-1], blob + b"\x00", blob[:5]):
            with self.assertRaises(ValueError):
                codec.decode_block(bad, (1, 2))
        for bits in SUPPORTED_BITS:
            codes = np.arange(2**bits, dtype=np.uint16)
            blob = codec.HEADER.pack(bits, 0, 1, len(codes)//2, 0, 1.) + codec.pack_codes(codes, bits)
            values = codec.decode_block(blob, (1, len(codes)//2))
            np.testing.assert_array_equal(values.view(np.float32).ravel(), codes.astype(np.float32)-(2**bits-1)/2)
            self.assertEqual(reference.decode(blob, 1, len(codes)//2), values.tobytes())

    def test_rounding_ties_zero_bias_and_strict_bias_fallback(self):
        a = np.tile(np.asarray([-7.5-7.5j, -1+0j, 0+1j, 7.5+7.5j], dtype="c8"), (4, 8))
        cfg = self.config(bits=[4])
        b = codec.encode_block(a, np.zeros(a.shape, bool), cfg, cfg.bounds())
        self.assertEqual(b.bits, 4)
        np.testing.assert_array_equal(b.decoded[0, :4], [-7.5-7.5j, -1.5+.5j, .5+.5j, 7.5+7.5j])
        a = np.full((8, 32), 7.5+0j, dtype="c8")
        b = codec.encode_block(a, np.zeros(a.shape, bool), cfg, cfg.bounds())
        self.assertEqual(b.metrics["bias_imag"], .5)
        strict = replace(cfg, max_complex_bias=.1)
        b = codec.encode_block(a, np.zeros(a.shape, bool), strict, strict.bounds())
        self.assertEqual(b.bits, 0)
        self.assertEqual(b.decoded.tobytes(), a.tobytes())

    def test_special_values_zero_blocks_overflow_and_stochastic_streams(self):
        rng = np.random.default_rng(34)
        a = (rng.normal(size=(8, 31)) + 1j*rng.normal(size=(8, 31))).astype("c8")
        a.view(np.uint32).flat[:6] = [0x80000000, 0, 0x7fc12345, 0xff800000, 0x7f800000, 0]
        a[-1, -1] = 1e9+2e9j
        flags = np.zeros(a.shape, bool); flags[0, 3] = True
        for bits in SUPPORTED_BITS:
            cfg = self.config(bits=[bits], scale_quantiles=[.95])
            block = codec.encode_block(a, flags, cfg, cfg.bounds())
            self.assertEqual(block.bits, bits)
            self.assertEqual(block.decoded[0, :4].tobytes(), a[0, :4].tobytes())
            self.assertEqual(block.decoded[-1, -1].tobytes(), a[-1, -1].tobytes())
            self.assertEqual(reference.decode(block.blob, *a.shape), block.decoded.tobytes())
        zeros = np.zeros((3, 17), dtype="c8")
        for flagged in (False, True):
            cfg = self.config()
            block = codec.encode_block(zeros, np.full(zeros.shape, flagged), cfg, cfg.bounds())
            self.assertEqual(block.decoded.tobytes(), zeros.tobytes())
        a = np.full((8, 32), np.finfo(np.float32).max, dtype="c8")
        cfg = self.config()
        with np.errstate(over="ignore"):
            block = codec.encode_block(a, np.zeros(a.shape, bool), cfg, cfg.bounds())
        self.assertEqual(block.bits, 0)
        cfg = self.config(rounding="stochastic-pcg64-v1", bits=[4])
        block = codec.encode_block(zeros, np.zeros(zeros.shape, bool), cfg, cfg.bounds())
        self.assertEqual(block.decoded.tobytes(), zeros.tobytes())
        a = np.full((64, 192), 7.5+0j, dtype="c8")
        def run(stream):
            return codec.encode_block(a, np.zeros(a.shape, bool), cfg, cfg.bounds(), stream=stream)
        one, repeat, other = run((0,)), run((0,)), run((1,))
        self.assertEqual(one.blob, repeat.blob)
        self.assertNotEqual(one.blob, other.blob)
        self.assertLess(abs(one.metrics["bias_imag"]), .02)

    @unittest.skipUnless(getattr(codec._native, "MIDRISE_V1", 0), "new native extension unavailable")
    def test_native_reference_equivalence_at_all_widths_and_bounds(self):
        rng = np.random.default_rng(81)
        for factor in (0., 1e-38, .001, 1., 1e25):
            a = ((rng.normal(size=(9, 31)) + 1j*rng.normal(size=(9, 31)))*factor).astype("c8")
            a.view(np.uint32).flat[:6] = [0x80000000, 0, 0x7fc12345, 0xff800000, 0x7f800000, 0]
            flags = rng.random(a.shape) < .1
            for bits in (*[(b,) for b in SUPPORTED_BITS], (3,4,6,8,12,16), (4,8,12,16)):
                cfg = self.config(bits=list(bits), compression="gzip", scale_quantiles=[1., .95],
                                  max_component_mse=max(1e-80, factor**2*.01), max_complex_error=max(1e-40, factor*.5))
                for p99 in (None, cfg.max_complex_error):
                    config = replace(cfg, max_p99_error=p99)
                    native = codec.encode_block(a, flags, config, config.bounds())
                    with patch.object(codec, "_native", None):
                        portable = codec.encode_block(a, flags, config, config.bounds())
                    self.assertEqual(native.blob, portable.blob)
                    self.assertEqual(native.decoded.tobytes(), portable.decoded.tobytes())
                    self.assertEqual(native.metrics, portable.metrics)

    def test_archive_export_protection_and_old_native_fallback(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); source = make_ms(root/"source.ms", channels=17)
            before = snapshot(source)
            cfg = self.config(baseline_group=2, channel_segment=8, workers=2)
            encode(source, root/"archive.h5", cfg)
            decode(root/"archive.h5", root/"restored.ms")
            with patch.object(codec, "_native", None):
                self.assertEqual(verify(root/"archive.h5", source, root/"restored.ms")["status"], "passed")
            self.assertEqual(before, snapshot(source))
        # A stale extension must never receive the extended native signature.
        with patch.object(codec, "_native", object()):
            a = np.ones((2, 17), dtype="c8")
            cfg = self.config()
            b = codec.encode_block(a, np.zeros(a.shape, bool), cfg, cfg.bounds())
            self.assertEqual(codec.decode_block(b.blob, a.shape).tobytes(), b.decoded.tobytes())
