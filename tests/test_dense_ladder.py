"""Precision tags, three-byte code spans and independent archive decoding."""
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from varquanms import codec
from varquanms.archive import encode, decode, verify
from varquanms.config import Config, DEFAULT_BITS, SUPPORTED_BITS
from test_inspection import make_ms
from test_reference_decoder import reference


class DenseLadderTests(unittest.TestCase):
    def test_independent_packing_vectors_and_three_byte_spans(self):
        for bits, expected in ((5, 'e0073800'), (10, '00fc1f00800300'),
                               (14, '00c0ff1f0000800300')):
            codes = np.array([0, 2**bits-1, 1, 2**(bits-1), 3], dtype='u2')
            self.assertEqual(codec.pack_codes(codes, bits).hex(), expected)
            np.testing.assert_array_equal(codec.unpack_codes(bytes.fromhex(expected), bits, 5), codes)
        rng = np.random.default_rng(501014)
        for bits in SUPPORTED_BITS:
            for count in (1, 3, 7, 8, 9, 31, 257):
                values = rng.integers(0, 2**bits, count, dtype='u2')
                expected = sum(int(v) << (i*bits) for i, v in enumerate(values))
                blob = expected.to_bytes((count*bits+7)//8, 'little')
                self.assertEqual(codec.pack_codes(values, bits), blob)
                np.testing.assert_array_equal(codec.unpack_codes(blob, bits, count), values)
                if count*bits % 8:
                    damaged = blob[:-1]+bytes([blob[-1] | 0x80])
                    with self.assertRaisesRegex(ValueError, 'padding'):
                        codec.unpack_codes(damaged, bits, count)

    def test_nibble_precisions_portable_native_and_reference(self):
        rng = np.random.default_rng(101014)
        values = (rng.normal(size=(3, 31))+1j*rng.normal(size=(3, 31))).astype('c8')
        values.view('u4').flat[:4] = [0x7fc12345, 0x80000000, 0x7f800000, 0xff800000]
        flags = np.zeros(values.shape, bool); flags[1, 3] = True
        for name in ('uniform-symmetric-v1', 'uniform-midrise-v1', 'power-midrise-v1'):
            for bits in SUPPORTED_BITS:
                cfg = Config(mode='fixed', codebook=name, bits=[bits], compression='none',
                             block_header='precision-nibble-v1', max_component_mse=10.,
                             max_complex_error=100., scale_quantiles=[1.])
                native = codec.encode_block(values, flags, cfg, cfg.bounds())
                with patch.object(codec, '_native', None):
                    portable = codec.encode_block(values, flags, cfg, cfg.bounds())
                self.assertEqual(native.bits, bits)
                self.assertEqual(native.blob, portable.blob)
                self.assertEqual(native.blob[0], 0xA0 | codec.PRECISION_TAG_BITS.index(bits))
                self.assertEqual(codec.unpack_header(native.blob)[0], bits)
                self.assertEqual(reference.decode(native.blob, *values.shape), native.decoded.tobytes())
                self.assertEqual(codec.decode_block(native.blob, values.shape).tobytes(), native.decoded.tobytes())
                self.assertEqual(native.decoded[0,:2].tobytes(), values[0,:2].tobytes())
                for tag in range(11, 16):
                    bad = bytes([0xA0 | tag])+native.blob[1:]
                    for decoder in (lambda b: codec.decode_block(b, values.shape),
                                    lambda b: reference.decode(b, *values.shape)):
                        with self.assertRaisesRegex(ValueError, 'precision tag'):
                            decoder(bad)
                for bad in (native.blob[:7], native.blob[:-1]):
                    with self.assertRaises(ValueError):
                        codec.decode_block(bad, values.shape)

    def test_exact_tag_and_older_native_fallback(self):
        values = np.array([[.111+.123j, .234+.937j]], dtype='c8')
        flags = np.zeros(values.shape, bool)
        cfg = Config(mode='adaptive', bits=list(DEFAULT_BITS), block_header='precision-nibble-v1',
                     max_component_mse=0., max_complex_error=0.)
        exact = codec.encode_block(values, flags, cfg, cfg.bounds())
        self.assertEqual(exact.bits, 0); self.assertEqual(exact.blob[0], 0xA0)
        self.assertEqual(reference.decode(exact.blob, *values.shape), values.tobytes())
        old = SimpleNamespace(MIDRISE_V1=1, POWER_V1=1, THREE_BIT_V1=1)
        for name in ('uniform-symmetric-v1', 'uniform-midrise-v1', 'power-midrise-v1'):
            for bits in (5, 10, 14):
                trial = replace(cfg, codebook=name, mode='fixed', bits=[bits], max_component_mse=1., max_complex_error=1.)
                with patch.object(codec, '_native', None):
                    expected = codec.encode_block(values, flags, trial, trial.bounds())
                with patch.object(codec, '_native', old):
                    actual = codec.encode_block(values, flags, trial, trial.bounds())
                self.assertEqual(actual.blob, expected.blob)
        with self.assertRaisesRegex(ValueError, 'block_header'):
            replace(cfg, block_header='unknown').validate()

    def test_nibble_archive_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); source = make_ms(root/'input.ms', channels=17)
            for bits in (5, 10, 14):
                path = root/f'{bits}.h5'; restored = root/f'{bits}.ms'
                cfg = Config(mode='fixed', bits=[bits], baseline_group=2, channel_segment=8,
                             codebook='power-midrise-v1', block_header='precision-nibble-v1',
                             max_component_mse=1., max_complex_error=2.)
                encode(source, path, cfg); decode(path, restored)
                self.assertEqual(verify(path, source, restored)['status'], 'passed')


if __name__ == '__main__':
    unittest.main()
