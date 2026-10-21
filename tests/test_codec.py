from dataclasses import replace
import struct
import unittest
import zlib

import numpy as np

from varquanms.codec import HEADER, decode_block, encode_block, error_metrics, pack_codes, unpack_codes, within_bounds
from varquanms.config import Config
from varquanms.ordering import amplitude_statistics, all_orders, make_order


class CodecTests(unittest.TestCase):
    def config(self, **changes):
        result = Config(mode="adaptive", max_component_mse=1.0, max_complex_error=10.0,
                        compression="none", baseline_group=8, channel_segment=32)
        result = replace(result, **changes)
        result.validate()
        return result

    def test_every_packing_width_and_literal_fixtures(self):
        fixtures = [(2, [0, 1, 2, 3], "e4"), (4, [1, 2, 15], "210f"),
                    (6, [0, 1, 63, 32], "40f083"), (8, [1, 255], "01ff")]
        for bits, codes, expected in fixtures:
            self.assertEqual(pack_codes(np.asarray(codes), bits).hex(), expected)
        rng = np.random.default_rng(123)
        for bits in (2, 4, 6, 8):
            for size in (0, 1, 3, 7, 31, 192):
                codes = rng.integers(0, 2**bits, size=size, dtype=np.uint8)
                encoded = pack_codes(codes, bits)
                self.assertEqual(len(encoded), (size * bits + 7) // 8)
                np.testing.assert_array_equal(unpack_codes(encoded, bits, size), codes)
                with self.assertRaises(ValueError):
                    unpack_codes(encoded + b"\x00", bits, size)

    def test_corrupt_padding_reserved_codes_and_truncated_blocks(self):
        with self.assertRaisesRegex(ValueError, "padding"):
            unpack_codes(b"\xf1", 2, 1)
        bad = HEADER.pack(2, 0, 0, 1, 0, 1.0) + pack_codes(np.asarray([3, 1]), 2)
        with self.assertRaisesRegex(ValueError, "Reserved"):
            decode_block(bad, (1, 1))
        config = self.config()
        values = np.ones((8, 32), dtype=np.complex64)
        block = encode_block(values, np.zeros_like(values, dtype=bool), config, config.bounds())
        for truncated in (block.blob[:3], block.blob[:-1]):
            with self.assertRaises((ValueError, zlib.error)):
                decode_block(truncated, values.shape)

    def test_zero_and_all_flagged_blocks_are_exact(self):
        config = self.config()
        values = np.zeros((3, 17), dtype=np.complex64)
        for flag in (False, True):
            block = encode_block(values, np.full(values.shape, flag), config, config.bounds())
            self.assertEqual(decode_block(block.blob, values.shape).tobytes(), values.tobytes())

    def test_flags_nonfinite_signed_zero_and_rfi_exceptions(self):
        config = self.config(scale_quantiles=[0.95], bits=[4], max_component_mse=0.02, max_complex_error=0.5)
        rng = np.random.default_rng(45)
        values = (rng.normal(size=(8, 32)) + 1j * rng.normal(size=(8, 32))).astype(np.complex64)
        values.view(np.uint32).flat[0:4] = [0x7FC12345, 0x80000000, 0x7F800000, 0xFF800000]
        values[-1, -1] = 1e9 + 2e9j
        flags = np.zeros(values.shape, dtype=bool)
        flags[0, 2] = True
        block = encode_block(values, flags, config, config.bounds())
        decoded = decode_block(block.blob, values.shape)
        for row, channel in ((0, 0), (0, 1), (0, 2), (7, 31)):
            self.assertEqual(values[row, channel].tobytes(), decoded[row, channel].tobytes())
        self.assertTrue(within_bounds(block.metrics, config.bounds()))

    def test_exact_escape_when_nonrepresentable_values_fail_zero_budget(self):
        config = self.config(bits=[2], max_component_mse=0, max_complex_error=0)
        values = np.asarray([[0.11 + 0.03j, 0.213 + 0.987j, 0.357 + 0.891j]], dtype=np.complex64)
        block = encode_block(values, np.zeros(values.shape, dtype=bool), config, config.bounds())
        self.assertEqual(block.bits, 0)
        self.assertEqual(block.decoded.tobytes(), values.tobytes())

    def test_stochastic_streams_are_reproducible_and_product_independent(self):
        config = self.config(rounding="stochastic-pcg64-v1", bits=[4], scale_quantiles=[1.0])
        values = np.linspace(-1, 1, 256, dtype=np.float32).reshape(8, 32).astype(np.complex64) * (1 + 1j)
        flags = np.zeros(values.shape, dtype=bool)
        one = encode_block(values, flags, config, config.bounds(), stream=(0, 0, 0, 0))
        repeat = encode_block(values, flags, config, config.bounds(), stream=(0, 0, 0, 0))
        other = encode_block(values, flags, config, config.bounds(), stream=(0, 1, 0, 0))
        self.assertEqual(one.blob, repeat.blob)
        self.assertNotEqual(one.blob, other.blob)

    def test_sorting_keeps_spectra_and_classes_and_is_invertible(self):
        values = np.asarray([3, 0, 1, 3, 8], dtype=np.complex64)[:, None, None] * np.ones((1, 17, 2), dtype=np.complex64)
        flags = np.zeros(values.shape, dtype=bool)
        flags[4] = True
        pairs = np.asarray([[0, 1], [0, 0], [1, 2], [0, 2], [2, 2]])
        uvw = np.column_stack([np.arange(5), np.zeros((5, 2))])
        for method in ("original", "length", "amplitude", "buckets"):
            config = self.config(ordering=method, baseline_group=3)
            orders, _, counts = all_orders(values, flags, pairs, uvw, config)
            self.assertEqual(counts[4].tolist(), [0, 0])
            for product, order in enumerate(orders):
                restored = np.empty_like(values[:, :, product])
                restored[order["permutation"]] = values[order["permutation"], :, product]
                self.assertEqual(restored.tobytes(), values[:, :, product].tobytes())
                for lo, hi in zip(order["group_boundaries"][:-1], order["group_boundaries"][1:]):
                    selected = pairs[order["permutation"][lo:hi]]
                    self.assertTrue(np.all(selected[:, 0] == selected[:, 1]) or np.all(selected[:, 0] != selected[:, 1]))

    def test_strict_config_and_no_implicit_noise_estimation(self):
        for changes in ({"noise": {"component_variance": 1, "epsilon": .01}}, {"bits": [3]},
                        {"mode": "fixed"}, {"rounding": "dither"}, {"max_complex_error": None}):
            with self.assertRaises(ValueError):
                self.config(**changes)

    def test_buckets_with_extreme_reference_do_not_overflow(self):
        config = self.config(ordering="buckets", bucket_reference=1e-320)
        result = make_order(np.asarray([1e30, 1e-30, 0.0]), np.ones(3), np.zeros(3, dtype=bool), np.ones(3), config)
        self.assertEqual(result["permutation"].tolist(), [2, 1, 0])
        for bucket in result["buckets"]:
            if bucket["kind"] == "positive":
                self.assertTrue(np.isfinite(bucket["lower"]) and np.isfinite(bucket["upper"]))


if __name__ == "__main__":
    unittest.main()
