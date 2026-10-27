from pathlib import Path
import tempfile
import unittest

from casacore.tables import table
import numpy as np

from vamp.archive import encode, decode, verify
from vamp.config import Config
from vamp.codec import unpack_codes
from vamp.weights import encode_chunk
from test_calibrated import calibrated_fixture


class WeightTests(unittest.TestCase):
    def test_three_value_dictionary_is_bitwise_exact_and_really_four_bit(self):
        bits = np.array([0, 0x80000000, 0x7fc00123, 0x3f800000, 0x7f800000], dtype='<u4')
        original = np.tile(bits, 5).view('<f4')
        book, codes, absolute, relative, lossy = encode_chunk(original, .01)
        self.assertEqual(len(codes), 13)
        result = book[unpack_codes(codes, 4, 25)]
        self.assertEqual(result.tobytes(), original.tobytes())
        self.assertEqual((absolute, relative, lossy), (0., 0., False))

    def test_relative_budget_zero_preservation_and_exact_fallback(self):
        original = np.linspace(1., 1.02, 101, dtype=np.float32)
        book, codes, _, relative, lossy = encode_chunk(original, .01)
        self.assertTrue(lossy); self.assertLess(relative, .01)
        self.assertEqual(len(codes), 51)
        self.assertIsNone(encode_chunk(original, 0))
        self.assertIsNone(encode_chunk(np.r_[np.float32(0), np.logspace(-5, 5, 32).astype(np.float32)], .01))

    def test_lossy_weight_round_trip_is_verified_against_original(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d); source = calibrated_fixture(root/'source.ms')
            # Weight cell sizes need not equal DATA correlations in this storage
            # fixture: the codec checks float32 array reconstruction, not science.
            with table(str(source), readonly=False, ack=False) as t:
                for row in range(3): t.putcell('WEIGHT', row, np.linspace(1., 1.02, 101, dtype=np.float32))
            cfg = Config(columns=['CORRECTED_DATA'], keep_calib_only=True)
            report = encode(source, root/'archive.h5', cfg)
            self.assertFalse(report['weights_quality']['WEIGHT']['original_bits_exact'])
            decode(root/'archive.h5', root/'restored.ms')
            result = verify(root/'archive.h5', source, root/'restored.ms')
            self.assertEqual(result['status'], 'passed')
            self.assertTrue(any(c.get('status') == 'weight_error_bound_passed' for c in result['column_checks']))


if __name__ == '__main__': unittest.main()
