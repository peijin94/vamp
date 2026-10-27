import importlib.util
import json
from pathlib import Path
import unittest

import numpy as np

from vamp.codec import decode_block, encode_block, pack_codes
from vamp.config import Config, SUPPORTED_BITS

root = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("reference_block_decoder", root / "scripts/reference_block_decoder.py")
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)


class ReferenceDecoderTests(unittest.TestCase):
    def test_literal_conformance_and_independent_decoder(self):
        fixtures = json.loads((root / "tests/fixtures/codec-v1.json").read_text())
        tagged = json.loads((root / "tests/fixtures/precision-nibble-v1.json").read_text())
        for key in ('packing', 'blocks'):
            fixtures[key].extend(tagged[key])
        for fixture in fixtures["packing"]:
            self.assertEqual(pack_codes(np.asarray(fixture["codes"]), fixture["bits"]).hex(), fixture["hex"])
        for fixture in fixtures["blocks"]:
            blob = bytes.fromhex(fixture["block_hex"])
            expected = bytes.fromhex(fixture["complex64_hex"])
            self.assertEqual(reference.decode(blob, fixture["rows"], fixture["channels"]), expected)
            self.assertEqual(decode_block(blob, (fixture["rows"], fixture["channels"])).tobytes(), expected)
        rng = np.random.default_rng(42)
        values = (rng.normal(size=(8, 31)) + 1j * rng.normal(size=(8, 31))).astype(np.complex64)
        values.view(np.uint32).flat[0:2] = [0x7FC12345, 0x80000000]
        flags = np.zeros(values.shape, dtype=bool)
        flags[2, 2] = True
        for bits in SUPPORTED_BITS:
            config = Config(mode="fixed", bits=[bits], max_component_mse=1.0, max_complex_error=5.0)
            block = encode_block(values, flags, config, config.bounds())
            self.assertEqual(block.bits, bits)
            self.assertEqual(reference.decode(block.blob, *values.shape), block.decoded.tobytes())


if __name__ == "__main__":
    unittest.main()
