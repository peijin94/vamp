from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

from casacore.tables import makearrcoldesc, makescacoldesc, maketabdesc, table
import h5py
import numpy as np

from vamp.archive import decode, encode, verify
from vamp.config import Config, SUPPORTED_BITS
from vamp.container import Reader
from vamp.visibility import read_visibility, read_visibility_block
from test_inspection import make_ms, snapshot


class QuantizedArchiveTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def source(self, special=True):
        pairs = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))
        source = make_ms(self.root / "source.ms", pairs=pairs, channels=17)
        rng = np.random.default_rng(876)
        values = (rng.normal(size=(6, 17, 2)) + 1j * rng.normal(size=(6, 17, 2))).astype(np.complex64)
        values[:, :, 0] *= np.float32(100)
        values[:, :, 1] *= np.float32(.01)
        flags = np.zeros(values.shape, dtype=bool)
        if special:
            values.view(np.uint32).flat[0:4] = [0x7FC12345, 0x80000000, 0x7F800000, 0xFF800000]
            values[1, 2] = np.complex64(1e6 + 2e6j)
            flags[2] = True
        with table(str(source), readonly=False, ack=False) as t:
            t.putcol("DATA", values)
            t.putcol("FLAG", flags)
        return source

    def config(self, **overrides):
        config = Config(mode="adaptive", baseline_group=3, channel_segment=16,
                        max_component_mse=.02, max_complex_error=.8, scale_quantiles=[1.0, .95])
        return replace(config, **overrides)

    def test_adaptive_orderings_and_exact_preserved_columns(self):
        source = self.source()
        before = snapshot(source)
        for index, method in enumerate(("original", "class-original", "length", "orientation", "length-orientation", "antenna", "amplitude", "buckets")):
            archive, restored = self.root / f"{index}.h5", self.root / f"{index}.ms"
            summary = encode(source, archive, self.config(ordering=method))
            self.assertEqual(sum(summary["stored_bytes_by_category"].values()), archive.stat().st_size)
            self.assertEqual(verify(archive, source)["status"], "passed")
            decode(archive, restored)
            self.assertEqual(verify(archive, source, restored)["status"], "passed")
        self.assertEqual(before, snapshot(source))

    def test_all_fixed_bit_widths_and_independent_last_block(self):
        source = self.source(special=False)
        for bits in SUPPORTED_BITS:
            archive = self.root / f"fixed{bits}.h5"
            encode(source, archive, self.config(mode="fixed", bits=[bits]))
            with Reader(archive) as reader:
                spec = reader.manifest["tables"][0]["columns"]["DATA"]
                product = spec["products"][1]
                rows, lo, hi, block = read_visibility_block(reader, spec, 1, product["block_count"] - 1)
                all_values = read_visibility(reader, spec)
                np.testing.assert_array_equal(block, all_values[rows, lo:hi, 1])
                self.assertTrue(set(product["component_bit_width_histogram"]).issubset({str(bits), "float32"}))
            self.assertEqual(verify(archive, source)["status"], "passed")

    def test_burn_flag_preserves_same_quantized_values_and_protected_samples(self):
        source = self.source()
        with table(str(source), readonly=False, ack=False) as t:
            t.putcell("FLAG_ROW", 3, True)
        original = snapshot(source)
        decoded = []
        summaries = []
        for burn in (False, True):
            archive, output = self.root / f"burn-{burn}.h5", self.root / f"burn-{burn}.ms"
            summary = encode(source, archive, self.config(burn_flag=burn))
            summaries.append(summary["visibility_quality"])
            decode(archive, output)
            self.assertEqual(verify(archive, source, output)["status"], "passed")
            with table(str(output), readonly=True, lockoptions="usernoread", ack=False) as t:
                decoded.append(t.getcol("DATA").tobytes())
                if burn:
                    self.assertFalse(t.getcol("FLAG").any())
                    self.assertFalse(t.getcol("FLAG_ROW").any())
        self.assertEqual(decoded[0], decoded[1])
        self.assertEqual(summaries[0], summaries[1])
        self.assertEqual(snapshot(source), original)

    def test_zero_budget_and_shared_stochastic_order(self):
        source = self.source()
        archive = self.root / "shared.h5"
        encode(source, archive, self.config(max_component_mse=0, max_complex_error=0,
                                          ordering="amplitude", shared_order=True, rounding="stochastic-pcg64-v1"))
        self.assertEqual(verify(archive, source)["status"], "passed")
        with Reader(archive) as reader:
            spec = reader.manifest["tables"][0]["columns"]["DATA"]
            decoded = read_visibility(reader, spec)
            with table(str(source), readonly=True, lockoptions="usernoread", ack=False) as t:
                self.assertEqual(decoded.tobytes(), t.getcol("DATA").tobytes())
            self.assertEqual(spec["products"][0]["permutation"], spec["products"][1]["permutation"])

    def test_unsupported_quantized_layout_fails_before_creating_output(self):
        source = make_ms(self.root / "source.ms", pairs=((0, 1), (0, 1)), channels=1, ddids=[0, 1])
        output = self.root / "unsupported.h5"
        with self.assertRaisesRegex(ValueError, "regular configuration"):
            encode(source, output, self.config())
        self.assertFalse(output.exists())
        self.assertEqual(list(self.root.glob("*.partial")), [])

    def test_block_corruption_is_detected_without_decoding_prior_blocks(self):
        source = self.source()
        archive = self.root / "corrupt.h5"
        encode(source, archive, self.config())
        with Reader(archive) as reader:
            spec = reader.manifest["tables"][0]["columns"]["DATA"]
            path = spec["products"][0]["payload"]
        with h5py.File(archive, "r+") as f:
            f[path][20] ^= 1
        with Reader(archive, verify=False) as reader:
            with self.assertRaisesRegex(ValueError, "CRC32"):
                read_visibility_block(reader, spec, 0, 0)

    def test_stokes_multiple_selected_and_exact_optional_columns(self):
        source = make_ms(self.root / "stokes.ms", channels=17, correlations=(1, 2, 3, 4))
        rng = np.random.default_rng(80)
        values = (rng.normal(size=(3, 17, 4)) + 1j * rng.normal(size=(3, 17, 4))).astype(np.complex64)
        with table(str(source), readonly=False, ack=False) as t:
            t.addcols(maketabdesc([makearrcoldesc("CORRECTED_DATA", 0j, ndim=2, valuetype="complex"),
                                  makearrcoldesc("FLOAT_DATA", 0., ndim=2, valuetype="float"),
                                  makearrcoldesc("WEIGHT_SPECTRUM", 0., ndim=2, valuetype="float"),
                                  makescacoldesc("CUSTOM_FLOAT", 0., valuetype="float")]))
            t.putcol("DATA", values)
            t.putcol("CORRECTED_DATA", values * np.float32(.001))
            floats = np.zeros(values.shape, dtype=np.float32)
            floats.view(np.uint32).flat[:3] = [0x7FC12345, 0x7F800000, 0x80000000]
            t.putcol("FLOAT_DATA", floats)
            t.putcol("WEIGHT_SPECTRUM", floats)
            t.putcol("CUSTOM_FLOAT", np.asarray([0x7FC12345, 0x80000000, 0x7F800000], dtype=np.uint32).view(np.float32))
        archive, restored = self.root / "stokes.h5", self.root / "restored.ms"
        encode(source, archive, self.config(columns=["DATA", "CORRECTED_DATA"],
                                           product_bounds={"Q": {"max_component_mse": 0, "max_complex_error": 0}}))
        decode(archive, restored)
        self.assertEqual(verify(archive, source, restored)["status"], "passed")
        with table(str(source), readonly=True, lockoptions="usernoread", ack=False) as a, table(str(restored), readonly=True, lockoptions="usernoread", ack=False) as b:
            for name in ("FLOAT_DATA", "WEIGHT_SPECTRUM", "CUSTOM_FLOAT"):
                self.assertEqual(a.getcol(name).tobytes(), b.getcol(name).tobytes())


if __name__ == "__main__":
    unittest.main()
