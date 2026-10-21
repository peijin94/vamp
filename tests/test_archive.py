import json
from pathlib import Path
import tempfile
import unittest

from casacore.tables import table
import h5py
import numpy as np

from varquanms.archive import compare_ms, decode, encode, verify
from varquanms.config import Config
from varquanms.container import Reader
from varquanms.metadata import equal, pack, unpack
from test_inspection import make_ms, snapshot


class LosslessArchiveTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def round_trip(self, path, deduplicate=True):
        original = snapshot(path)
        archive = self.root / "archive.h5"
        output = self.root / "restored.ms"
        summary = encode(path, archive, Config(deduplicate=deduplicate))
        self.assertEqual(summary["archive_bytes"], sum(summary["stored_bytes_by_category"].values()))
        self.assertEqual(verify(archive, source=path)["status"], "passed")
        self.assertEqual(decode(archive, output)["python_casacore_reopen"], "passed")
        self.assertTrue(compare_ms(path, output))
        self.assertEqual(verify(archive, source=path, restored=output)["status"], "passed")
        self.assertEqual(original, snapshot(path))
        return archive, output

    def test_lossless_regular_exact_metadata_and_optional_column(self):
        path = make_ms(self.root / "source.ms")
        with table(str(path), readonly=False, ack=False) as t:
            t.putkeyword("nested", {"example": np.asarray([-0.0, 2.0]), "label": "μJy"})
            values = np.zeros((3, 2), dtype=np.complex64)
            values.view(np.uint32).flat[:4] = [0x7FC12345, 0x80000000, 0x7F800000, 0xFF800000]
            t.putcell("DATA", 1, values)
            t.putcell("FLAG_ROW", 1, True)
        self.round_trip(path)

    def test_lossless_without_structural_dedup(self):
        self.round_trip(make_ms(self.root / "source.ms"), deduplicate=False)

    def test_export_writes_more_than_one_io_chunk(self):
        source = make_ms(self.root / "source.ms", pairs=((0, 1),) * 1025, channels=1)
        self.round_trip(source)

    def test_lossless_explicit_irregular_cells(self):
        path = make_ms(self.root / "source.ms", shapes=[(3, 2), (1, 2), (3, 2)])
        self.round_trip(path)

    def test_dictionary_of_multidimensional_flag_cells(self):
        path = make_ms(self.root / "source.ms", pairs=((0, 1),) * 12, channels=17)
        with table(str(path), readonly=False, ack=False) as t:
            flags = np.zeros((12, 17, 2), dtype=bool)
            flags[2:6] = True
            t.putcol("FLAG", flags)
        self.round_trip(path)

    def test_lossless_single_channel_duplicate_and_mixed_contexts(self):
        path = make_ms(self.root / "source.ms", pairs=((0, 1), (0, 1), (0, 1), (1, 0)), channels=1, ddids=[0, 0, 1, 1], fields=[0, 0, 1, 1])
        self.round_trip(path)

    def test_corruption_incomplete_and_truncation_rejected(self):
        path = make_ms(self.root / "source.ms")
        archive, _ = self.round_trip(path)
        with h5py.File(archive, "r+") as f:
            f["tables/t0/c0/values"][0] += 1
        with self.assertRaisesRegex(ValueError, "checksum"):
            verify(archive)
        with h5py.File(archive, "r+") as f:
            f.attrs["complete"] = False
        with self.assertRaisesRegex(ValueError, "incomplete"):
            Reader(archive)
        raw = archive.read_bytes()
        truncated = self.root / "truncated.h5"
        truncated.write_bytes(raw[:len(raw) // 2])
        with self.assertRaises((ValueError, OSError, RuntimeError)):
            Reader(truncated)

    def test_existing_source_or_output_is_never_overwritten(self):
        source = make_ms(self.root / "source.ms")
        archive, output = self.round_trip(source)
        for destination in (source, archive, source / "new.h5"):
            with self.assertRaises(ValueError):
                encode(source, destination)
        with self.assertRaises(ValueError):
            decode(archive, output)
        with self.assertRaises(ValueError):
            Config.load(None, invented=True)

    def test_numeric_metadata_preserves_bits_and_scalar_kind(self):
        values = {"float": -0.0, "array": np.asarray([1, 0x7FC12345], dtype=np.uint32).view(np.float32),
                  "complex": complex(2, -0.0), "numpy_scalar": np.float32(-0.0),
                  "strings": ["µ", "", "abc"], "int": 2**64}
        restored = unpack(json.loads(json.dumps(pack(values))))
        self.assertTrue(equal(values, restored))
        self.assertIs(type(restored["float"]), float)
        self.assertIs(type(restored["numpy_scalar"]), np.float32)


if __name__ == "__main__":
    unittest.main()
