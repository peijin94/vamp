from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
import io
import json
from pathlib import Path
import tempfile
import unittest

from casacore.tables import table
import numpy as np

from varquanms.archive import compare_ms, decode, encode, verify
from varquanms.cli import main
from varquanms.config import Config
from varquanms.container import Reader
from varquanms.flags import ENCODING, FLAG_COLUMNS, POLICY_KEYWORD, burned_columns
from test_inspection import make_ms, snapshot


class BurnFlagTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def source(self, irregular=False):
        source = make_ms(self.root / "source.ms", shapes=[(3, 2), (1, 2), (3, 2)] if irregular else None)
        with table(str(source), readonly=False, ack=False) as t:
            for row in range(t.nrows()):
                shape = t.getcell("FLAG", row).shape
                t.putcell("FLAG", row, np.full(shape, row != 1, dtype=bool))
                t.putcell("FLAG_ROW", row, row == 1)
                # Keep one undefined FLAG_CATEGORY cell; an empty defined cell differs.
                if row != 1:
                    t.putcell("FLAG_CATEGORY", row, np.ones((0 if row == 2 else 1, *shape), dtype=bool))
            values = t.getcell("DATA", 0)
            values.view(np.uint32).flat[:4] = [0x7FC12345, 0x80000000, 0x7F800000, 0xFF800000]
            t.putcell("DATA", 0, values)
        with table(str(source / "ANTENNA"), readonly=False, ack=False) as t:
            t.putcell("FLAG_ROW", 0, True)
        return source

    def test_default_preserves_flags_and_burn_is_strict_boolean(self):
        self.assertFalse(Config().burn_flag)
        for invalid in (0, 1, "true", "false", None, []):
            with self.assertRaisesRegex(ValueError, "burn_flag must be boolean"):
                Config(burn_flag=invalid).validate()
        source = self.source()
        archive, output = self.root / "exact.h5", self.root / "exact.ms"
        summary = encode(source, archive)
        decode(archive, output)
        self.assertTrue(summary["exact_ms_round_trip"])
        self.assertTrue(compare_ms(source, output))
        with table(str(output), readonly=True, lockoptions="usernoread", ack=False) as t:
            self.assertNotIn(POLICY_KEYWORD, t.keywordnames())

    def test_drop_values_keep_shapes_presence_weights_metadata_and_source(self):
        for irregular in (False, True):
            with self.subTest(irregular=irregular), tempfile.TemporaryDirectory() as directory:
                self.root = Path(directory)
                source = self.source(irregular)
                original = snapshot(source)
                archive, output = self.root / "burn.h5", self.root / "burn.ms"
                summary = encode(source, archive, Config(burn_flag=True, io_chunk_rows=2))
                self.assertFalse(summary["exact_ms_round_trip"])
                self.assertTrue(summary["flag_policy"]["burn_flag"])
                self.assertEqual(summary["archive_bytes"], sum(summary["stored_bytes_by_category"].values()))
                with Reader(archive) as reader:
                    columns = reader.manifest["tables"][0]["columns"]
                    for name in FLAG_COLUMNS:
                        spec = columns[name]
                        self.assertEqual(spec["encoding"], ENCODING)
                        self.assertNotIn("values", spec)
                        for chunk in spec.get("shape_chunks", []):
                            self.assertTrue(all(shape is None or isinstance(shape, list) for shape in reader.json(chunk["path"])))
                    self.assertEqual(reader.manifest["provenance"]["flag_policy"], summary["flag_policy"])
                self.assertEqual(decode(archive, output)["flag_policy"], summary["flag_policy"])
                result = verify(archive, source, output)
                self.assertEqual(result["status"], "passed")
                discarded = [check for check in result["column_checks"] if check.get("status") == "intentionally_discarded_shapes_preserved"]
                self.assertEqual({check["column"] for check in discarded}, set(FLAG_COLUMNS))
                with table(str(source), readonly=True, lockoptions="usernoread", ack=False) as a, table(str(output), readonly=True, lockoptions="usernoread", ack=False) as b:
                    self.assertEqual(b.getkeyword(POLICY_KEYWORD)["policy"], ENCODING)
                    for name in FLAG_COLUMNS:
                        for row in range(a.nrows()):
                            self.assertEqual(a.iscelldefined(name, row), b.iscelldefined(name, row))
                            if a.iscelldefined(name, row):
                                value = np.asarray(b.getcell(name, row))
                                self.assertEqual(value.shape, np.asarray(a.getcell(name, row)).shape)
                                self.assertFalse(value.any())
                with table(str(output / "ANTENNA"), readonly=True, lockoptions="usernoread", ack=False) as t:
                    self.assertTrue(t.getcell("FLAG_ROW", 0))
                self.assertEqual(snapshot(source), original)
                with table(str(output), readonly=False, ack=False) as t:
                    t.putcell("FLAG_ROW", 0, True)
                with self.assertRaisesRegex(ValueError, "flag reconstruction differs: ./FLAG_ROW"):
                    verify(archive, source, output)

    def test_existing_provenance_is_not_silently_overwritten(self):
        source = self.source()
        archive, output = self.root / "first.h5", self.root / "first.ms"
        encode(source, archive, Config(burn_flag=True))
        decode(archive, output)
        # An identical marker survives another archive generation.
        encode(output, self.root / "second.h5", Config(burn_flag=True))
        with table(str(source), readonly=False, ack=False) as t:
            t.putkeyword(POLICY_KEYWORD, {"user_owned": True})
        with self.assertRaisesRegex(ValueError, "Cannot overwrite"):
            encode(source, self.root / "conflict.h5", Config(burn_flag=True))
        self.assertFalse((self.root / "conflict.h5").exists())
        self.assertEqual(list(self.root.glob("*.partial")), [])

    def test_regular_flag_arrays_have_no_payload(self):
        source = make_ms(self.root / "source.ms")
        with table(str(source), readonly=False, ack=False) as t:
            t.putcol("FLAG_CATEGORY", np.ones((3, 1, 3, 2), dtype=bool))
        archive = self.root / "burn.h5"
        encode(source, archive, Config(burn_flag=True, deduplicate=False))
        with Reader(archive) as reader:
            main_table = reader.manifest["tables"][0]
            for name in FLAG_COLUMNS:
                spec = main_table["columns"][name]
                self.assertEqual(set(spec), {"encoding", "row_count", "dtype", "cell_shape"})
                index = main_table["column_names"].index(name)
                self.assertFalse(any(path.startswith(f"tables/t0/c{index}/") for path in reader.entries))
            changed = deepcopy(reader.manifest)
            changed["config"]["burn_flag"] = False
            with self.assertRaisesRegex(ValueError, "policy"):
                burned_columns(changed)
            changed = deepcopy(reader.manifest)
            changed["tables"][0]["columns"]["UVW"]["encoding"] = ENCODING
            with self.assertRaisesRegex(ValueError, "policy"):
                burned_columns(changed)

    def test_cli_and_config_opt_in_warn_and_default_does_not_override_config(self):
        source = self.source()
        config_path = self.root / "config.json"
        config_path.write_text(json.dumps({"burn_flag": True}))
        for index, args in enumerate((["--burn-flag"], ["--burn_flag"], ["--config", str(config_path)])):
            out, err = io.StringIO(), io.StringIO()
            with redirect_stdout(out), redirect_stderr(err):
                result = main(["encode", str(source), str(self.root / f"cli{index}.h5"), *args])
            self.assertEqual(result, 0, err.getvalue())
            self.assertIn("WARNING", err.getvalue())
            self.assertTrue(json.loads(out.getvalue())["flag_policy"]["burn_flag"])


if __name__ == "__main__":
    unittest.main()
