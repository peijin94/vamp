"""Small MS fixtures exercise metadata, row identity and special-value handling."""

import base64
from contextlib import redirect_stderr, redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest

from casacore.tables import default_ms, makearrcoldesc, maketabdesc, table
import numpy as np

from varquanms.cli import main
from varquanms.inspection import inspect_ms, json_value, numeric_counts, unique_bits


def make_ms(path, pairs=((0, 0), (0, 1), (1, 1)), channels=3,
            ddids=None, fields=None, shapes=None, correlations=(12, 9)):
    """Generate an actual casacore MS, including required referenced subtables."""
    nrow = len(pairs)
    ddids = [0] * nrow if ddids is None else ddids
    fields = [0] * nrow if fields is None else fields
    ncorr = len(correlations)
    desc = maketabdesc(makearrcoldesc("DATA", 0j, ndim=2, valuetype="complex"))
    with default_ms(str(path), desc) as t:
        t.addrows(nrow)
        for row, pair in enumerate(pairs):
            for col, value in {
                "ANTENNA1": pair[0], "ANTENNA2": pair[1],
                "FEED1": 0, "FEED2": 0, "ARRAY_ID": 0,
                "DATA_DESC_ID": ddids[row], "FIELD_ID": fields[row],
                "SCAN_NUMBER": 1, "OBSERVATION_ID": 0,
                "PROCESSOR_ID": -1, "STATE_ID": -1,
                "TIME": 5295265208.539846, "TIME_CENTROID": 5295265213.555357,
                "INTERVAL": 10.0, "EXPOSURE": 10.0, "FLAG_ROW": False,
            }.items():
                t.putcell(col, row, value)
            shape = (channels, ncorr) if shapes is None else shapes[row]
            t.putcell("DATA", row, np.zeros(shape, dtype=np.complex64))
            t.putcell("FLAG", row, np.zeros(shape, dtype=bool))
            t.putcell("WEIGHT", row, np.ones(ncorr, dtype=np.float32))
            t.putcell("SIGMA", row, np.ones(ncorr, dtype=np.float32))
            t.putcell("UVW", row, np.asarray([0.0, -0.0, float(pair[1] - pair[0])]))
    for name, count in (("ANTENNA", max(max(pair) for pair in pairs) + 1), ("FIELD", max(fields) + 1),
                        ("DATA_DESCRIPTION", max(ddids) + 1),
                        ("SPECTRAL_WINDOW", max(ddids) + 1), ("POLARIZATION", 1)):
        with table(str(path / name), readonly=False, ack=False) as t:
            t.addrows(count)
            if name == "POLARIZATION":
                t.putcell("NUM_CORR", 0, ncorr)
                t.putcell("CORR_TYPE", 0, np.asarray(correlations, dtype=np.int32))
                t.putcell("CORR_PRODUCT", 0, np.zeros((ncorr, 2), dtype=np.int32))
            elif name == "DATA_DESCRIPTION":
                for row in range(count):
                    t.putcell("SPECTRAL_WINDOW_ID", row, row)
                    t.putcell("POLARIZATION_ID", row, 0)
            elif name == "SPECTRAL_WINDOW":
                for row in range(count):
                    t.putcell("NUM_CHAN", row, channels)
                    t.putcell("CHAN_FREQ", row, np.arange(channels, dtype=np.float64) * 24_000 + 59e6 + row * 24_000)
                    t.putcell("CHAN_WIDTH", row, np.full(channels, 24_000.0))
                    t.putcell("MEAS_FREQ_REF", row, 5)
    return path


def snapshot(path):
    return {str(p.relative_to(path)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in path.rglob("*") if p.is_file()}


class InspectionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_regular_ms_product_order_and_no_writes(self):
        path = make_ms(self.root / "input.ms")
        before = snapshot(path)
        report = inspect_ms(path, chunk_rows=2)
        self.assertEqual(before, snapshot(path))
        self.assertEqual(report["status"], "complete")
        self.assertEqual(report["main"]["row_count"], 3)
        self.assertEqual(report["layout"]["organization"], "channels_inside_cells")
        self.assertTrue(report["layout"]["regular_single_configuration_candidate"])
        config = report["layout"]["configurations"][0]
        self.assertEqual(config["correlation_names"], ["YY", "XX"])
        self.assertEqual(config["missing_cross_pairs_relative_to_full_antenna_table"], 0)
        self.assertEqual(report["main"]["columns"]["DATA"]["logical_array_bytes"], 3 * 3 * 2 * 8)
        self.assertEqual(report["main"]["columns"]["FLAG_CATEGORY"]["undefined_cells"], 3)
        self.assertEqual(report["uvw"]["bitwise_distinct_geometry_records"], 2)
        self.assertTrue(report["layout"]["archive_implementation_available"])
        json.dumps(report, allow_nan=False)

    def test_single_channel_records_keep_ddids_distinct(self):
        path = make_ms(self.root / "channels.ms", pairs=((0, 1), (0, 1)), channels=1, ddids=[0, 1])
        report = inspect_ms(path, chunk_rows=1)
        self.assertEqual(report["layout"]["organization"], "single_channel_data_description_records")
        self.assertEqual(report["row_identity"]["row_key_duplicate_count"], 0)
        self.assertEqual(len(report["layout"]["configurations"]), 2)
        self.assertFalse(report["layout"]["regular_single_configuration_candidate"])

    def test_duplicates_mixed_fields_and_ordered_antennas(self):
        path = make_ms(self.root / "identities.ms", pairs=((0, 1), (0, 1), (1, 0), (0, 1)), fields=[0, 0, 0, 1])
        report = inspect_ms(path)
        identities = report["row_identity"]
        self.assertEqual(identities["row_key_duplicate_count"], 1)
        self.assertEqual(identities["exact_duplicate_rows_beyond_first"], 1)
        self.assertEqual(identities["duplicate_row_groups"], [[0, 1]])
        self.assertEqual(identities["antenna1_greater_than_antenna2_rows"], 1)
        self.assertEqual(identities["context_count_excluding_pair_and_data_description"], 2)
        self.assertFalse(report["layout"]["regular_single_configuration_candidate"])

    def test_heterogeneous_shapes_are_reported(self):
        path = make_ms(self.root / "shapes.ms", shapes=[(3, 2), (1, 2), (3, 2)])
        report = inspect_ms(path)
        self.assertEqual(report["main"]["columns"]["DATA"]["cell_shapes"], {"[3, 2]": 2, "[1, 2]": 1})
        self.assertFalse(report["layout"]["regular_single_configuration_candidate"])

    def test_flags_nonfinite_and_signed_zero_bits(self):
        path = make_ms(self.root / "special.ms", pairs=((0, 1),), channels=2)
        values = np.zeros((2, 2), dtype=np.complex64)
        values.view(np.uint32).flat[0] = 0x7FC12345
        values.view(np.uint32).flat[1] = 0x80000000
        values.view(np.uint32).flat[2] = 0x7F800000
        with table(str(path), readonly=False, ack=False) as t:
            t.putcell("DATA", 0, values)
            t.putcell("FLAG_ROW", 0, True)
        report = inspect_ms(path)
        stats = report["main"]["columns"]["DATA"]["statistics"]
        self.assertEqual(stats["nonfinite_values"], 2)
        self.assertEqual(stats["negative_zero_components"], 1)
        vis = report["visibility"]["DATA"]
        self.assertEqual(vis["flagged_samples_including_flag_row"], 4)
        self.assertEqual(vis["all_flagged_rows"], 1)
        self.assertEqual(vis["nonfinite_flagged_samples"], 2)
        encoded = json_value(values)
        restored = np.frombuffer(base64.b64decode(encoded["data_base64"]), dtype=encoded["dtype"]).reshape(encoded["shape"])
        self.assertEqual(restored.tobytes(), values.tobytes())

    def test_bitwise_geometry_keeps_signed_zeros_distinct(self):
        values = np.asarray([[0.0, 0.0, 0.0], [-0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        unique, mapping, counts = unique_bits(values)
        self.assertEqual(len(unique), 2)
        self.assertEqual(mapping.tolist(), [0, 1, 0])
        self.assertEqual(counts.tolist(), [2, 1])
        self.assertEqual(unique[mapping].tobytes(), values.tobytes())

    def test_output_never_overwrites_or_writes_into_source(self):
        path = make_ms(self.root / "safe.ms")
        output = self.root / "inspection.json"
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            self.assertEqual(main(["inspect", str(path), "--report", str(output)]), 0)
            original = output.read_bytes()
            self.assertEqual(main(["inspect", str(path), "--report", str(output)]), 1)
            self.assertEqual(main(["inspect", str(path), "--report", str(path / "new.json")]), 1)
        self.assertEqual(original, output.read_bytes())
        self.assertFalse((path / "new.json").exists())

    def test_missing_selected_column_fails_before_final_report(self):
        path = make_ms(self.root / "missing.ms")
        output = self.root / "missing.json"
        with redirect_stderr(io.StringIO()):
            self.assertEqual(main(["inspect", str(path), "--report", str(output), "--columns", "CORRECTED_DATA"]), 1)
        self.assertFalse(output.exists())

    def test_unknown_options_are_rejected(self):
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            main(["inspect", "input.ms", "--report", "r.json", "--imaginary-option"])
        self.assertEqual(error.exception.code, 2)

    def test_stokes_names_are_preserved(self):
        path = make_ms(self.root / "stokes.ms", correlations=(1, 2, 3, 4))
        report = inspect_ms(path)
        self.assertEqual(report["layout"]["configurations"][0]["correlation_names"], ["I", "Q", "U", "V"])


if __name__ == "__main__":
    unittest.main()
