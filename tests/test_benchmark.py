import json
from pathlib import Path
import tempfile
import unittest

from vamp.benchmark import load_matrix, run
from test_inspection import make_ms


class BenchmarkTests(unittest.TestCase):
    def test_two_isolated_end_to_end_runs_and_invalid_matrix(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = make_ms(root / "source.ms", channels=17)
            matrix = {"schema": "vamp.benchmark/1", "base": {}, "export_ms": True,
                      "cache_methodology": "synthetic fixture; no cache eviction",
                      "experiments": [{"name": "lossless", "config": {}},
                                      {"name": "adaptive", "config": {"mode": "adaptive", "baseline_group": 3,
                                                                          "channel_segment": 16, "max_component_mse": .01,
                                                                          "max_complex_error": 1.0}}]}
            matrix_path = root / "matrix.json"
            matrix_path.write_text(json.dumps(matrix))
            legacy = dict(matrix, schema="varquanms.benchmark/1")
            matrix_path.write_text(json.dumps(legacy))
            self.assertEqual(len(load_matrix(matrix_path)[1]), 2)
            matrix_path.write_text(json.dumps(matrix))
            result = run(source, matrix_path, root / "results")
            self.assertEqual(result["status"], "complete")
            self.assertEqual([e["verification_status"] for e in result["experiments"]], ["passed", "passed"])
            matrix["experiments"][1]["name"] = "../unsafe"
            matrix_path.write_text(json.dumps(matrix))
            with self.assertRaises(ValueError):
                run(source, matrix_path, root / "unsafe-results")
            self.assertFalse((root / "unsafe-results").exists())


if __name__ == "__main__":
    unittest.main()
