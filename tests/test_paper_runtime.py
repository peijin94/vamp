import copy
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location("summarize_runtime", Path(__file__).resolve().parents[1] / "scripts/summarize_runtime.py")
runtime = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runtime)


class RuntimeTests(unittest.TestCase):
    def encoding(self):
        return {"timings_seconds": {"inspection_preflight": 2, "before_integrity_finalization": 19},
                "elapsed_seconds": 20, "visibility_block_statistics": {"DATA": {"timings_seconds": {
                    "codec_total": 10, "packing_and_secondary_compression": 4,
                    "quantization_and_scale_estimation": 3, "candidate_reconstruction_and_error_checks": 2}}},
                "archive_bytes": 100, "source_physical_file_bytes": 400,
                "selected_visibility_logical_bytes": 320, "peak_process_rss_bytes": 5000}

    def test_nested_timers_not_double_counted(self):
        values = runtime.stages(self.encoding())
        self.assertEqual(sum(values.values()), 20)
        self.assertNotIn("codec_total", values)
        self.assertEqual(values["other_exact_columns_and_uninstrumented_overhead"], 8)
        invalid = self.encoding()
        invalid["elapsed_seconds"] = 1
        invalid["timings_seconds"]["before_integrity_finalization"] = 1
        with self.assertRaises(ValueError):
            runtime.stages(invalid)

    def test_repeats_require_verified_identical_config_and_process_time(self):
        record = {"name": "repeat", "configuration": {"mode": "adaptive", "bits": [4, 8, 12, 16]},
                  "verification_status": "passed", "encoding": self.encoding(),
                  "process_timings_seconds": {"encode_process": 21, "verification_process": 3}}
        benchmark = {"status": "complete", "source": "input.ms", "cache_methodology": "warm", "notes": "test",
                     "experiments": [record, copy.deepcopy(record)]}
        result = runtime.summarize(benchmark)
        self.assertEqual(result["adaptive_summary"]["encode_process_seconds"]["median"], 21)
        benchmark["experiments"][1]["configuration"]["bits"] = [4]
        with self.assertRaises(ValueError):
            runtime.summarize(benchmark)
        benchmark["experiments"][1]["verification_status"] = "failed"
        with self.assertRaises(ValueError):
            runtime.summarize(benchmark)


if __name__ == "__main__":
    unittest.main()
