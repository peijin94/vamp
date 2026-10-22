import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

from varquanms.config import Config

_spec = importlib.util.spec_from_file_location("summarize_grouping", Path(__file__).resolve().parents[1] / "scripts" / "summarize_grouping.py")
summary = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(summary)


class SummaryTests(unittest.TestCase):
    def test_rejects_mixed_imaging_controls_before_aggregating(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "case").mkdir()
            matrix = {"schema": "varquanms.benchmark/1", "base": {}, "export_ms": True,
                      "cache_methodology": "fixture", "experiments": [{"name": "case", "config": {}}]}
            (root / "matrix.json").write_text(json.dumps(matrix))
            (root / "case" / "complete.json").write_text(json.dumps({"configuration": Config().to_dict(), "verification_status": "passed"}))
            repeat = {"parameters": {"threads": 1}, "container_id": "same", "reference": "reference"}
            (root / "repeat.json").write_text(json.dumps(repeat))
            image = {**repeat, "status": "measured", "parameters": {"threads": 8}}
            (root / "wsclean-case.json").write_text(json.dumps(image))
            with self.assertRaisesRegex(ValueError, "Imaging controls differ"):
                summary.summarize(root / "matrix.json", [root], root, root / "repeat.json")


if __name__ == "__main__":
    unittest.main()
