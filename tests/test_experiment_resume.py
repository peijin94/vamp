"""Concurrent experiment stages must never expose partial configuration JSON."""
from concurrent.futures import ThreadPoolExecutor
import importlib.util
import json
from pathlib import Path
import tempfile
import threading
import unittest

spec=importlib.util.spec_from_file_location('run_error_limits',
    Path(__file__).resolve().parents[1]/'scripts/run_error_limits.py')
runner=importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class ExperimentResumeTests(unittest.TestCase):
    def test_concurrent_initialization_and_unchanged_resume(self):
        config={'epsilon_D':.18,'bits':[3,4,5,6,8,10,12,14,16],
                'metadata':'x'*100000}
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'case.json'
            barrier=threading.Barrier(12)
            def run(_):
                barrier.wait()
                for _ in range(10):
                    runner.ensure_case_config(path,config)
                    self.assertEqual(json.loads(path.read_text()),config)
            with ThreadPoolExecutor(max_workers=12) as pool:
                list(pool.map(run,range(12)))
            before=path.stat()
            runner.ensure_case_config(path,config)
            after=path.stat()
            self.assertEqual((before.st_ino,before.st_mtime_ns),
                             (after.st_ino,after.st_mtime_ns))
            self.assertEqual(list(Path(d).iterdir()),[path])

    def test_conflicting_concurrent_config_does_not_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'case.json'
            barrier=threading.Barrier(2)
            def run(value):
                barrier.wait()
                try:
                    runner.ensure_case_config(path,{'epsilon_D':value})
                    return value,True
                except ValueError:
                    return value,False
            with ThreadPoolExecutor(max_workers=2) as pool:
                results=list(pool.map(run,[.018,.18]))
            winner=[value for value,success in results if success]
            self.assertEqual(len(winner),1)
            self.assertEqual(json.loads(path.read_text()),{'epsilon_D':winner[0]})
            self.assertEqual(list(Path(d).iterdir()),[path])


if __name__=='__main__':unittest.main()
