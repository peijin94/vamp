from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from vamp import codec, visibility
from vamp.archive import encode, verify
from vamp.config import Config, SUPPORTED_BITS
from vamp.container import Reader
from test_inspection import make_ms


@unittest.skipIf(codec._native is None, "optional native extension unavailable")
class AccelerationTests(unittest.TestCase):
    def test_native_buffer_bounds_and_worker_validation(self):
        native = codec._native
        a = np.ones((2,3), dtype=np.complex64)
        flags = np.zeros(a.shape, dtype=bool)
        for workers in (0, -1, 65, True, 1.5):
            with self.assertRaises(ValueError):
                Config(workers=workers).validate()
        for bad in (b"", b"\0", b"\0"*7):
            with self.assertRaises(ValueError):
                native.candidate(a, bad, 8, 1., 3)
            with self.assertRaises(ValueError):
                native.quality(a, a, bad, 2, 3, 1)
        with self.assertRaises(ValueError):
            native.search(a, flags, b"\x09", np.asarray([1.],dtype="f4"), np.ones(4), 3, 1)
        with self.assertRaises(ValueError):
            native.search(a, flags, b"\x08", b"", np.ones(4), 3, 1)
        for channels in (0, 2, -1):
            with self.assertRaises(ValueError):
                native.quality(a, a, flags, 2, channels, 1)
        for flagged in (False, True):
            value = np.asarray([[[1+2j]]],dtype="c8")
            summary, arrays = visibility.quality(value,value,np.full(value.shape,flagged),["I"])
            self.assertEqual(summary["by_product"]["I"]["component_mse"],0)
            self.assertEqual(arrays["adjacent_channel_error_covariance"].shape,(0,1))
        wide = np.zeros((1,65537,1), dtype="c8")
        summary, _ = visibility.quality(wide,wide,np.zeros(wide.shape,bool),["I"])
        self.assertNotIn("calculation", summary)

    def test_native_and_numpy_choose_identical_serialized_candidates(self):
        rng = np.random.default_rng(119)
        for factor in (0., 1e-30, .001, 1., 1e25):
            a = ((rng.normal(size=(9, 31)) + 1j*rng.normal(size=(9, 31))) * factor).astype(np.complex64)
            a.view(np.uint32).flat[:6] = [0x80000000, 0, 0x7fc12345, 0xff800000, 0x7f800000, 0]
            flags = rng.random(a.shape) < .13
            for bits in (*[(b,) for b in SUPPORTED_BITS], (3,4,6,8,12,16), (4,8,12,16)):
                config = Config(mode="adaptive", bits=list(bits), scale_quantiles=[1., .95],
                                max_component_mse=max(1e-70, factor**2*.001),
                                max_complex_error=max(1e-35, factor*.1))
                native = codec.encode_block(a, flags, config, config.bounds())
                with patch.object(codec, "_native", None):
                    numpy = codec.encode_block(a, flags, config, config.bounds())
                self.assertEqual(native.blob, numpy.blob)
                self.assertEqual(native.decoded.tobytes(), numpy.decoded.tobytes())
                self.assertEqual(native.metrics, numpy.metrics)

    def test_near_boundaries_and_p99_are_checked_strictly(self):
        rng = np.random.default_rng(17)
        a = (rng.normal(size=(8,32)) + 1j*rng.normal(size=(8,32))).astype(np.complex64)
        flags = np.zeros(a.shape, bool)
        config = Config(mode="fixed", bits=[8], scale_quantiles=[1.], max_component_mse=1., max_complex_error=1.)
        reference = codec.encode_block(a, flags, config, config.bounds())
        for bound, metric in (("max_component_mse","component_mse"), ("max_complex_bias","complex_bias_abs"),
                              ("max_complex_error","max_complex_error"), ("max_baseline_component_mse","max_baseline_component_mse"),
                              ("max_p99_error","p99_error")):
            for direction in (0., np.inf):
                candidate = replace(config, **{bound: float(np.nextafter(reference.metrics[metric], direction))})
                native = codec.encode_block(a, flags, candidate, candidate.bounds())
                with patch.object(codec, "_native", None):
                    numpy = codec.encode_block(a, flags, candidate, candidate.bounds())
                self.assertEqual(native.blob, numpy.blob)
                self.assertTrue(codec.within_bounds(native.metrics, candidate.bounds()))

    def test_quality_moments_match_reference_and_workers_are_deterministic(self):
        rng = np.random.default_rng(3)
        a = (rng.normal(size=(2100, 9, 4)) + 1j*rng.normal(size=(2100, 9, 4))).astype(np.complex64)
        b = (a + (rng.normal(size=a.shape) + 1j*rng.normal(size=a.shape))*.001).astype(np.complex64)
        flags = rng.random(a.shape)<.15
        a.view(np.uint32).flat[:4] = [0x80000000,0,0x7fc12345,0xff800000]
        protected = visibility.protected_mask(a, flags); b[protected] = a[protected]
        names = ["XX","XY","YX","YY"]
        single, arrays = visibility.quality(a,b,flags,names,workers=1)
        parallel, parallel_arrays = visibility.quality(a,b,flags,names,workers=4)
        self.assertEqual(single, parallel)
        for key in arrays:
            self.assertEqual(arrays[key].tobytes(), parallel_arrays[key].tobytes())
        reference, reference_arrays = visibility.quality(a,b,flags,names,method="numpy-v1")
        self.assertEqual(single["protected_sample_count"], reference["protected_sample_count"])
        for name in names:
            for key, value in reference["by_product"][name].items():
                np.testing.assert_allclose(single["by_product"][name][key],value,rtol=1e-12,atol=1e-20)
        for key in arrays:
            np.testing.assert_allclose(arrays[key], reference_arrays[key], rtol=1e-12, atol=1e-20)
        changed = b.copy(); changed[protected] = 0
        with self.assertRaises(ValueError):
            visibility.quality(a, changed, flags, names)

    def test_parallel_archives_and_portable_verification(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = make_ms(root/"input.ms", channels=17)
            config = Config(mode="adaptive", baseline_group=2, channel_segment=8,
                            max_component_mse=1., max_complex_error=1.)
            outputs = []
            for workers in (1, 4):
                path = root/f"w{workers}.h5"
                encode(source, path, replace(config,workers=workers))
                self.assertEqual(verify(path, source)["status"], "passed")
                with Reader(path) as reader:
                    outputs.append({name: entry["sha256"] for name,entry in reader.entries.items()})
                with patch.object(visibility,"_native",None), patch.object(codec,"_native",None):
                    self.assertEqual(verify(path, source)["status"], "passed")
            self.assertEqual(*outputs)


if __name__ == "__main__":
    unittest.main()
