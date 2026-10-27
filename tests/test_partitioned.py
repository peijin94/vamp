import json
from pathlib import Path
import tempfile
import unittest

from casacore.tables import table, makearrcoldesc, maketabdesc, tabledefinehypercolumn
import numpy as np

from vamp.archive import compare_ms
from vamp.config import Config
from vamp.partitioned import context_runs, decode_partitioned, encode_partitioned
from vamp.imaging import file_sha256
from test_inspection import make_ms, snapshot


class PartitionedTests(unittest.TestCase):
    def test_independent_lossless_reconstruction_and_corrupt_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = make_ms(root / 'source.ms', pairs=((0, 1), (1, 0), (0, 1), (1, 0)))
            with table(str(source), readonly=False, ack=False) as t:
                desc = maketabdesc(makearrcoldesc('CORRECTED_DATA', 0j, ndim=2, valuetype='complex'))
                tabledefinehypercolumn(desc, 'TiledCorrected', 3, ['CORRECTED_DATA'])
                t.addcols(desc)
                t.putcol('TIME', np.array([1., 1., 2., 2.]))
                t.putcol('DATA', np.arange(24, dtype=np.float32).reshape(4, 3, 2).astype(np.complex64))
                # Different UVW, signed zero, and a flagged nonfinite payload survive.
                t.putcell('UVW', 2, np.array([3., -0., 5.]))
                a = t.getcell('DATA', 3); a.real[1, 1] = np.nan; t.putcell('DATA', 3, a)
                f = t.getcell('FLAG', 3); f[1, 1] = True; t.putcell('FLAG', 3, f)
            before = snapshot(source)
            self.assertEqual(context_runs(source)[1], [(0, 2), (2, 4)])
            bundle = root / 'data.bundle'
            encode_partitioned(source, bundle)
            self.assertEqual(snapshot(source), before)
            # Removing source access must not affect decoding.
            moved = root / 'hidden.ms'; source.rename(moved)
            decode_partitioned(bundle, root / 'restored.ms')
            self.assertTrue(compare_ms(moved, root / 'restored.ms'))
            path = bundle / 'manifest.json'; manifest = json.loads(path.read_text())
            manifest['partitions'][1]['start'] = 1
            path.write_text(json.dumps(manifest))
            (bundle/'manifest.sha256').write_text(file_sha256(path)+'\n')
            with self.assertRaisesRegex(ValueError, 'mapping'):
                decode_partitioned(bundle, root / 'corrupt.ms')
            self.assertFalse((root / 'corrupt.ms').exists())

    def test_quantized_context_split_and_flag_preservation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = make_ms(root / 'source.ms', pairs=((0, 1), (1, 0), (0, 1), (1, 0)), fields=[0, 0, 1, 1])
            with table(str(source), readonly=False, ack=False) as t:
                t.putcol('DATA', np.random.default_rng(2).normal(size=(4, 3, 2)).astype(np.complex64))
                flags = t.getcol('FLAG'); flags[2, 0, 0] = True; t.putcol('FLAG', flags)
            cfg = Config(mode='adaptive', max_component_mse=.01, max_complex_error=.5,
                         codebook='power-midrise-v1', baseline_group=2)
            bundle = root / 'data.bundle'
            encode_partitioned(source, bundle, cfg)
            decode_partitioned(bundle, root / 'restored.ms')
            checks = compare_ms(source, root / 'restored.ms', approximate_columns=['DATA'])
            self.assertTrue(checks)
            with table(str(source), ack=False) as a, table(str(root / 'restored.ms'), ack=False) as b:
                self.assertEqual(a.getcell('DATA', 2)[0, 0], b.getcell('DATA', 2)[0, 0])
            part = bundle / 'partition-0001.h5'
            with part.open('ab') as stream: stream.write(b'corrupted')
            with self.assertRaisesRegex(ValueError, 'checksum'):
                decode_partitioned(bundle, root / 'bad.ms')


if __name__ == '__main__':
    unittest.main()
