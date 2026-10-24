import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import unittest

from casacore.tables import table
import h5py
import numpy as np

from varquanms.archive import compare_ms, decode
from varquanms.config import Config
from varquanms.container import Reader
from varquanms.metadata import dumps
from varquanms.partitioned import encode_partitioned, decode_partitioned, partition_reader
from varquanms.shared import compact_bundle, load_shared
from test_inspection import make_ms, snapshot
from test_calibrated import calibrated_fixture


class SharedTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def make_bundle(self, calibrated=False):
        source = self.root/'source.ms'
        if calibrated:
            calibrated_fixture(source)
            cfg = Config(mode='adaptive', columns=['CORRECTED_DATA'], keep_calib_only=True,
                         max_component_mse=.01, max_complex_error=.5, codebook='power-midrise-v1')
        else:
            make_ms(source)
            cfg = Config()
        with table(str(source), ack=False, readonly=False) as t:
            t.putcol('TIME', np.array([1., 1., 2.]))
            if not calibrated:
                values = t.getcol('DATA')
                values.real[0, 0, 0] = -0.
                values.imag[2, 0, 0] = np.inf
                values.view(np.uint32).reshape(-1)[3] = 0x7fc00012  # Preserve a NaN payload.
                t.putcol('DATA', values)
                t.putcell('UVW', 2, np.array([-0., 1., 3.]))
        bundle = self.root/'original.bundle'
        encode_partitioned(source, bundle, cfg)
        return source, bundle

    def test_exact_independent_round_trip_and_shared_storage(self):
        source, bundle = self.make_bundle()
        before = snapshot(source)
        compact = self.root/'shared.h5'
        result = compact_bundle(bundle, compact)
        self.assertEqual(sum(result['stored_bytes_by_category'].values()), compact.stat().st_size)
        self.assertGreater(result['storage_compaction']['counts']['shared_array_references'], 0)
        original = json.loads((bundle/'manifest.json').read_text())
        shared = load_shared(compact)
        # Partition scientific manifests and every logical-array digest are identical.
        for a, b in zip(original['partitions'], shared['partitions']):
            with Reader(bundle/a['archive']) as x, partition_reader(compact, b) as y:
                self.assertEqual(x.h5['manifest'][...].tobytes(), y.h5['manifest'][...].tobytes())
                self.assertEqual(x.entries, y.entries)
        self.assertEqual(snapshot(source), before)
        source.rename(self.root/'hidden.ms')
        shutil.rmtree(bundle)
        decode_partitioned(compact, self.root/'restored.ms')
        self.assertTrue(compare_ms(self.root/'hidden.ms', self.root/'restored.ms'))
        # A later partition is directly addressable without decoding earlier ones.
        decode(compact, self.root/'last.ms', group=shared['partitions'][-1]['group'])
        with table(str(self.root/'last.ms'), ack=False) as t:
            self.assertEqual(t.nrows(), 1)
            self.assertEqual(t.getcell('TIME', 0), 2.)

    def test_calibrated_quantization_export_is_bitwise_unchanged(self):
        source, bundle = self.make_bundle(calibrated=True)
        decode_partitioned(bundle, self.root/'ordinary.ms')
        result = compact_bundle(bundle, self.root/'shared.h5')
        self.assertGreater(result['storage_compaction']['counts']['payload_allocation_bytes_removed'], 0)
        shared = load_shared(self.root/'shared.h5')
        for part in shared['partitions']:
            with partition_reader(self.root/'shared.h5', part) as reader:
                for entry in reader.entries.values():
                    if entry['role'] == 'visibility_block_payload':
                        ds = reader.h5[entry['path']]
                        self.assertEqual(ds.id.get_storage_size(), entry['logical_bytes'])
        decode_partitioned(self.root/'shared.h5', self.root/'shared.ms')
        self.assertTrue(compare_ms(self.root/'ordinary.ms', self.root/'shared.ms'))
        with table(str(self.root/'shared.ms'), ack=False) as t:
            self.assertNotIn('MODEL_DATA', t.colnames())
            self.assertEqual(t.getcol('DATA').tobytes(), t.getcol('CORRECTED_DATA').tobytes())
        # Contiguous payloads retain codec CRCs and inventory SHA256 protection.
        with h5py.File(self.root/'shared.h5', 'r+') as h5:
            group = h5[shared['partitions'][0]['group']]
            spec = json.loads(group['manifest'][...].tobytes())
            entry = next(e for e in spec['inventory'] if e['role'] == 'visibility_block_payload')
            group[entry['path']][0] ^= np.uint8(1)
        with self.assertRaisesRegex(ValueError, 'checksum'):
            load_shared(self.root/'shared.h5')

    def test_corruption_and_unsafe_links_fail_before_export(self):
        _, bundle = self.make_bundle()
        compact = self.root/'shared.h5'
        compact_bundle(bundle, compact)
        for failure in ('payload', 'mapping', 'external', 'virtual', 'cycle', 'missing', 'verification', 'incomplete', 'truncated'):
            with self.subTest(failure=failure):
                damaged = self.root/(failure+'.h5')
                shutil.copyfile(compact, damaged)
                with h5py.File(damaged, 'r+') as h5:
                    m = json.loads(h5['manifest'][...].tobytes())
                    group = m['partitions'][0]['group']
                    if failure == 'payload':
                        spec = json.loads(h5[group+'/manifest'][...].tobytes())
                        entry = next(e for e in spec['inventory'] if e['role']=='subtables' and e['logical_bytes'])
                        ds = h5[group+'/'+entry['path']]
                        values = ds[...]; values.reshape(-1)[0] = values.reshape(-1)[0]+1
                        ds[...] = values
                    elif failure == 'mapping':
                        m['partitions'][1]['start'] = 0
                        raw = dumps(m)
                        del h5['manifest']
                        h5.create_dataset('manifest', data=np.frombuffer(raw, dtype='u1'))
                        h5.attrs['manifest_sha256'] = hashlib.sha256(raw).hexdigest()
                    elif failure == 'external':
                        h5['escape'] = h5py.ExternalLink('outside.h5', '/')
                    elif failure == 'virtual':
                        layout = h5py.VirtualLayout(shape=(1,), dtype='u1')
                        layout[:] = h5py.VirtualSource('outside.h5', 'data', shape=(1,))
                        h5.create_virtual_dataset('virtual', layout)
                    elif failure == 'cycle':
                        h5['partitions']['cycle'] = h5['partitions']
                    elif failure == 'missing':
                        del h5[group+'/manifest']
                    elif failure == 'verification':
                        ds = h5[m['partitions'][0]['verification_dataset']]
                        ds[0] = 0
                    elif failure == 'incomplete':
                        h5.attrs['complete'] = False
                if failure == 'truncated':
                    with damaged.open('r+b') as stream:
                        stream.truncate(damaged.stat().st_size//2)
                out = self.root/(failure+'.ms')
                with self.assertRaises((ValueError, KeyError, OSError)):
                    decode_partitioned(damaged, out)
                self.assertFalse(out.exists())


if __name__ == '__main__':
    unittest.main()
