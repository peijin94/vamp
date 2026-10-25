"""Dimensionless budgets, reproducible sampling and independent reconstruction."""
from dataclasses import replace
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from casacore.tables import table

from varquanms import codec
from varquanms.archive import encode, verify, compare_ms
from varquanms.config import Config
from varquanms.container import Reader
from varquanms.error_budget import block_references, normalized_limits, normalized_policy
from varquanms.partitioned import encode_partitioned, decode_partitioned, _load
from varquanms.references import sample_references, Sampler, validate_reference_manifest
from varquanms.shared import compact_bundle
from test_inspection import make_ms, snapshot


def configuration(**eps):
    return Config(mode='adaptive',codebook='power-midrise-v1',baseline_group=2,channel_segment=2,
                  normalized_bounds=normalized_policy(**dict(dict(epsilon_G=.02,epsilon_S=.02,
                                                      epsilon_D=.03,epsilon_R=.02),**eps)))


class NormalizedBudgetTests(unittest.TestCase):
    def test_exact_statistics_and_zero_mad(self):
        a=np.array([[1+10j,3+14j,100+200j,complex(float('nan'),0)]],dtype='c8')
        f=np.array([[False,False,True,False]])
        s=block_references(a,f)
        self.assertEqual(s['valid_samples'],2)
        self.assertEqual(s['median_real'],2);self.assertEqual(s['median_imag'],12)
        self.assertAlmostEqual(s['S'],np.sqrt(2.5));self.assertAlmostEqual(s['D'],np.sqrt(2.5))
        self.assertAlmostEqual(s['R'],np.sqrt(76.5))
        p=normalized_policy(epsilon_G=.1,epsilon_S=.2,epsilon_D=.3,epsilon_R=.4)
        b=normalized_limits(p,s,2.)
        self.assertAlmostEqual(b['max_component_mse'],.04)
        self.assertAlmostEqual(b['max_complex_error'],2.)
        sparse=np.zeros((2,32),dtype='c8');sparse[0,0]=100
        s=block_references(sparse,np.zeros(sparse.shape,bool))
        self.assertEqual(s['D'],0);self.assertGreater(s['S'],0)
        self.assertEqual(normalized_limits(p,s,2.)['max_component_mse'],0)
        # No 1e-6 or other absolute floor for all-invalid / all-flagged blocks.
        self.assertEqual(normalized_limits(p,block_references(a,np.ones(a.shape,bool)),None)['max_component_mse'],0)

    def test_gain_scaling_native_portable_and_exact_fallback(self):
        rng=np.random.default_rng(9)
        a=(rng.normal(size=(8,16))+1j*rng.normal(size=(8,16))).astype('c8')
        a[0,0]=complex(-0.,1);f=np.zeros(a.shape,bool);f[0,1]=True
        cfg=replace(configuration(),compression='none',scale_quantiles=[1.])
        choices=[]
        for gain in (2.**-16,1.,2.**16):
            v=a*gain;g=2*gain
            result=codec.encode_block(v,f,cfg,cfg.bounds(),global_std=g)
            with patch.object(codec,'_native',None):
                portable=codec.encode_block(v,f,cfg,cfg.bounds(),global_std=g)
            self.assertEqual(result.blob,portable.blob)
            self.assertTrue(codec.within_bounds(result.metrics,normalized_limits(cfg.normalized_bounds,block_references(v,f),g)))
            choices.append(result.bits)
            self.assertEqual(result.decoded[f].tobytes(),v[f].tobytes())
        self.assertEqual(len(set(choices)),1)
        cfg=replace(cfg,normalized_bounds=normalized_policy(epsilon_D=0.))
        result=codec.encode_block(a,f,cfg,cfg.bounds())
        self.assertEqual(result.decoded.tobytes(),a.tobytes())
        with self.assertRaisesRegex(ValueError,'Missing'):
            codec.encode_block(a,f,configuration(),configuration().bounds())

    def test_unknown_mixed_nonfinite_and_empty_policy_rejected(self):
        cfg=configuration();cfg.validate()
        for p in [dict(reference='distribution-v1'),dict(reference='distribution-v1',epsilon_G=-1),
                  dict(reference='distribution-v1',epsilon_G=float('nan')),
                  dict(reference='distribution-v1',epsilon_G=True),
                  dict(reference='distribution-v1',epsilon_G=.1,epsilon_Z=.1),
                  dict(reference='distribution-v1',epsilon_G=.1,sampling_fraction=0)]:
            with self.assertRaises(ValueError):replace(cfg,normalized_bounds=p).validate()
        for update in [dict(max_component_mse=1),dict(max_complex_error=1),dict(relative_bounds={'reference':'block-rms-v1','max_component_rmse':.1})]:
            with self.assertRaisesRegex(ValueError,'do not mix'):replace(cfg,**update).validate()

    def test_priority_sampling_chunk_independence_minimum_and_gain(self):
        ids=np.arange(10000,dtype='<u8');a=(ids+1j*ids[::-1]).astype('c8')
        results=[]
        for chunk in (19,10000):
            s=Sampler(.02,128,71)
            for start in range(0,len(a),chunk):s.add(ids[start:start+chunk],a[start:start+chunk])
            results.append(s.finish())
        self.assertEqual(results[0]['sample_sha256'],results[1]['sample_sha256'])
        self.assertAlmostEqual(results[0]['component_std'],results[1]['component_std'],places=9)
        s=Sampler(.000001,128,71);s.add(ids,a);r=s.finish()
        self.assertTrue(r['minimum_sample_used']);self.assertEqual(r['sampled_samples'],128)
        s=Sampler(1.,128,71);s.add(ids[:20],a[:20]);r=s.finish()
        self.assertEqual(r['sampled_samples'],20)

    def test_source_strata_sampling_metadata_and_partition_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            pairs=((0,0),(0,1),(1,1),(1,0))*2
            source=make_ms(root/'source.ms',pairs=pairs,channels=5,correlations=(12,9,10,11))
            rng=np.random.default_rng(123)
            with table(str(source),readonly=False,ack=False) as t:
                a=(rng.normal(size=(8,5,4))+1j*rng.normal(size=(8,5,4))).astype('c8')
                a[::2]*=100;a[:,:,2]*=.001;a[0,0,0]=complex(float('nan'),0)
                f=np.zeros(a.shape,bool);f[0,1]=True
                t.putcol('DATA',a);t.putcol('FLAG',f)
                t.putcol('TIME',np.r_[np.ones(4),np.full(4,2.)])
                # Concatenated time slots may have separate observation IDs;
                # unchanged physical configuration must share global references.
                t.putcol('OBSERVATION_ID',np.r_[np.zeros(4,dtype='i4'),np.ones(4,dtype='i4')])
            cfg=configuration();before=snapshot(source)
            refs=sample_references(source,cfg)
            self.assertEqual(len(refs['groups']),1)
            group=next(iter(refs['groups'].values()))['columns']['DATA']
            self.assertEqual(set(group),{'XX','YY','XY','YX'})
            self.assertGreater(group['XX']['auto']['component_std'],50*group['XX']['cross']['component_std'])
            self.assertEqual(group['YY']['auto']['eligible_samples'],18)
            bundle=root/'archive.bundle';encode_partitioned(source,bundle,cfg)
            _,manifest=_load(bundle)
            self.assertEqual(manifest['distribution_references']['groups'],refs['groups'])
            for part in manifest['partitions']:
                with Reader(bundle/part['archive']) as reader:
                    self.assertEqual(reader.manifest['distribution_references'],manifest['distribution_references'])
            compact_bundle(bundle,root/'shared.h5')
            decode_partitioned(root/'shared.h5',root/'restored.ms')
            self.assertTrue(compare_ms(source,root/'restored.ms',approximate_columns=['DATA']))
            self.assertEqual(before,snapshot(source))

    def test_single_archive_and_mixed_configuration_references(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=make_ms(root/'source.ms',pairs=((0,0),(0,1)))
            cfg=configuration();path=root/'single.h5'
            encode(source,path,cfg);self.assertEqual(verify(path,source)['status'],'passed')
            with Reader(path) as reader:
                original=reader.manifest
            for failure in ('missing','scale','count','binding'):
                broken=deepcopy(original)
                if failure=='missing':del broken['distribution_references']
                elif failure in ('scale','count'):
                    g=next(iter(broken['distribution_references']['groups'].values()))
                    item=g['columns']['DATA']['XX']['auto']
                    if failure=='scale':item['component_std']=-1
                    else:item['sampled_samples']=item['eligible_samples']+1
                else:
                    broken['tables'][0]['columns']['DATA']['distribution_reference_key']='invalid'
                with self.assertRaises(ValueError):validate_reference_manifest(broken)
            mixed=make_ms(root/'mixed.ms',pairs=((0,0),(0,1),(0,0),(0,1)),fields=[0,0,1,1])
            refs=sample_references(mixed,cfg)
            self.assertEqual(len(refs['groups']),2)

    def test_parallel_partitions_keep_payloads_and_export_identical(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);source=make_ms(root/'source.ms',pairs=((0,0),(0,1),(0,0),(0,1)),channels=5)
            with table(str(source),readonly=False,ack=False) as t:
                t.putcol('TIME',np.array([1.,1.,2.,2.]))
            before=snapshot(source);cfg=configuration()
            for n in (1,2):
                execution_cfg=cfg if n==1 else replace(cfg,workers=2)
                with patch('varquanms.partitioned.os.cpu_count',return_value=4):
                    encode_partitioned(source,root/f'archive{n}',execution_cfg,partition_workers=n)
                decode_partitioned(root/f'archive{n}',root/f'export{n}.ms')
            self.assertTrue(compare_ms(root/'export1.ms',root/'export2.ms'))
            _,a=_load(root/'archive1');_,b=_load(root/'archive2')
            self.assertEqual(b['execution']['partition_workers'],2)
            for pa,pb in zip(a['partitions'],b['partitions']):
                with Reader(root/'archive1'/pa['archive']) as ra,Reader(root/'archive2'/pb['archive']) as rb:
                    self.assertEqual(ra.manifest['source']['row_window']['start'],pa['start'])
                    self.assertIsNone(pa['summary']['whole_archive_ratio'])
                    ia={e['path']:e for e in ra.manifest['inventory']}
                    ib={e['path']:e for e in rb.manifest['inventory']}
                    self.assertEqual(ia,ib)
                    for path in ia:self.assertEqual(ra.h5[path][...].tobytes(),rb.h5[path][...].tobytes())
            self.assertEqual(before,snapshot(source))
            with patch('varquanms.partitioned.os.cpu_count',return_value=2):
                with self.assertRaisesRegex(ValueError,'CPU budget'):
                    encode_partitioned(source,root/'bad',replace(cfg,workers=2),partition_workers=2)
            with self.assertRaises(ValueError):encode_partitioned(source,root/'bad',cfg,partition_workers=0)


if __name__=='__main__':unittest.main()
