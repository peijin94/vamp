"""Matched-block max/std gamma experiment on first/middle/last integrations."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
import json
from pathlib import Path
import time

import numpy as np

from varquanms.codec import encode_block
from varquanms.config import Config
from varquanms.inspection import _open, CORRELATIONS
from varquanms.ordering import all_orders
from varquanms.partitioned import context_runs
from varquanms.power import maximum_std_ratio


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('ms',type=Path); p.add_argument('output',type=Path)
    p.add_argument('--groups-per-product',type=int,default=64)
    args=p.parse_args()
    if not 1 <= args.groups_per_product <= 4096:
        raise ValueError('groups-per-product must be in [1,4096]')
    args.output.mkdir(exist_ok=False)
    base=Config.load('configs/calibrated-test-amplitude256.json')
    _,runs=context_runs(args.ms)
    selected_times=sorted({0,len(runs)//2,len(runs)-1})
    jobs=[]
    for integration in selected_times:
        start,stop=runs[integration]
        with _open(args.ms) as t:
            data=t.getcol('CORRECTED_DATA',start,stop-start)
            flags=t.getcol('FLAG',start,stop-start)|t.getcol('FLAG_ROW',start,stop-start)[:,None,None]
            uvw=t.getcol('UVW',start,stop-start)
            pairs=np.column_stack([t.getcol(c,start,stop-start) for c in ('ANTENNA1','ANTENNA2')])
            ddids=np.unique(t.getcol('DATA_DESC_ID',start,stop-start))
        if len(ddids)!=1: raise ValueError('Screen requires one configuration per context')
        if data.shape[1] != base.channel_segment:
            raise ValueError('This matched-block experiment requires 48 channels per DATA cell')
        with _open(args.ms/'DATA_DESCRIPTION') as t: polid=t.getcell('POLARIZATION_ID',int(ddids[0]))
        with _open(args.ms/'POLARIZATION') as t: names=[CORRELATIONS[int(v)] for v in t.getcell('CORR_TYPE',int(polid))]
        data[flags]=0j
        orders,_,_=all_orders(data,flags,pairs,uvw,base)
        for product,order in enumerate(orders):
            bounds=order['group_boundaries']
            indices=np.unique(np.linspace(0,len(bounds)-2,min(args.groups_per_product,len(bounds)-1)).astype(int))
            for index in indices:
                rows=order['permutation'][bounds[index]:bounds[index+1]]
                a=data[rows,:,product];f=flags[rows,:,product]
                components=np.stack([a.real,a.imag],axis=-1).astype('f8')
                protected=f|~np.isfinite(a)|np.any((components==0)&np.signbit(components),axis=-1)
                eligible=components[~protected]
                center=np.median(eligible,axis=0).astype('f4') if len(eligible) else np.zeros(2,dtype='f4')
                stat=maximum_std_ratio(eligible-center.astype('f8'))
                jobs.append(dict(integration=integration,product=names[product],group=int(index),values=a,flags=f,statistics=stat))
    policies=[('symmetric',dict(codebook='uniform-symmetric-v1')),
              ('median-uniform',dict(codebook='power-midrise-v1',power_gammas=[1.])),
              ('power-search2',dict(codebook='power-midrise-v1',power_gammas=[1.,.5])),
              ('power-search3',dict(codebook='power-midrise-v1',power_gammas=[1.,.5,.25]))]
    for reference in (3.,4.,6.):
        for policy in ('max-std-v1','max-std-inlier-v1'):
            policies.append((f'{policy}-r{reference:g}',dict(codebook='power-midrise-v1',power_gammas=[1.,.5,.25],
                                                           power_gamma_policy=policy,power_ratio_reference=reference)))
    records=[]
    for budget in (64.,256.):
        for label,options in policies:
            cfg=replace(base,**options,max_component_mse=budget,max_baseline_component_mse=2*budget,
                        max_complex_error=float(10*np.sqrt(budget)),max_complex_bias=float(np.sqrt(budget)),workers=8)
            cfg.validate()
            stamp=time.perf_counter()
            def run(job):
                block=encode_block(job['values'],job['flags'],cfg,cfg.bounds(job['product']))
                return dict(integration=job['integration'],product=job['product'],group=job['group'],
                            samples=block.decoded.size,bytes=len(block.blob),bits=block.bits,
                            gamma=block.parameters.get('gamma'),metrics=block.metrics,timings=block.timings,
                            exceptions=block.exceptions)
            with ThreadPoolExecutor(max_workers=8) as pool: blocks=list(pool.map(run,jobs))
            elapsed=time.perf_counter()-stamp
            count=sum(b['metrics']['valid_samples'] for b in blocks)
            record=dict(name=f'{label}-b{budget:g}',config=cfg.to_dict(),sampled_blocks=len(blocks),
                        payload_bytes=sum(b['bytes'] for b in blocks),
                        component_mse=sum(b['metrics']['component_mse']*b['metrics']['valid_samples'] for b in blocks)/count,
                        max_complex_error=max(b['metrics']['max_complex_error'] for b in blocks),
                        block_bits=dict(Counter(str(b['bits']) if b['bits'] else 'float32' for b in blocks)),
                        gamma_blocks=dict(Counter(str(b['gamma']) for b in blocks if b['gamma'] is not None)),
                        elapsed_seconds=elapsed,blocks=blocks)
            records.append(record)
            print(json.dumps({k:v for k,v in record.items() if k not in ('config','blocks')}),flush=True)
    report=dict(schema='varquanms.gamma-ratio-screen/1',source=str(args.ms),integrations=selected_times,
                blocks=[{k:v for k,v in job.items() if k not in ('values','flags')} for job in jobs],
                experiments=records,
                limitations=['Same exact sampled blocks for all profiles; payload bytes exclude whole-archive metadata.',
                             'Population std is signal spread, not noise; max/std depends on sample count and outliers also inflate std.',
                             'Wall times use eight workers on a shared host without cache eviction.'])
    (args.output/'screen.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')


if __name__=='__main__': main()
