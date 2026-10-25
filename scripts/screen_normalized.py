"""Screen dimensionless epsilon combinations and block shapes at a byte ceiling.

Each actual packed candidate is encoded once, then its measured error/byte
frontier is reused for coefficient searches. Finalists still require complete
archives and source/export verification before CLEAN. No image selects a pilot
threshold. This study driver uses the calibrated 48-channel observation.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
import hashlib
import itertools
import json
import math
from pathlib import Path
import time

import numpy as np

from varquanms import codec
from varquanms.config import Config
from varquanms.error_budget import EPSILONS, block_references, normalized_limits
from varquanms.power import _parameters
from screen_error_limits import load_blocks

LAYOUTS=[('amplitude128x48',128,48,'amplitude','max-std-v1'),
         ('baseline1x48',1,48,'original','max-std-v1'),
         ('amplitude8x24',8,24,'amplitude','max-std-v1'),
         ('amplitude32x24',32,24,'amplitude','max-std-v1'),
         ('amplitude32x48',32,48,'amplitude','max-std-v1'),
         ('amplitude128x16',128,16,'amplitude','max-std-v1'),
         ('amplitude128x24',128,24,'amplitude','max-std-v1'),
         ('amplitude256x48',256,48,'amplitude','max-std-v1'),
         ('search128x48',128,48,'amplitude','search')]


def recipes():
    result=[];seen=set()
    def add(name,weights):
        identity=tuple(weights.get(k) for k in 'GSDR')
        if identity not in seen:
            result.append(dict(name=name,weights=weights));seen.add(identity)
    for size in range(1,5):
        for subset in itertools.combinations('GSDR',size):
            name=''.join(subset);add(name,{k:1. for k in subset})
    for a,b in itertools.combinations('GSDR',2):
        for ratio in (.25,4.):add(f'{a}{b}-{a}{ratio:g}',{a:ratio,b:1.})
    for g,d in itertools.product((.25,1.,4.),(.5,1.,2.)):
        add(f'GSDR-G{g:g}-D{d:g}',dict(G=g,S=1.,D=d,R=1.))
    add('GSDR-relative-strong',dict(G=1.,S=4.,D=4.,R=1.))
    add('GSDR-spread-strong',dict(G=4.,S=1.,D=2.,R=.25))
    return result


def policy_for(base,recipe,tau):
    return {**base.normalized_bounds,**{k:None for k in EPSILONS},
            **{'epsilon_'+k:float(tau*v) for k,v in recipe['weights'].items()}}


def frontier(job,config,references):
    a=np.ascontiguousarray(job['values'],dtype='<c8');f=np.ascontiguousarray(job['flags'],dtype=bool)
    source=references['groups'][job['configuration_key']]['columns']['CORRECTED_DATA'][job['product']]
    global_std=source['auto' if job['auto'] else 'cross']['component_std']
    stats=block_references(a,f)
    if not stats['valid_samples']:global_std=0.
    _,_,_,widths,params=_parameters(a,f,config)
    options=[];seen=set()
    compression,body=codec._compress(a.tobytes(),config)
    exact=codec.HEADER.pack(0,compression,0,a.size,0,0.)+body
    metrics=codec.error_metrics(a,a,~f&np.isfinite(a),percentile=False)
    options.append(dict(bytes=len(exact),bits=0,metrics=metrics,required_threshold=0.,gamma=None))
    seen.add(exact)
    for width,param in zip(widths,params):
        blob,decoded,_,_=codec._native.power_search(a,f,np.asarray([width],dtype='u1'),
            np.asarray([param],dtype='<f4'),np.full(4,np.inf,dtype='<f8'),a.shape[1],
            config.compression_level if config.compression=='gzip' else -1)
        if blob in seen:continue
        seen.add(blob)
        restored=np.frombuffer(decoded,dtype='<c8').reshape(a.shape)
        m=codec.error_metrics(a,restored,~f&np.isfinite(a),percentile=False)
        p=config.normalized_bounds
        required=max(math.sqrt(m['component_mse']),m['max_complex_error']/p['max_error_factor'],
            m['complex_bias_abs']/p['bias_factor'],math.sqrt(m['max_baseline_component_mse']/p['baseline_mse_factor']))
        bits=codec.HEADER.unpack_from(blob)[0]
        options.append(dict(bytes=len(blob),bits=bits,metrics=m,required_threshold=required,
                            gamma=float(param[3]) if bits else None))
    # Stable payload order; exact fallback wins equal-size ties as in production.
    options.sort(key=lambda c:c['bytes'])
    return dict(**{k:v for k,v in job.items() if k not in ('values','flags')},
                samples=a.size,global_std=global_std,references=stats,candidates=options)


class Matrix:
    def __init__(self,records):
        self.records=records;n=len(records);m=max(len(r['candidates']) for r in records)
        self.required=np.full((n,m),np.inf);self.size=np.full((n,m),np.iinfo('i8').max,dtype='i8')
        self.mse=np.zeros((n,m));self.bits=np.zeros((n,m),dtype='i4')
        self.samples=np.array([r['samples'] for r in records]);self.valid=np.array([r['references']['valid_samples'] for r in records])
        self.parallel=np.array([not r['auto'] and r['product'] in ('XX','YY') for r in records])
        self.refs={k:np.array([r['global_std'] if k=='G' else r['references'][k] for r in records]) for k in 'GSDR'}
        for i,r in enumerate(records):
            for j,c in enumerate(r['candidates']):
                self.required[i,j]=c['required_threshold'];self.size[i,j]=c['bytes']
                self.mse[i,j]=c['metrics']['component_mse'];self.bits[i,j]=c['bits']
        self.original_bytes=int(8*self.samples.sum())

    def evaluate(self,recipe,tau):
        ref=np.minimum.reduce([v*self.refs[k] for k,v in recipe['weights'].items()])
        good=self.required<=tau*ref[:,None]
        sizes=np.where(good,self.size,np.iinfo('i8').max)
        choices=np.argmin(sizes,axis=1);rows=np.arange(len(choices))
        chosen_sizes=sizes[rows,choices];chosen_mse=self.mse[rows,choices];chosen_bits=self.bits[rows,choices]
        assert np.all(chosen_sizes<np.iinfo('i8').max)
        total=int(chosen_sizes.sum());weights=self.valid*self.parallel
        return dict(tau=float(tau),payload_bytes=total,payload_fraction=total/self.original_bytes,
            cross_parallel_component_mse=float(np.sum(chosen_mse*weights)/max(1,weights.sum())),
            component_mse=float(np.sum(chosen_mse*self.valid)/max(1,self.valid.sum())),
            sample_bits={str(int(b)):int(self.samples[chosen_bits==b].sum()) for b in np.unique(chosen_bits)},
            choices=choices.tolist())


def select_threshold(matrix,recipe,target):
    lo,hi=1e-7,16.;trace=[]
    upper=matrix.evaluate(recipe,hi);trace.append(upper)
    if upper['payload_fraction']>target:
        return dict(status='infeasible_at_tested_maximum',**{k:v for k,v in upper.items() if k!='choices'}),trace
    for _ in range(16):
        mid=math.sqrt(lo*hi);trial=matrix.evaluate(recipe,mid);trace.append(trial)
        if trial['payload_fraction']>target:lo=mid
        else:hi=mid
    chosen=min((t for t in trace if t['payload_fraction']<=target),key=lambda t:t['tau'])
    return dict(status='feasible',**chosen),trace


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--references',type=Path,required=True)
    p.add_argument('--target',type=float,default=.16);p.add_argument('--pilot-target',type=float,default=.1575)
    p.add_argument('--workers',type=int,default=16)
    p.add_argument('--layouts',nargs='*')
    args=p.parse_args();args.output.mkdir(exist_ok=True)
    if not 0<args.pilot_target<args.target<1:raise ValueError('Invalid byte targets')
    if codec._native is None or not getattr(codec._native,'THREE_BIT_V1',False):
        raise ValueError('This screening driver requires the tested native power encoder')
    base=Config.load('configs/calibrated-normalized.json')
    references=json.loads(args.references.read_text())
    assert references['sampling']['fraction']==base.normalized_bounds['sampling_fraction']
    assert references['sampling']['seed']==base.normalized_bounds['sampling_seed']
    all_recipes=recipes();entries=[];layout_reports=[];started=time.perf_counter()
    layouts=LAYOUTS if not args.layouts else [x for x in LAYOUTS if x[0] in args.layouts]
    if args.layouts and set(args.layouts)-{x[0] for x in LAYOUTS}:raise ValueError('Unknown layouts')
    for layout,group,segment,ordering,gamma in layouts:
        cfg=replace(base,baseline_group=group,channel_segment=segment,ordering=ordering,power_gamma_policy=gamma)
        jobs=load_blocks(args.source,cfg,groups=64 if group>1 else 512)
        cache=args.output/(layout+'-frontier.json')
        if cache.exists():
            cached=json.loads(cache.read_text());assert cached['config']==cfg.to_dict()
            assert cached['references_sha256']==hashlib.sha256(args.references.read_bytes()).hexdigest()
            records=cached['blocks']
        else:
            stamp=time.perf_counter()
            with ThreadPoolExecutor(max_workers=args.workers) as pool:
                records=list(pool.map(lambda j:frontier(j,cfg,references),jobs))
            cached=dict(config=cfg.to_dict(),references_sha256=hashlib.sha256(args.references.read_bytes()).hexdigest(),
                blocks=records,seconds=time.perf_counter()-stamp)
            cache.write_text(json.dumps(cached,indent=2,allow_nan=False)+'\n')
        matrix=Matrix(records);trials=[];local=[]
        for recipe in all_recipes:
            chosen,trace=select_threshold(matrix,recipe,args.pilot_target)
            entry=dict(layout=layout,family=recipe['name'],weights=recipe['weights'],**chosen)
            entry.pop('choices',None)
            if chosen['status']=='feasible':
                final=replace(cfg,normalized_bounds=policy_for(cfg,recipe,chosen['tau']))
                final.validate();entry['config']=final.to_dict()
                # Check the cached candidate decision against the actual encoder
                # at auto, cross, partial and special groups for every recipe.
                for index in sorted({0,len(jobs)//3,len(jobs)//2,len(jobs)-1}):
                    job=jobs[index];expected=records[index]['candidates'][chosen['choices'][index]]
                    actual=codec.encode_block(job['values'],job['flags'],final,final.bounds(),global_std=records[index]['global_std'])
                    if len(actual.blob)!=expected['bytes']:
                        raise ValueError(f'Frontier/encoder disagreement: {layout}/{recipe["name"]}/{index}')
                entry['frontier_encoder_checks']=4
            entries.append(entry);local.append(entry)
            trials.append(dict(family=recipe['name'],trials=[{k:v for k,v in t.items() if k!='choices'} for t in trace]))
        report=dict(layout=layout,config=cfg.to_dict(),blocks=len(records),
            original_sample_bytes=matrix.original_bytes,candidate_frontier_seconds=cached['seconds'],
            experiments=local,traces=trials)
        (args.output/(layout+'-screen.json')).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
        layout_reports.append({k:v for k,v in report.items() if k!='traces'})
        print(json.dumps(dict(layout=layout,blocks=len(records),recipes=len(local),
            feasible=sum(e['status']=='feasible' for e in local),seconds=cached['seconds'])),flush=True)
    feasible=[e for e in entries if e['status']=='feasible'];selected=[]
    def nominate(options,reason):
        options=[e for e in options if e not in selected]
        if not options:return
        item=min(options,key=lambda e:e['cross_parallel_component_mse'])
        item['nomination']=reason;selected.append(item)
    main_cases=[e for e in feasible if e['layout']=='amplitude128x48']
    for family in 'GSDR':nominate([e for e in main_cases if e['family']==family],'single-reference control')
    nominate([e for e in main_cases if len(e['weights'])==4],'best four-reference main pilot')
    nominate([e for e in main_cases if 1<len(e['weights'])<4],'best pair/triple main pilot')
    nominate([e for e in feasible if e['layout']=='baseline1x48'],'mandatory unsorted per-baseline control')
    nominate([e for e in feasible if e['layout']=='amplitude32x24'],'smaller grouped control')
    nominate([e for e in feasible if e['layout'] not in ('amplitude128x48','baseline1x48','amplitude32x24','search128x48')],
             'best other block-shape pilot')
    nominate([e for e in feasible if e['layout']=='search128x48'],'full gamma-search alternative')
    plan=[]
    for index,e in enumerate(selected):
        name=f'eps16-{index:02d}-{e["layout"]}-{e["family"].lower().replace(".","p")}'
        plan.append(dict(name=name,family=e['family'],layout=e['layout'],weights=e['weights'],
            parameter=e['tau'],nomination=e['nomination'],payload_fraction_target=args.target,
            partition_workers=6 if e['config']['baseline_group']==1 else 1,
            pilot_payload_fraction=e['payload_fraction'],config=dict(e['config'],
                workers=1 if e['config']['baseline_group']==1 else e['config']['workers'])))
    (args.output/'plan.json').write_text(json.dumps(plan,indent=2,allow_nan=False)+'\n')
    result=dict(schema='varquanms.normalized-screen/1',source=str(args.source),target=args.target,
        pilot_target=args.pilot_target,references_path=str(args.references),references=references,
        recipes=all_recipes,layouts=layout_reports,experiments=entries,selected=plan,
        seconds=time.perf_counter()-started,
        selection='Payload ceiling determines tau; measured pilot cross XX/YY MSE nominates alternatives. No image used for threshold selection.',
        limitations=['Different layouts sample different rank blocks; pilot MSE is not an imaging ranking.',
            'One shared tau multiplies the recorded dimensionless recipe weights; zero-reference blocks require exact values.',
            'Each actual candidate includes packed codes, factors, headers, exceptions and secondary compression; whole archive overhead is separate.',
            'This study uses first/middle/last integrations and requires a single 48-channel configuration within each integration.'])
    (args.output/'screen.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status='complete',experiments=len(entries),finalists=len(plan))),flush=True)


if __name__=='__main__':main()
