"""Select imaging finalists at a specified visibility-payload fraction.

Defaults reproduce the earlier 24.8% pilot and 25% ceiling; command-line options
set other rates and candidate widths. Grouping pilots sample different blocks,
so their metrics nominate finalists rather than establish an image ranking.
"""
import argparse
from dataclasses import replace
import json
from pathlib import Path

import numpy as np

from screen_error_limits import configuration, evaluate, load_blocks
from vamp.config import Config


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('ms',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--target-fraction',type=float,default=.25)
    p.add_argument('--pilot-fraction',type=float)
    p.add_argument('--bits',type=int,nargs='+',help='Explicit candidate widths; omit to reproduce the stored base configuration')
    p.add_argument('--layouts',nargs='+')
    p.add_argument('--max-finalists',type=int,choices=(4,5),default=5)
    p.add_argument('--prefix',default='quarter')
    args=p.parse_args()
    if not 0<args.target_fraction<1:raise ValueError('Target fraction must be between zero and one')
    pilot=args.pilot_fraction if args.pilot_fraction is not None else .992*args.target_fraction
    if not 0<pilot<args.target_fraction:raise ValueError('Pilot fraction must be below the target')
    args.output.mkdir(exist_ok=False)
    base=Config.load('configs/calibrated-gamma-ratio.json')
    if args.bits is not None:base=replace(base,bits=args.bits)
    base.validate()
    layouts=[('amplitude128x48',128,48,'amplitude','max-std-v1'),
             ('amplitude32x24',32,24,'amplitude','max-std-v1'),
             ('amplitude8x48',8,48,'amplitude','max-std-v1'),
             ('baseline1x48',1,48,'original','max-std-v1'),
             ('search128x48',128,48,'amplitude','search')]
    if args.layouts:
        if set(args.layouts)-{x[0] for x in layouts}:raise ValueError('Unknown layout')
        layouts=[x for x in layouts if x[0] in args.layouts]
        if 'amplitude128x48' not in args.layouts or len(layouts)<2:raise ValueError('Require main layout and an alternative')
    entries=[];trace=[];selected=[]
    for layout,group,segment,ordering,gamma in layouts:
        cfg=replace(base,baseline_group=group,channel_segment=segment,ordering=ordering,power_gamma_policy=gamma)
        jobs=load_blocks(args.ms,cfg,groups=64 if group>1 else 512)
        pivot=float(np.median([j['signal_component_power'] for j in jobs if not j['auto'] and j['signal_component_power']>0]))
        total=sum(j['values'].size for j in jobs)
        target=pilot*8*total
        for family in ('absolute','relative','hybrid'):
            low,high=((1e-6,1.) if family=='relative' else (1e-4,1e8))
            trials=[]
            for _ in range(14):
                value=float(np.sqrt(low*high));c=configuration(cfg,family,value,pivot)
                r=evaluate(jobs,c)
                parallel=[b for j,b in zip(jobs,r['blocks']) if not j['auto'] and j['product'] in ('XX','YY')]
                count=sum(b['metrics']['valid_samples'] for b in parallel)
                r.update(layout=layout,family=family,parameter=value,sampled_samples=total,
                    payload_ratio=r['payload_bytes']/(8*total),
                    cross_parallel_component_mse=sum(b['metrics']['component_mse']*b['metrics']['valid_samples'] for b in parallel)/max(1,count),
                    pivot_component_power=pivot)
                r.pop('blocks')
                trials.append(r);trace.append(r)
                if r['payload_bytes']>target:low=value
                else:high=value
            under=[r for r in trials if r['payload_bytes']<=target]
            if not under:
                entries.append(dict(name=layout+'-'+family,layout=layout,family=family,
                    status='no_feasible_trial',minimum_tried_payload_fraction=min(r['payload_ratio'] for r in trials)))
                print(json.dumps(entries[-1]),flush=True)
                continue
            # Tightest tried threshold that fits. Do not choose on image error.
            chosen=min(under,key=lambda r:r['parameter'])
            name=layout+'-'+family
            record=dict(name=name,**chosen);entries.append(record)
            print(json.dumps({k:v for k,v in record.items() if k!='config'}),flush=True)
            if layout=='amplitude128x48':selected.append(record)
    # One additional layout/codec finalist nominated by parallel-hand MSE only.
    alternatives=[e for e in entries if e['layout']!='amplitude128x48' and 'config' in e]
    if alternatives:selected.append(min(alternatives,key=lambda e:e['cross_parallel_component_mse']))
    # Include an actual smaller-group image control even when its visibility
    # MSE is worse: unweighted visibility MSE does not determine Briggs images.
    smaller=[e for e in entries if e['layout'] in ('amplitude32x24','amplitude8x48') and 'config' in e]
    if smaller and args.max_finalists==5:
        smaller=min(smaller,key=lambda e:e['cross_parallel_component_mse'])
        if smaller not in selected:selected.append(smaller)
    plan=[dict(name=args.prefix+'-'+e['name'],family=e['family'],layout=e['layout'],parameter=e['parameter'],
               payload_fraction_target=args.target_fraction,pilot_payload_fraction=e['payload_ratio'],config=e['config']) for e in selected]
    (args.output/'plan.json').write_text(json.dumps(plan,indent=2,allow_nan=False)+'\n')
    (args.output/'screen.json').write_text(json.dumps(dict(schema='vamp.fixed-rate-screen/1',
        source=str(args.ms),payload_fraction_target=args.target_fraction,pilot_fraction_target=pilot,
        layouts=layouts,experiments=entries,trials=trace,selected=plan,
        selection='Feasible error rules at 128x48, the alternative layout/codec with lowest sampled cross XX/YY MSE, and optionally a smaller-group control. No candidate image used for selection.',
        limitations=['Grouping variants sample different original blocks at regular rank intervals. All products/times use original values, but this is a pilot, not a matched-cell image claim.',
                     'Payload includes packed codes, stored factors, exact exceptions, headers and conventional compression. Whole archive overhead is separate.',
                     'Only actual whole-data payload within the configured ceiling qualifies. Final image comparisons use two CLEAN runs and report dirty/bright/faint/peak metrics.']),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
