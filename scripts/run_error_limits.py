"""Run the preselected matched-rate plan through archives, export and CLEAN.

Two pipelines may overlap on the shared host. Timings are recorded, not used
as isolated speed comparisons. Completed stages can be resumed; partial outputs
are retained for diagnosis and never passed off as complete.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,default=Path('results/calibrated-error-limits/screen-final/plan.json'))
    p.add_argument('--root',type=Path,default=Path('results/calibrated-error-limits'))
    p.add_argument('--source',type=Path,default=Path('data/working/lwa-15slot-flat.ms'))
    p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--cases',nargs='*')
    stages=p.add_mutually_exclusive_group()
    stages.add_argument('--repeat-only',action='store_true',help='Repeat CLEAN on already verified exports')
    stages.add_argument('--image-only',action='store_true',help='First CLEAN on already verified exports')
    stages.add_argument('--archive-only',action='store_true',help='Encode, verify/export and measure visibility errors')
    args=p.parse_args()
    if not 1<=args.jobs<=2:raise ValueError('Use one or two pipelines')
    if args.root.resolve().is_relative_to(args.source.resolve()):
        raise ValueError('Experiment outputs must be outside the source MS')
    args.root.mkdir(parents=True,exist_ok=True)
    plan=json.loads(args.plan.read_text())
    if not isinstance(plan,list) or not plan or any(not isinstance(e,dict) or
            not re.fullmatch(r'[a-z0-9][a-z0-9-]*',e.get('name','')) for e in plan):
        raise ValueError('Plan requires simple, nonempty case names')
    if len({e['name'] for e in plan})!=len(plan):raise ValueError('Duplicate case names')
    if args.cases:
        if set(args.cases)-{x['name'] for x in plan}:raise ValueError('Unknown case')
        plan=[x for x in plan if x['name'] in args.cases]
    env={**os.environ,'PYTHONPATH':'src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MPLCONFIGDIR':'/tmp/varquanms-mpl'}
    def run(case):
        name=case['name']; base=args.root/name
        cfg=args.root/(name+'-config.json')
        if cfg.exists() and json.loads(cfg.read_text())!=case['config']:raise ValueError('Existing configuration differs')
        cfg.write_text(json.dumps(case['config'],indent=2)+'\n')
        def stage(label,script,options,completed):
            if completed.exists():
                data=json.loads(completed.read_text())
                if data.get('status') in ('failed','running'):raise ValueError('Prior stage failed')
                return
            stamp=time.perf_counter()
            command=[sys.executable,'scripts/'+script,*map(str,options)]
            with (args.root/(name+'-'+label+'.log')).open('a') as log:
                result=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
            print(json.dumps(dict(case=name,stage=label,seconds=time.perf_counter()-stamp,returncode=result.returncode)),flush=True)
            if result.returncode or not completed.exists():raise RuntimeError(f'{name}: {label} failed; see stage log')
        bundle=Path(str(base)+'.bundle');shared=Path(str(base)+'.h5');ms=Path(str(base)+'.ms')
        report=lambda stage:args.root/(name+'-'+stage+'.json')
        if args.repeat_only or args.image_only:
            verified=json.loads(report('export').read_text())
            if verified['status']!='passed':raise ValueError('Repeat requires verified export')
            images=args.root/(name+('-repeat-images' if args.repeat_only else '-images'))
            label='repeat-comparison' if args.repeat_only else 'comparison'
            stage('image-repeat' if args.repeat_only else 'image','run_wsclean.py',[ms,'--output',images,'--data-column','CORRECTED_DATA',
                  '--clean-profile','calibrated','--threads','8','--beam-from',
                  'results/calibrated-template/reference/image-image.fits'],images/'complete.json')
            stage(label,'compare_wsclean.py',['results/calibrated-template/reference',images,
                  '--mask','results/calibrated-template/above10deg.fits','--report',report(label)],report(label))
            return name
        stage('encode','timeslot_archive.py',['encode',args.source,bundle,'--config',cfg,
              '--partition-workers',case.get('partition_workers',1),'--report',report('encode')],report('encode'))
        if 'payload_fraction_target' in case:
            manifest=json.loads((bundle/'manifest.json').read_text())
            payload=samples=0
            for part in manifest['partitions']:
                for product in part['summary']['visibility_block_statistics']['CORRECTED_DATA']['products']:
                    payload+=product['serialized_block_bytes']
                    samples+=sum(product['complex_sample_bit_width_histogram'].values())
            rate=dict(payload_bytes=payload,original_complex64_bytes=8*samples,
                      fraction=payload/(8*samples),target=case['payload_fraction_target'],
                      within_target=payload<=case['payload_fraction_target']*8*samples)
            report('rate').write_text(json.dumps(rate,indent=2)+'\n')
            if not rate['within_target']:
                raise ValueError(f'{name}: actual payload exceeds target; retain archive and retune before imaging')
        stage('compact','timeslot_archive.py',['compact',bundle,shared,'--report',report('compact')],report('compact'))
        if 'payload_fraction_target' in case:
            compact=json.loads(report('compact').read_text())
            categories=compact['stored_bytes_by_category']
            rate['physical_payload_bytes']=categories['visibility_block_payload']
            rate['physical_payload_within_target']=rate['physical_payload_bytes']<=rate['target']*rate['original_complex64_bytes']
            rate['payload_plus_permutations_and_indexes_bytes']=sum(categories[k] for k in
                ('visibility_block_payload','permutations','visibility_indexes'))
            report('rate').write_text(json.dumps(rate,indent=2)+'\n')
            if not rate['physical_payload_within_target']:
                raise ValueError(f'{name}: physical payload exceeds target')
        stage('export','validate_partition_export.py',[shared,args.source,ms,'--report',report('export')],report('export'))
        stage('visibility','analyze_visibility_errors.py',[args.source,ms,'--report',report('visibility')],report('visibility'))
        if args.archive_only:return name
        images=args.root/(name+'-images')
        stage('image','run_wsclean.py',[ms,'--output',images,'--data-column','CORRECTED_DATA','--clean-profile','calibrated',
              '--threads','8','--beam-from','results/calibrated-template/reference/image-image.fits'],images/'complete.json')
        stage('comparison','compare_wsclean.py',['results/calibrated-template/reference',images,'--mask',
              'results/calibrated-template/above10deg.fits','--report',report('comparison')],report('comparison'))
        return name
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for name in pool.map(run,plan):print(json.dumps(dict(case=name,status='complete')),flush=True)


if __name__=='__main__':main()
