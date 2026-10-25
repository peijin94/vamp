"""Overlap two archive pipelines with two fixed-template CLEAN workers."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('results/calibrated-fixed16'))
    p.add_argument('--plan',type=Path,default=Path('results/calibrated-fixed16/screen/plan.json'))
    p.add_argument('--source',type=Path,default=Path('data/working/lwa-15slot-flat.ms'))
    p.add_argument('--monitor-only',action='store_true',help='Image verified exports produced by separately running archive jobs')
    p.add_argument('--stokes',nargs='+',choices=('I','V'),default=['I','V'])
    p.add_argument('--reference-root',type=Path,default=Path('results/calibrated-template'))
    p.add_argument('--comparison-mask',type=Path,default=Path('results/calibrated-template/above10deg.fits'))
    args=p.parse_args();plan=json.loads(args.plan.read_text())
    args.root.mkdir(parents=True,exist_ok=True)
    # Start the many-small-block control early so it does not become a late tail.
    archive_plan=args.root/'archive-plan.json'
    archive_plan.write_text(json.dumps(sorted(plan,key=lambda e:e['config']['baseline_group']!=1),indent=2)+'\n')
    env={**os.environ,'PYTHONPATH':'src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
    base=[sys.executable,'scripts/run_error_limits.py','--root',str(args.root),'--plan',str(args.plan),
          '--source',str(args.source),'--reference-root',str(args.reference_root),
          '--comparison-mask',str(args.comparison_mask),'--stokes',*args.stokes]
    with (args.root/'driver-archives.log').open('a') as log:
        archives=(None if args.monitor_only else subprocess.Popen(
            [*base,'--plan',str(archive_plan),'--archive-only','--jobs','2'],env=env,stdout=log,stderr=subprocess.STDOUT))
        submitted=set();futures={};failed=[]
        def image(name,repeat):
            with (args.root/(name+('-repeat-driver.log' if repeat else '-image-driver.log'))).open('a') as output:
                command=[*base,'--cases',name,'--repeat-only' if repeat else '--image-only','--jobs','1']
                return subprocess.run(command,env=env,stdout=output,stderr=subprocess.STDOUT).returncode
        with ThreadPoolExecutor(max_workers=2) as pool:
            while True:
                for case in plan:
                    name=case['name'];report=args.root/(name+'-export.json')
                    if name in submitted or not report.exists():continue
                    try:ready=json.loads(report.read_text())['status']=='passed'
                    except (json.JSONDecodeError,KeyError):continue
                    if ready:
                        for repeat in (False,True):futures[pool.submit(image,name,repeat)]=(name,repeat)
                        submitted.add(name)
                        print(json.dumps(dict(case=name,event='two_images_queued')),flush=True)
                for future in list(futures):
                    if not future.done():continue
                    name,repeat=futures.pop(future);code=future.result()
                    print(json.dumps(dict(case=name,repeat=repeat,returncode=code)),flush=True)
                    if code:failed.append((name,repeat))
                if not futures and ((archives is None and len(submitted)==len(plan)) or
                                    (archives is not None and archives.poll() is not None)):break
                time.sleep(2)
        missing=[e['name'] for e in plan if e['name'] not in submitted]
        code=archives.returncode if archives is not None else None
        result=dict(status='passed' if not failed and not missing and code in (None,0) else 'failed',
                    archive_returncode=code,external_archive_jobs=args.monitor_only,
                    failed_images=failed,missing_exports=missing)
        (args.root/'study-status.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result),flush=True)
        if result['status']!='passed':raise SystemExit(1)


if __name__=='__main__':main()
