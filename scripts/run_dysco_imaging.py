"""Image verified Dysco cases in Stokes I and V with the fixed CLEAN template."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from run_dysco_benchmark import case_name


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--normalization',choices=('AF','RF'),required=True)
    p.add_argument('--bits',nargs='+',type=int,required=True)
    p.add_argument('--truncation',type=float,choices=(1.5,2.5,3.5),default=2.5)
    p.add_argument('--flagged-input',choices=('zero','nan'),default='zero')
    p.add_argument('--stokes',nargs='+',choices=('I','V'),default=['I','V'])
    p.add_argument('--repeat',action='store_true')
    p.add_argument('--jobs',type=int,choices=(1,2),default=2)
    args=p.parse_args()
    root=args.root.resolve()
    identity=root/'input-identity.json'
    if json.loads(identity.read_text()).get('status')!='passed':
        raise ValueError('Input identity was not verified')
    env={**os.environ,'PYTHONPATH':'src:scripts','OPENBLAS_NUM_THREADS':'1',
         'OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    jobs=[]
    for bits in args.bits:
        case=case_name(args.normalization,bits,args.truncation,args.flagged_input)
        archive=root/(case+'.ms')
        verification=root/(case+'-verification-v2.json')
        info=json.loads(verification.read_text())
        if info.get('status')!='passed' or Path(info['decoded']).resolve()!=archive:
            raise ValueError(f'{case}: verified Dysco output required')
        for pol in args.stokes:
            stem=case+'-'+pol.lower()+('-repeat' if args.repeat else '')
            images=root/(stem+'-images')
            comparison=root/(stem+'-comparison.json')
            if images.exists() or comparison.exists():
                raise ValueError(f'Existing or incomplete imaging output: {stem}')
            jobs.append((stem,archive,verification,pol,images,comparison))
    def run(job):
        stem,archive,verification,pol,images,comparison=job
        ref=Path('results/calibrated-template')/('reference' if pol=='I' else 'reference-v')
        image=[sys.executable,'scripts/run_wsclean.py',str(archive),'--output',str(images),
               '--data-column','DATA','--clean-profile','calibrated','--threads','8',
               '--polarizations',pol,'--beam-from',str(ref/'image-image.fits')]
        compare=[sys.executable,'scripts/compare_dysco_images.py',str(ref),str(images),
                 '--input-identity',str(identity),'--verification',str(verification),
                 '--mask','results/calibrated-template/above10deg.fits','--report',str(comparison)]
        started=time.perf_counter()
        with (root/(stem+'-driver.log')).open('x') as log:
            for command in (image,compare):
                result=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
                if result.returncode:
                    raise RuntimeError(f'{stem} failed: see {log.name}')
        print(json.dumps({'case':stem,'status':'complete','wall_seconds':time.perf_counter()-started}),flush=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(run,jobs))


if __name__=='__main__':
    main()
