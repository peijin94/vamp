"""Run isolated-in-sequence fresh-process workflow benchmarks and a separate profile."""
import argparse
import importlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time


def child(args):
    Config=importlib.import_module(args.package+'.config').Config
    encode_partitioned=importlib.import_module(args.package+'.partitioned').encode_partitioned
    compact_bundle=importlib.import_module(args.package+'.shared').compact_bundle
    cfg=Config.load(args.config)
    before=time.perf_counter()
    bundle=encode_partitioned(args.source,args.output.with_suffix('.bundle'),cfg,
                              partition_workers=args.partition_workers)
    middle=time.perf_counter()
    compact=compact_bundle(args.output.with_suffix('.bundle'),args.output.with_suffix('.h5'))
    record=dict(bundle=bundle,compact=compact,
                bundle_wall_seconds=middle-before,total_workflow_seconds=time.perf_counter()-before,
                peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                cpu_user_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime,
                cpu_system_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_stime,
                child_cpu_user_seconds=resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime,
                child_cpu_system_seconds=resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime,
                largest_child_peak_rss_bytes=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*1024,
                memory_scope='Parent and largest terminated child RSS peaks; simultaneous aggregate process-tree peak is not measured')
    args.output.with_suffix('.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--code-root',type=Path,default=Path('.'))
    p.add_argument('--partition-workers',type=int,default=1)
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--profile-first',action='store_true')
    p.add_argument('--child',action='store_true')
    p.add_argument('--package',choices=('vamp','varquanms'),default='vamp',help=argparse.SUPPRESS)
    args=p.parse_args()
    if args.child:return child(args)
    package=next((name for name in ('vamp','varquanms')
                  if (args.code_root/'src'/name/'__init__.py').is_file()),None)
    if package is None:
        raise ValueError('code-root must contain a vamp or legacy source tree')
    if args.repeats<0 or args.repeats==0 and not args.profile_first:
        raise ValueError('Select positive repeats or a separate profile')
    args.output.mkdir(parents=True,exist_ok=True)
    env={**os.environ,'PYTHONPATH':str(args.code_root.resolve()/'src'),
         'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
    schedule=([('profile',True)] if args.profile_first else [])+[(f'run-{i+1}',False) for i in range(args.repeats)]
    records=[]
    for name,profile in schedule:
        output=args.output/name
        if output.with_suffix('.json').exists():raise ValueError('Benchmark output already exists')
        command=[sys.executable]
        if profile:command+=['-m','cProfile','-o',str(output.with_suffix('.prof'))]
        command += [str(Path(__file__).resolve()),'--child','--source',str(args.source.resolve()),
                    '--config',str(args.config.resolve()),'--output',str(output.resolve()),
                    '--partition-workers',str(args.partition_workers),'--package',package]
        start=time.perf_counter()
        with output.with_suffix('.log').open('w') as log:
            completed=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
        record=dict(name=name,profiled=profile,command=command,returncode=completed.returncode,
                    process_wall_seconds=time.perf_counter()-start)
        records.append(record)
        (args.output/'processes.json').write_text(json.dumps(records,indent=2)+'\n')
        print(json.dumps(record),flush=True)
        if completed.returncode:raise SystemExit(completed.returncode)


if __name__=='__main__':main()
