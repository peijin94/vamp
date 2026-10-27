"""Time one fixed source context with precomputed full-MS references for ablation.

This is not the end-to-end MS conversion benchmark. Run under PYTHONPATH for
the desired frozen implementation, with no concurrent tests or CLEAN jobs.
"""
import argparse
import json
from pathlib import Path
import resource
import time

from vamp.config import Config
from vamp.partitioned import _encode_context


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--config',type=Path,required=True)
    parser.add_argument('--reference-bundle',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--index',type=int,default=0)
    args=parser.parse_args()
    manifest=json.loads((args.reference_bundle/'manifest.json').read_text())
    config=Config.load(args.config)
    if config.to_dict()!=manifest['config']:raise ValueError('Ablation scientific configuration differs')
    part=manifest['partitions'][args.index]
    args.output.mkdir()
    task=(args.source.resolve(),args.output.resolve(),config,manifest.get('distribution_references'),
          args.index,part['start'],part['stop'])
    start=time.perf_counter()
    record,_=_encode_context(task)
    result=dict(context=record,wall_seconds=time.perf_counter()-start,
                peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                cpu_user_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime,
                cpu_system_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_stime,
                scope='One source context, copy/view + inspection + encode + verification; shared reference sampling and whole-archive compaction excluded')
    (args.output/'timing.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
