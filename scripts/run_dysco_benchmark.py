"""Run one-pass Dysco MS encodes sequentially in the pinned LINC container."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time


def case_name(normalization, bits, truncation, flagged_input='zero'):
    suffix = '' if truncation == 2.5 else f'-t{round(10*truncation):02d}'
    suffix += '-nan' if flagged_input == 'nan' else ''
    return f'dysco-{normalization.lower()}-{bits:02d}{suffix}'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--normalization', choices=('AF', 'RF'), required=True)
    p.add_argument('--bits', type=int, nargs='+', required=True)
    p.add_argument('--truncation', type=float, choices=(1.5,2.5,3.5), default=2.5)
    p.add_argument('--flagged-input', choices=('zero','nan'), default='zero')
    p.add_argument('--cpus', type=float, default=None, help='Optional Docker CPU quota')
    p.add_argument('--image', default='astronrd/linc:5.2.1')
    args = p.parse_args()
    source, root = args.input.resolve(), args.output.resolve()
    script = Path(__file__).with_name('dysco_encode.py').resolve()
    if not source.is_dir() or root.is_relative_to(source):
        raise ValueError('Input must be an MS and output must be outside it')
    if (any(b not in (2,3,4,6,8,10,12,16) for b in args.bits)
            or len(set(args.bits)) != len(args.bits) or (args.cpus is not None and args.cpus <= 0)):
        raise ValueError('Unsupported or duplicate bit rate')
    root.mkdir(parents=True, exist_ok=True)
    image = json.loads(subprocess.check_output(['docker', 'image', 'inspect', args.image], text=True))[0]
    for bits in args.bits:
        name = case_name(args.normalization, bits, args.truncation, args.flagged_input)
        report, log, process = (root / (name + suffix) for suffix in ('-encode.json','-encode.log','-process.json'))
        if any(path.exists() for path in (report, log, process, root/(name+'.ms'), root/(name+'.ms.partial'))):
            raise ValueError(f'Existing or incomplete output: {name}')
        command = ['docker','run','--rm','--network','none']
        if args.cpus is not None:
            command += ['--cpus',str(args.cpus)]
        command += ['--user',f'{os.getuid()}:{os.getgid()}',
                   '-e','OPENBLAS_NUM_THREADS=1','-e','OMP_NUM_THREADS=1',
                   '-v',f'{source}:/input.ms:ro','-v',f'{root}:/work',
                   '-v',f'{script}:/encode.py:ro','-w','/work',image['Id'],
                   'python3','-u','/encode.py','--input','/input.ms',
                   '--output',f'/work/{name}.ms','--report',f'/work/{name}-encode.json',
                   '--bits',str(bits),'--normalization',args.normalization,
                   '--truncation',str(args.truncation),'--flagged-input',args.flagged_input]
        started = time.perf_counter()
        with log.open('x') as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
        record = {'schema':'vamp.dysco-process/1','case':name,'command':command,
                  'container_id':image['Id'],'container_repo_digests':image.get('RepoDigests',[]),
                  'process_wall_seconds':time.perf_counter()-started,'returncode':result.returncode,
                  'timing_scope':'Docker launch through process exit; source preparation and subsequent verification/imaging are excluded',
                  'cache_policy':'No cache eviction; warm/shared-host conditions; sequential encodes in this invocation',
                  'cpu_quota':args.cpus}
        process.write_text(json.dumps(record,indent=2)+'\n')
        print(json.dumps(record), flush=True)
        if result.returncode or not report.exists() or json.loads(report.read_text()).get('status') != 'complete':
            raise RuntimeError(f'{name} failed; inspect {log} and the partial output')


if __name__ == '__main__':
    main()
