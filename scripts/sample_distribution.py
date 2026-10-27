"""Compare sampled global reference stability without modifying the MS."""
import argparse
from dataclasses import replace
import json
from pathlib import Path

from vamp.config import Config
from vamp.references import sample_references


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('output',type=Path)
    args=p.parse_args();args.output.mkdir(exist_ok=False)
    base=Config.load('configs/calibrated-normalized.json')
    records=[]
    for fraction,seed in ((.02,20261024),(.03,20261024),(.05,20261024),(.02,20261025)):
        name=f'fraction{int(fraction*100)}-seed{seed}'
        cfg=replace(base,normalized_bounds=dict(base.normalized_bounds,sampling_fraction=fraction,sampling_seed=seed))
        record=sample_references(args.source,cfg)
        (args.output/(name+'.json')).write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
        print(json.dumps(dict(name=name,seconds=record['elapsed_seconds'],
            configurations=len(record['groups']),strata=sum(len(classes) for g in record['groups'].values()
                    for column in g['columns'].values() for classes in column.values()))),flush=True)
        records.append(dict(name=name,references=record))
    (args.output/'stability.json').write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
