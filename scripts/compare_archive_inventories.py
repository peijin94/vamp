"""Verify complete partition payload/quality identity across encoder implementations."""
import argparse
import json
from pathlib import Path

from vamp.partitioned import _load, _signature, partition_reader


def compare(baseline,candidate):
    _,a=_load(baseline);_,b=_load(candidate)
    ac=dict(a['config']);bc=dict(b['config'])
    ac.pop('workers',None);bc.pop('workers',None)
    if ac!=bc:raise ValueError('Scientific configurations differ')
    if a['row_count']!=b['row_count'] or len(a['partitions'])!=len(b['partitions']):
        raise ValueError('Partition domains differ')
    fields=('path','shape','dtype','sha256','logical_bytes','role')
    count=logical=0
    for pa,pb in zip(a['partitions'],b['partitions']):
        if (pa['start'],pa['stop'])!=(pb['start'],pb['stop']):raise ValueError('Row mapping differs')
        with partition_reader(baseline,pa) as ra,partition_reader(candidate,pb) as rb:
            if _signature(ra.manifest)!=_signature(rb.manifest):
                raise ValueError('Logical schemas or subtable contents differ')
            ia={e['path']:{k:e[k] for k in fields} for e in ra.manifest['inventory']}
            ib={e['path']:{k:e[k] for k in fields} for e in rb.manifest['inventory']}
            if ia!=ib:
                changed=sorted(k for k in ia.keys()|ib.keys() if ia.get(k)!=ib.get(k))
                raise ValueError('Scientific dataset inventory differs: '+str(changed))
            for ta,tb in zip(ra.manifest['tables'],rb.manifest['tables']):
                for name,sa in ta['columns'].items():
                    sb=tb['columns'][name]
                    for quality in ('quality_summary','quality'):
                        if sa.get(quality)!=sb.get(quality):raise ValueError('Stored quality differs: '+name)
            count+=len(ia);logical+=sum(e['logical_bytes'] for e in ia.values())
    return dict(status='passed',baseline=str(baseline),candidate=str(candidate),
                partitions=len(a['partitions']),dataset_references=count,logical_bytes_compared=logical,
                all_inventory_checksums_verified=True,scientific_dataset_hashes_identical=True,
                logical_schemas_subtables_and_quality_identical=True,
                scope='All scientific dataset shapes/dtypes/hashes, logical schemas, subtable values and quality summaries; execution timings, storage-manager provenance and physical HDF5 layout may differ')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('baseline',type=Path);parser.add_argument('candidate',type=Path)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    if args.report.exists():raise ValueError('Report must be new')
    result=compare(args.baseline,args.candidate)
    args.report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__':main()
