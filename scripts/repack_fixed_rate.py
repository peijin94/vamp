"""Remove streaming-chunk padding without changing any benchmark codec bytes."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

from vamp.container import Reader
from vamp.imaging import file_sha256
from vamp.shared import compact_bundle, load_shared


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path('results/calibrated-fixed25'))
    p.add_argument('--plan', type=Path)
    p.add_argument('--cases', nargs='*')
    p.add_argument('--jobs', type=int, choices=(1,2), default=1)
    args = p.parse_args()
    plan=json.loads((args.plan or args.root/'screen/plan.json').read_text())
    if args.cases:
        if set(args.cases)-{e['name'] for e in plan}:raise ValueError('Unknown case')
        plan=[e for e in plan if e['name'] in args.cases]
    def repack(case):
        name = case['name']
        output = args.root/(name+'-tight.h5')
        report = args.root/(name+'-tight-compact.json')
        if report.exists():
            return
        result = compact_bundle(args.root/(name+'.bundle'), output)
        current = load_shared(output)
        previous = load_shared(args.root/(name+'.h5'))
        checked = 0
        for old, new in zip(previous['partitions'], current['partitions'], strict=True):
            with Reader(args.root/(name+'.h5'), group=old['group']) as a, Reader(output, group=new['group']) as b:
                assert a.h5['manifest'][...].tobytes() == b.h5['manifest'][...].tobytes()
                assert a.entries == b.entries
                checked += len(a.entries)
        payload=sum(product['serialized_block_bytes'] for part in current['partitions']
                    for product in part['summary']['visibility_block_statistics']['CORRECTED_DATA']['products'])
        physical = result['stored_bytes_by_category']['visibility_block_payload']
        assert physical <= payload
        if 'payload_fraction_target' in case:
            rate = json.loads((args.root/(name+'-rate.json')).read_text())
            assert payload == rate['payload_bytes'] <= case['payload_fraction_target'] * rate['original_complex64_bytes']
            result['physical_payload_budget_passed'] = True
        result['matches_imaged_archive'] = dict(
            previous_archive=str(args.root/(name+'.h5')),
            previous_archive_sha256=file_sha256(args.root/(name+'.h5')),
            all_scientific_manifests_bitwise_equal=True,
            verified_array_references=checked,
            all_logical_array_checksums_equal=True,
            previous_archive_bytes=(args.root/(name+'.h5')).stat().st_size)
        result['archive_sha256'] = file_sha256(output)
        report.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
        print(name, result['archive_bytes'], physical, flush=True)
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        list(pool.map(repack, plan))


if __name__ == '__main__':
    main()
