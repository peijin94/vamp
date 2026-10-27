"""Decode an ordered bundle independently, then compare the complete original MS."""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from run_wsclean import snapshot
from vamp.archive import compare_ms, read_column
from vamp.cli import _write_new
from vamp.inspection import _open
from vamp.metadata import equal
from vamp.partitioned import decode_partitioned, _load, partition_reader, partition_verification


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle', type=Path); p.add_argument('source', type=Path)
    p.add_argument('restored', type=Path); p.add_argument('--report', type=Path, required=True)
    p.add_argument('--verify-existing', action='store_true', help='verify a previously decoded MS without writing it again')
    data_export = p.add_mutually_exclusive_group()
    data_export.add_argument('--ms-v2-data', action='store_true', default=None)
    data_export.add_argument('--no-ms-v2-data', dest='ms_v2_data', action='store_false', default=None)
    args = p.parse_args()
    if args.verify_existing and args.ms_v2_data is not None:
        raise ValueError('DATA creation options apply only to new exports')
    if args.report.exists(): raise ValueError('Report must be new')
    started = time.perf_counter()
    before = snapshot(args.source)
    # The decoder receives no source path.
    if args.verify_existing:
        if not args.restored.is_dir(): raise ValueError('No existing restored MS')
        _, manifest = _load(args.bundle)
        exported = {'status':'existing_export_verified_below','output':str(args.restored),
                    'elapsed_seconds':None,'note':'Export occurred in an earlier decoder invocation; its wall time was not retained.'}
    else:
        exported = decode_partitioned(args.bundle, args.restored, ms_v2_data=args.ms_v2_data)
        _, manifest = _load(args.bundle)
    selected = manifest['config']['columns'] if manifest['config']['mode'] != 'lossless' else []
    selected = list(selected) + sorted({c for part in manifest['partitions'] for c,q in part['summary'].get('weights_quality', {}).items() if not q['original_bits_exact']})
    checks = compare_ms(args.source, args.restored, approximate_columns=selected,
                        keep_calib_only=manifest['config'].get('keep_calib_only', False), calib_config=manifest['config'])
    reference_check = None
    if manifest['config'].get('normalized_bounds') is not None:
        from vamp.config import Config
        from vamp.references import sample_references
        fresh = sample_references(args.source, Config(**manifest['config']))
        recorded = manifest['distribution_references']
        if any(fresh[k] != recorded[k] for k in ('schema','source_rows','sampling','groups')):
            raise ValueError('Global distribution references differ from source resampling')
        reference_check = dict(status='passed', source_resampling_identical=True,
                               seconds=fresh['elapsed_seconds'])
    with _open(args.restored) as t:
        for part in manifest['partitions']:
            validated = partition_verification(args.bundle, part)
            if validated['status'] != 'passed': raise ValueError('Partition did not pass source verification')
            with partition_reader(args.bundle, part) as reader:
                spec = next(s for s in reader.manifest['tables'] if s['path'] == '.')
                for name in selected:
                    values = read_column(reader, spec['columns'][name], 0, part['stop']-part['start'])
                    for start in range(0, len(values), 1024):
                        wanted = values[start:start+1024]
                        actual = np.asarray(t.getcol(name, startrow=part['start']+start, nrow=len(wanted)))
                        if not equal(actual, wanted): raise ValueError('Export differs from archive decoder')
    if before != snapshot(args.source): raise ValueError('Source files changed during export verification')
    result = {'schema': 'vamp.partition-export-verification/1', 'status': 'passed',
              'source': str(args.source), 'restored': str(args.restored), 'export': exported,
              'column_checks': checks, 'quantized_columns_match_archive_decoder_bitwise': selected,
              'partitions_verified_against_source_at_encode': len(manifest['partitions']),
              'keep_calib_only': manifest['config'].get('keep_calib_only', False),
              'source_files_unchanged': True, 'source_file_sha256': before,
              'distribution_reference_verification': reference_check,
              'elapsed_seconds_including_export_and_verification': time.perf_counter()-started,
              'scope': 'All declared exact metadata and auxiliary columns compared bitwise against the full source; selected columns match the independently decoded archives. CASA and imaging are separate checks.'}
    _write_new(args.report, json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'exact_columns': len(checks), 'export': exported}, indent=2))


if __name__ == '__main__': main()
