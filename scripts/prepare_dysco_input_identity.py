"""Bind a losslessly verified calibrated-only MS to the original reference."""
import argparse
import json
from pathlib import Path

from run_wsclean import snapshot
from vamp.archive import compare_ms
from vamp.imaging import file_sha256
from vamp.partitioned import _load


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('original',type=Path)
    parser.add_argument('prepared',type=Path)
    parser.add_argument('--archive',type=Path,required=True)
    parser.add_argument('--verification',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    if args.report.exists():
        raise ValueError('Input identity report must be new')
    verified=json.loads(args.verification.read_text())
    checks=verified.get('column_checks',[])
    exact_visibility=any(c.get('table')=='.' and c.get('column')=='CORRECTED_DATA'
                         and c.get('export')=='bitwise_exact' for c in checks)
    derived=any(c.get('table')=='.' and c.get('column')=='DATA'
                and c.get('export')=='derived_calibrated_copy_matches_CORRECTED_DATA_bits' for c in checks)
    if (verified.get('status')!='passed' or verified.get('keep_calib_only') is not True
            or verified.get('source_files_unchanged') is not True
            or verified.get('quantized_columns_match_archive_decoder_bitwise')
            or not verified.get('export',{}).get('derived_calibrated_DATA')
            or not exact_visibility or not derived
            or Path(verified['source']).resolve()!=args.original.resolve()
            or Path(verified['restored']).resolve()!=args.prepared.resolve()):
        raise ValueError('A lossless calibrated-only export with a verified DATA alias is required')
    original_hashes=snapshot(args.original)
    if original_hashes != verified['source_file_sha256']:
        raise ValueError('Original source differs from the export verification')
    _,manifest=_load(args.archive)
    if manifest['config']['mode']!='lossless' or not manifest['config'].get('keep_calib_only'):
        raise ValueError('The reference archive must be lossless and calibrated-only')
    prepared_hashes=snapshot(args.prepared)
    fresh_checks=compare_ms(args.original,args.prepared,keep_calib_only=True,calib_config=manifest['config'])
    if snapshot(args.original)!=original_hashes or snapshot(args.prepared)!=prepared_hashes:
        raise ValueError('The original or prepared MS changed during the independent comparison')
    result={'schema':'vamp.dysco-input-identity/1','status':'passed',
            'original_ms':str(args.original.resolve()),'prepared_ms':str(args.prepared.resolve()),
            'verification':str(args.verification),'verification_sha256':file_sha256(args.verification),
            'lossless_archive':str(args.archive.resolve()),'lossless_archive_sha256':file_sha256(args.archive),
            'independent_exact_column_checks':len(fresh_checks),
            'original_sha256':original_hashes,'prepared_sha256':prepared_hashes}
    args.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':'passed','original_files':len(result['original_sha256']),
                      'prepared_files':len(result['prepared_sha256'])}),flush=True)


if __name__=='__main__':
    main()
