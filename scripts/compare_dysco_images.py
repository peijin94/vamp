"""Compare Dysco DATA imaging with a verified calibrated CORRECTED_DATA reference."""
import argparse
import json
from pathlib import Path

from vamp.imaging import compare_run, file_sha256


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('reference',type=Path)
    p.add_argument('candidate',type=Path)
    p.add_argument('--input-identity',type=Path,required=True)
    p.add_argument('--verification',type=Path,required=True)
    p.add_argument('--mask',type=Path,required=True)
    p.add_argument('--report',type=Path,required=True)
    args=p.parse_args()
    if args.report.exists():
        raise ValueError('Comparison report must be new')
    identity=json.loads(args.input_identity.read_text())
    verification=json.loads(args.verification.read_text())
    ref=json.loads((args.reference/'complete.json').read_text())
    candidate=json.loads((args.candidate/'complete.json').read_text())
    if identity.get('status')!='passed' or verification.get('status')!='passed':
        raise ValueError('Matching input and Dysco verification are required')
    if (ref['parameters'].get('data_column') != 'CORRECTED_DATA'
            or candidate['parameters'].get('data_column') != 'DATA'):
        raise ValueError('Expected calibrated CORRECTED_DATA and verified Dysco DATA imaging')
    if file_sha256(identity['verification']) != identity['verification_sha256']:
        raise ValueError('Calibrated input verification has changed')
    if (Path(verification['source']).resolve()!=Path(identity['prepared_ms']).resolve()
            or Path(verification['decoded']).resolve()!=Path(candidate['source_ms']).resolve()
            or ref['source_sha256']!=identity['original_sha256']
            or candidate['source_sha256']!=verification['decoded_sha256']):
        raise ValueError('Imaged sources do not match the verified calibrated input and Dysco output')
    result=compare_run(args.reference,args.candidate,args.mask,equivalent_calibrated_columns=True)
    final='Cleaning finished because the final threshold was reached.' in (args.candidate/'wsclean.log').read_text()
    result.update(final_threshold_reached=final,input_identity=str(args.input_identity),
                  input_identity_sha256=file_sha256(args.input_identity),
                  dysco_verification=str(args.verification),dysco_verification_sha256=file_sha256(args.verification),
                  matching_scope='DATA contains the verified Dysco-compressed calibrated visibility stream; the reference images original CORRECTED_DATA. All other imaging parameters and the container match.')
    if not final:
        result['clean_warning']='The configured final CLEAN threshold was not reached; image differences are measured but convergence was not established.'
    args.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    n=result['residual_normalized']['image-image.fits']
    print(json.dumps({'status':result['status'],'difference_rms':n['difference_rms'],
                      'difference_over_reference_residual_rms':n['difference_over_reference_residual_rms']}),flush=True)


if __name__=='__main__':
    main()
