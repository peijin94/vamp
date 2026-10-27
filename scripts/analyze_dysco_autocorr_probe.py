"""Audit the autocorrelations in the existing two-integration Dysco AF probe.

This is a diagnostic, not a complete-MS archive or imaging comparison.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from varquanms.inspection import _open


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source', type=Path)
    p.add_argument('probe', type=Path)
    p.add_argument('--report', type=Path, required=True)
    args = p.parse_args()
    if args.report.exists():
        raise ValueError('Report must be new')
    names = ('TIME', 'TIME_CENTROID', 'ANTENNA1', 'ANTENNA2', 'UVW', 'FLAG',
             'FLAG_ROW', 'WEIGHT', 'SIGMA', 'FIELD_ID', 'DATA_DESC_ID',
             'OBSERVATION_ID', 'SCAN_NUMBER')
    totals = np.zeros((4,3), dtype=np.float64)
    with _open(args.source) as src, _open(args.probe) as probe:
        n = probe.nrows()
        if n != 124256 or src.nrows() <= n:
            raise ValueError('Expected the documented two-integration prefix')
        manager = probe.getdminfo('DATA')
        spec = manager['SPEC']
        if (manager['TYPE'] != 'DyscoStMan' or spec['normalization'] != 'AF' or spec['dataBitCount'] != 8
                or spec['distribution'] != 'TruncatedGaussian' or spec['distributionTruncation'] != 2.5):
            raise ValueError('Expected the Dysco AF 8-bit probe')
        for start in range(0,n,2048):
            count = min(2048,n-start)
            for name in names:
                a, b = np.asarray(src.getcol(name,start,count)), np.asarray(probe.getcol(name,start,count))
                if a.dtype != b.dtype or a.shape != b.shape or a.tobytes() != b.tobytes():
                    raise ValueError(f'Probe prefix metadata differs: {name}, row {start}')
            auto = src.getcol('ANTENNA1',start,count) == src.getcol('ANTENNA2',start,count)
            if not auto.any():
                continue
            flags = src.getcol('FLAG',start,count)[auto] | src.getcol('FLAG_ROW',start,count)[auto,None,None]
            original = src.getcol('CORRECTED_DATA',start,count)[auto]
            decoded = probe.getcol('DATA',start,count)[auto]
            for k in range(4):
                valid = ~flags[:,:,k] & np.isfinite(original[:,:,k])
                x = original[:,:,k][valid].astype('c16')
                y = decoded[:,:,k][valid].astype('c16')
                if not np.isfinite(y).all():
                    raise ValueError('Unflagged reconstructed autocorrelation is not finite')
                totals[k] += (len(x),np.sum(np.abs(x)**2),np.sum(np.abs(y-x)**2))
    products = {}
    for name, (count,power,error) in zip(('XX','XY','YX','YY'),totals):
        products[name] = {'unflagged_samples':int(count),
                          'relative_complex_rms_error':float(np.sqrt(error/power))}
    result = {'schema':'varquanms.dysco-autocorr-probe/1','status':'exploratory',
              'source':str(args.source),'probe':str(args.probe),'rows':n,'integrations':2,
              'verified_prefix_columns':list(names),'flagged_input':'zero control',
              'normalization':'AF','bits_per_component':8,'distribution':'TruncatedGaussian',
              'distribution_truncation':2.5,'products':products,
              'scope':'First two integrations only; not an archive acceptance or imaging test.'}
    args.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()
