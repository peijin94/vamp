"""Compare autocorrelations on exactly the same first two integrations."""
import argparse
from contextlib import ExitStack
import json
from pathlib import Path

import numpy as np
from vamp.inspection import _open, CORRELATIONS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--candidate',nargs=3,action='append',metavar=('LABEL','MS','COLUMN'),required=True)
    parser.add_argument('--rows',type=int,required=True)
    parser.add_argument('--report',type=Path,required=True)
    args=parser.parse_args()
    if args.report.exists() or args.rows < 1 or len({x[0] for x in args.candidate}) != len(args.candidate):
        raise ValueError('Use a new report, positive row count and unique candidate labels')
    with ExitStack() as stack:
        src=stack.enter_context(_open(args.source))
        if args.rows > src.nrows():
            raise ValueError('Row limit exceeds source')
        times=np.unique(src.getcol('TIME',0,args.rows))
        if (len(times)!=2 or (args.rows < src.nrows()
                and src.getcell('TIME',args.rows) == src.getcell('TIME',args.rows-1))):
            raise ValueError('The prefix must contain exactly two complete integrations')
        data_desc=np.unique(src.getcol('DATA_DESC_ID',0,args.rows))
        if len(data_desc)!=1:
            raise ValueError('Expected a single data description')
        with _open(args.source/'DATA_DESCRIPTION') as table:
            polid=table.getcell('POLARIZATION_ID',int(data_desc[0]))
        with _open(args.source/'POLARIZATION') as table:
            names=[CORRELATIONS[int(x)] for x in table.getcell('CORR_TYPE',int(polid))]
        if set(names) != {'XX','XY','YX','YY'}:
            raise ValueError('Expected the four documented linear products')
        tables=[]
        for label,path,column in args.candidate:
            t=stack.enter_context(_open(Path(path)))
            if t.nrows() < args.rows or column not in t.colnames():
                raise ValueError(f'{label}: missing requested prefix or visibility column')
            tables.append((label,path,column,t))
        totals={label:np.zeros((len(names),3),dtype=np.float64) for label,_,_,_ in tables}
        keys=('TIME','ANTENNA1','ANTENNA2','UVW','FLAG','FLAG_ROW','WEIGHT','SIGMA','DATA_DESC_ID','FIELD_ID')
        for start in range(0,args.rows,2048):
            count=min(2048,args.rows-start)
            expected={key:np.asarray(src.getcol(key,start,count)) for key in keys}
            auto=expected['ANTENNA1']==expected['ANTENNA2']
            flags=expected['FLAG'][auto] | expected['FLAG_ROW'][auto,None,None]
            original=src.getcol('CORRECTED_DATA',start,count)[auto]
            for label,path,column,t in tables:
                for key in keys:
                    actual=np.asarray(t.getcol(key,start,count)); reference=expected[key]
                    if actual.dtype!=reference.dtype or actual.shape!=reference.shape or actual.tobytes()!=reference.tobytes():
                        raise ValueError(f'{label}: prefix {key} differs at row {start}')
                decoded=t.getcol(column,start,count)[auto]
                for k in range(len(names)):
                    valid=~flags[:,:,k] & np.isfinite(original[:,:,k])
                    x=original[:,:,k][valid].astype('c16'); y=decoded[:,:,k][valid].astype('c16')
                    if not np.isfinite(y).all():
                        raise ValueError(f'{label}: unflagged nonfinite sample')
                    totals[label][k] += (len(x),np.sum(np.abs(x)**2),np.sum(np.abs(y-x)**2))
    cases=[]
    for label,path,column,_ in tables:
        products={name:{'unflagged_samples':int(count),'relative_complex_rms_error':float(np.sqrt(err/power))}
                  for name,(count,power,err) in zip(names,totals[label])}
        cases.append({'label':label,'ms':path,'column':column,'products':products})
    result={'schema':'vamp.autocorrelation-prefix-comparison/1','status':'measured',
            'source':str(args.source),'rows':args.rows,'integrations':2,'verified_prefix_columns':list(keys),
            'time_values':times.tolist(),
            'scope':'Unflagged autocorrelations in the same input row prefix; no imaging or full-MS acceptance is implied.',
            'cases':cases}
    args.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'cases':len(cases),'rows':args.rows}),flush=True)


if __name__ == '__main__':
    main()
