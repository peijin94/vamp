"""Streaming visibility error diagnostics, with autocorrelations kept separate."""
import argparse
import json
from pathlib import Path

import numpy as np

from varquanms.inspection import _open, CORRELATIONS


def analyze(source, decoded, source_column='CORRECTED_DATA', decoded_column='CORRECTED_DATA'):
    amplitude_edges = np.r_[0., np.logspace(-2,8,21), np.inf]
    length_edges = np.array([0.,30.,100.,300.,1000.,3000.,np.inf])
    moments = {}; amplitude = {}; lengths = {}; channel = {}
    def add(name, a, e, valid, baseline):
        x, err = a[valid].astype('c16'), e[valid]
        n = len(x)
        if not n: return
        power = (x.real*x.real+x.imag*x.imag)/2
        squared = (err.real*err.real+err.imag*err.imag)/2
        nonzero = power > 0
        row = moments.setdefault(name,np.zeros(9,dtype='f8'))
        row += [n,power.sum(),squared.sum(),err.real.sum(),err.imag.sum(),
                np.sum(squared[nonzero]/power[nonzero]),np.count_nonzero(nonzero),
                np.count_nonzero(~nonzero),np.count_nonzero((~nonzero)&(squared > 0))]
        for dest, bins, coordinate in [(amplitude,amplitude_edges,np.abs(x)),
                                        (lengths,length_edges,np.broadcast_to(baseline[:,None],a.shape)[valid])]:
            accum = dest.setdefault(name,np.zeros((len(bins)-1,3),dtype='f8'))
            index = np.searchsorted(bins,coordinate,side='right')-1
            for c,values in enumerate((np.ones(n),power,squared)):
                accum[:,c] += np.bincount(index,weights=values,minlength=len(bins)-1)
        accum=channel.setdefault(name,np.zeros((a.shape[1],3),dtype='f8'))
        ee=np.where(valid,np.abs(e)**2/2,0.)
        pp=np.where(valid,np.abs(a.astype('c16'))**2/2,0.)
        accum += np.column_stack([valid.sum(axis=0),pp.sum(axis=0),ee.sum(axis=0)])
    with _open(source) as src, _open(decoded) as dst:
        if src.nrows()!=dst.nrows(): raise ValueError('Row counts differ')
        ids=np.unique(src.getcol('DATA_DESC_ID'))
        if len(ids)!=1: raise ValueError('Study requires a single data description')
        with _open(source/'DATA_DESCRIPTION') as t: polid=t.getcell('POLARIZATION_ID',int(ids[0]))
        with _open(source/'POLARIZATION') as t: names=[CORRELATIONS[int(v)] for v in t.getcell('CORR_TYPE',int(polid))]
        if set(names)!= {'XX','XY','YX','YY'}: raise ValueError('Study expects four linear-feed products')
        xx,yy=names.index('XX'),names.index('YY')
        for start in range(0,src.nrows(),2048):
            count=min(2048,src.nrows()-start)
            a=src.getcol(source_column,start,count)
            b=dst.getcol(decoded_column,start,count)
            f=src.getcol('FLAG',start,count)|src.getcol('FLAG_ROW',start,count)[:,None,None]
            valid=~f & np.isfinite(a)
            if not np.all(np.isfinite(b[valid])): raise ValueError('Nonfinite reconstruction')
            error=np.zeros(a.shape,dtype='c16')
            error[valid]=b[valid].astype('c16')-a[valid].astype('c16')
            auto=src.getcol('ANTENNA1',start,count)==src.getcol('ANTENNA2',start,count)
            baseline=np.linalg.norm(src.getcol('UVW',start,count),axis=1)
            for p,name in enumerate(names):
                for kind,rows in [('cross',~auto),('auto',auto)]:
                    add(kind+'/'+name,a[:,:,p],error[:,:,p],valid[:,:,p]&rows[:,None],baseline)
            paired=valid[:,:,xx]&valid[:,:,yy]&~auto[:,None]
            stokes=(a[:,:,xx].astype('c16')+a[:,:,yy].astype('c16'))/2
            add('cross/I-paired',stokes,(error[:,:,xx]+error[:,:,yy])/2,paired,baseline)
    def finish(row):
        n,power,squared,re,im,relative,nz,z,ze=row
        return dict(samples=int(n),signal_component_power=power/n,component_mse=squared/n,
                    relative_rmse=float(np.sqrt(squared/power)) if power else None,
                    bias_real=re/n,bias_imag=im/n,
                    mean_sample_relative_mse=relative/nz if nz else None,
                    zero_amplitude_samples=int(z),zero_amplitude_changed_samples=int(ze))
    def bin_records(raw):
        return [dict(samples=int(n),signal_component_power=p/n if n else None,
                     component_mse=e/n if n else None,
                     relative_rmse=float(np.sqrt(e/p)) if p else None) for n,p,e in raw]
    return dict(schema='varquanms.visibility-error-diagnostics/1',source=str(source),decoded=str(decoded),
        source_column=source_column,decoded_column=decoded_column,
        products={k:finish(v) for k,v in moments.items()},
        amplitude_edges=[float(v) if np.isfinite(v) else None for v in amplitude_edges],
        length_edges_m=[float(v) if np.isfinite(v) else None for v in length_edges],
        by_amplitude={k:bin_records(v) for k,v in amplitude.items()},
        by_baseline_length={k:bin_records(v) for k,v in lengths.items()},
        by_channel={k:bin_records(v) for k,v in channel.items()},
        scope='Unflagged finite samples; I=(XX+YY)/2 only where both are valid. Cross/auto separated. No imager weight approximation. Relative sample metric excludes zero amplitude and reports zero changes separately.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source',type=Path);p.add_argument('decoded',type=Path);p.add_argument('--report',type=Path,required=True)
    p.add_argument('--source-column',default='CORRECTED_DATA')
    p.add_argument('--decoded-column',default='CORRECTED_DATA')
    args=p.parse_args()
    if args.report.exists():raise ValueError('Report must be new')
    args.report.write_text(json.dumps(analyze(args.source,args.decoded,args.source_column,args.decoded_column),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
