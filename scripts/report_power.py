"""Collect median/power-law controls, reconstruction checks and paper figures."""
import argparse
from collections import Counter
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from casacore.tables import table

from report_midrise import save
from run_wsclean import snapshot
from varquanms.container import Reader
from varquanms.power import reconstruct
from varquanms.visibility import quality

ROOT=Path(__file__).resolve().parents[1]
plt.rcParams['svg.hashsalt']='varquanms-power-v1'


def figures(report):
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    codes=np.arange(16);colors=['#2564a0','#c17a20','#338960']
    for row,(gamma,color) in enumerate(zip([1.,.5,.25],colors)):
        levels=reconstruct(codes,4,1.,0.,gamma)
        axes[0].hlines(-row,-1,1,color='#cccccc')
        axes[0].scatter(levels,np.full(16,-row),color=color,s=30,label=f'γ = {gamma:g}')
        axes[0].text(-1,-row+.15,f'γ = {gamma:g}',color=color)
        axes[1].plot(codes,levels,marker='o',ms=3,color=color,label=f'γ = {gamma:g}')
    axes[0].axvline(0,color='gray',ls=':')
    axes[0].set(ylim=(-2.4,.6),yticks=[],xlabel='Residual / stored range t',title='4-bit reconstruction levels around the median')
    axes[1].set(xlabel='Stored unsigned code q',ylabel='(Reconstruction − median) / t',title='Inverse signed power transform')
    axes[1].legend();save(fig,'power-levels')
    records=report['experiments']
    names=['symmetric','zero-uniform','median-uniform','zero-power','median-power-1','symmetric-comparable']
    labels=['Symmetric','Zero\nuniform','Median\nuniform','Zero\npower','Median\npower','Symmetric\ncomparable']
    chosen=[next(r for r in records if r['name']==n) for n in names]
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),constrained_layout=True)
    colors=['#4b6480','#6b91b0','#419b9a','#b6924a','#d66b2a','#756499']
    axes[0].bar(np.arange(6),[r['encoding']['archive_bytes']/1e6 for r in chosen],color=colors)
    axes[0].set(xticks=np.arange(6),xticklabels=labels,ylabel='Whole archive (MB)',title='Metadata and exact exceptions included')
    for label,r,color in zip(labels,chosen,colors):
        axes[1].scatter(r['component_mse'],r['encoding']['archive_bytes']/1e6,color=color,s=65,label=label.replace('\n',' '))
    axes[1].set(xlabel='Measured component MSE',ylabel='Whole archive (MB)',title='Measured rate–distortion controls')
    axes[1].ticklabel_format(axis='x',style='sci',scilimits=(0,0));axes[1].legend(fontsize=8)
    save(fig,'power-comparison')


def summarize(entry):
    e=entry['encoding'];q=e['visibility_quality']['DATA']['by_product']
    hist,sample,gamma,accounting=Counter(),Counter(),Counter(),Counter()
    for p in e['visibility_block_statistics']['DATA']['products']:
        hist.update(p['component_bit_width_histogram']);sample.update(p['complex_sample_bit_width_histogram'])
        gamma.update(p.get('power_gamma_block_histogram',{}));accounting.update(p['logical_bytes_before_secondary_compression'])
    with Reader(e['archive'],verify=False) as reader:
        environment=reader.manifest['environment']
        spec=next(t for t in reader.manifest['tables'] if t['path']=='.')['columns']['DATA']
        arrays=spec['quality_arrays']
        cross=reader.h5[arrays['cross_product_error_covariance']][:]
        adjacent=reader.h5[arrays['adjacent_channel_error_covariance']][:]
        covariance={'cross_product_real':cross.real.tolist(),'cross_product_imag':cross.imag.tolist(),
                    'maximum_absolute_adjacent_channel_by_product':np.max(np.abs(adjacent),axis=0).tolist()}
    count=sum(v['valid_samples'] for v in q.values())
    return {**entry,'environment':environment,'component_mse':sum(v['component_mse']*v['valid_samples'] for v in q.values())/count,
            'block_histogram':dict(hist),'sample_histogram':dict(sample),'gamma_histogram':dict(gamma),
            'logical_byte_accounting':dict(accounting),'error_covariance':covariance}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT/'results/power-20261023')
    parser.add_argument('--figures-only',action='store_true')
    args=parser.parse_args();destination=ROOT/'reports/power-evaluation.json'
    if args.figures_only:figures(json.loads(destination.read_text()));return
    initial=json.loads((args.root/'benchmark/results.json').read_text())
    final=json.loads((args.root/'final-benchmark/results.json').read_text())
    fidelity=json.loads((args.root/'fidelity-benchmark/results.json').read_text())
    entries=initial['experiments']+final['experiments']+fidelity['experiments']
    if any(e['verification_status']!='passed' for e in entries):raise ValueError('Every archive must be verified')
    report={'schema':'varquanms.power-study/1','date_utc':'2026-10-23','source':initial['source'],
            'methodology':initial['notes'],'cache_methodology':initial['cache_methodology'],
            'experiments':[summarize(e) for e in entries],
            'primary_sources':['https://arxiv.org/html/1609.02019v1','https://www.audiolabs-erlangen.de/resources/MIR/FMP/C2/C2S2_DigitalSignalQuantization.html'],
            'formula':{'center':'componentwise median, stored float32 before residual estimation',
                       'encode_coordinate':'M/2 * (1 + sign(x-C)*abs((x-C)/t)**gamma); nearest reconstructed bracketing level',
                       'decode':'u=(2*q-M)/M; x_hat=float32(C+t*sign(u)*abs(u)**(1/gamma)); fixed binary64 operations',
                       'gammas':[1.,.5,.25],'header_bytes':28}}
    selected=[e for e in final['experiments'] if e['name'].startswith('median-power-fast-')]
    if len(selected)!=3:raise ValueError('Require three final timing repetitions')
    hashes=[]
    for e in selected+[next(e for e in initial['experiments'] if e['name']=='median-power-1')]:
        with Reader(e['encoding']['archive'],verify=False) as reader:hashes.append({k:v['sha256'] for k,v in reader.entries.items()})
    if any(h!=hashes[0] for h in hashes[1:]):raise ValueError('Acceleration/repetition changed stored datasets')
    times=[e['process_timings_seconds']['encode_process'] for e in selected]
    report['final_timing']={'process_seconds':times,'median_seconds':statistics.median(times),
                            'minimum_seconds':min(times),'maximum_seconds':max(times),'scientific_payloads_unchanged':True,
                            'scope':'Fresh process, startup+inspection+complete encode+integrity publication. Excludes verification/export/imaging. Twenty-four workers; shared host; no cold-cache claim.'}
    report['export_verification']=json.loads((args.root/'export-verification.json').read_text())
    source=Path(initial['source'])
    with table(str(source),readonly=True,lockoptions='usernoread',ack=False) as t:
        original=t.getcol('DATA');flags=t.getcol('FLAG')|t.getcol('FLAG_ROW')[:,None,None]
        auto=t.getcol('ANTENNA1')==t.getcol('ANTENNA2')
    with table(str(args.root/'restored.ms'),readonly=True,lockoptions='usernoread',ack=False) as t:decoded=t.getcol('DATA')
    names=list(selected[0]['encoding']['visibility_quality']['DATA']['by_product'])
    report['baseline_classes']={}
    for name,mask in (('autocorrelation',auto),('cross_correlation',~auto)):
        report['baseline_classes'][name]=quality(original[mask],decoded[mask],flags[mask],names)[0]
    before=json.loads((args.root/'source-before.json').read_text())
    if before!=snapshot(source):raise ValueError('Source files changed')
    report.update(source_files_unchanged=True,source_file_sha256=before)
    report['imaging']=json.loads((ROOT/'reports/wsclean-power.json').read_text())
    run=json.loads((args.root/'wsclean/complete.json').read_text())
    report['imaging_run']={k:run[k] for k in ('status','container_id','wsclean_version','parameters','elapsed_seconds','source_unchanged')}
    report['limitations']='One integration, absolute engineering bounds, no validated noise model; Stokes I only, CLEAN iteration cap. Calibration, weak polarization, solar peak, narrow features and temporal behavior not accepted by this test. Equal configured bounds are not equal measured distortion.'
    destination.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');figures(report)


if __name__=='__main__':main()
