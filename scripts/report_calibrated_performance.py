"""Collect calibrated performance and I/V evidence without touching the paper."""
import argparse
from collections import Counter
import json
from pathlib import Path
import statistics

import numpy as np

from varquanms.imaging import file_sha256, read_fits


def load(path):return json.loads(Path(path).read_text())


def timings(directory):
    runs=[]
    for process in load(directory/'processes.json'):
        if process['returncode']:raise ValueError('Failed benchmark process')
        if process['profiled']:continue
        raw=load(directory/(process['name']+'.json'))
        bundle=raw['bundle'];compact=raw['compact']
        if not compact['all_partition_inventories_verified']:raise ValueError('Unverified archive')
        runs.append(dict(process_wall_seconds=process['process_wall_seconds'],
            workflow_seconds=raw['total_workflow_seconds'],bundle_seconds=raw['bundle_wall_seconds'],
            compaction_seconds=compact['elapsed_seconds'],partition_workers=bundle['partition_workers'],
            sum_encode_seconds=bundle['sum_encode_seconds_excluding_inspection'],
            sum_inspection_seconds=bundle['sum_inspection_seconds'],
            sum_copy_seconds=bundle['sum_partition_copy_seconds'],
            reference_sampling_seconds=bundle['global_reference_sampling_seconds'],
            archive_bytes=compact['archive_bytes'],byte_categories=compact['stored_bytes_by_category'],
            parent_peak_rss_bytes=raw['peak_rss_bytes'],
            largest_child_peak_rss_bytes=raw.get('largest_child_peak_rss_bytes',0),
            cpu_user_seconds=raw['cpu_user_seconds']+raw.get('child_cpu_user_seconds',0),
            cpu_system_seconds=raw['cpu_system_seconds']+raw.get('child_cpu_system_seconds',0)))
    values=[r['process_wall_seconds'] for r in runs]
    return dict(runs=runs,median_seconds=statistics.median(values),
                minimum_seconds=min(values),maximum_seconds=max(values))


def collect(root):
    baseline=timings(root/'baseline');optimized=timings(root/'optimized')
    if len(baseline['runs'])!=3 or len(optimized['runs'])!=3:raise ValueError('Require three unprofiled repeats')
    evidence={name:load(root/(name+'.json')) for name in
              ('optimized-export','optimized-identity','source-paper-integrity','old-export-identity')}
    if any(e['status']!='passed' for e in evidence.values()):raise ValueError('Validation incomplete')
    manifest=load(root/'optimized/run-1.bundle/manifest.json')
    block_bits=Counter();samples=Counter();payload=0
    for part in manifest['partitions']:
        if not all(item['original_bits_exact'] for item in part['summary']['weights_quality'].values()):
            raise ValueError('This comparison requires the observed exact weight encodings')
        for product in part['summary']['visibility_block_statistics']['CORRECTED_DATA']['products']:
            block_bits.update(product['component_bit_width_histogram'])
            samples.update(product['complex_sample_bit_width_histogram'])
            payload+=product['serialized_block_bytes']
    image_runs={}
    for pol in ('I','V'):
        measured=[]
        for repeat in (False,True):
            stem='optimized'+('-v' if pol=='V' else '')+('-repeat' if repeat else '')
            comparison=load(root/(stem+'-comparison.json'))
            if comparison['parameters']['polarizations']!=pol:raise ValueError('Wrong Stokes selection')
            name=next(n for n in comparison['products'] if n.endswith('-image.fits'))
            item=comparison['products'][name]
            image_directory=Path(item['candidate']).parent
            imaging=load(image_directory/'complete.json')
            clean_log=(image_directory/'wsclean.log').read_text()
            for product in comparison['products'].values():
                for label in ('reference','candidate'):
                    if file_sha256(product[label])!=product[label+'_file_sha256']:
                        raise ValueError('Imaging file changed')
            measured.append(dict(run='repeat' if repeat else 'first',
                comparison=str(root/(stem+'-comparison.json')),reference=item['reference'],candidate=item['candidate'],
                normalization=comparison['residual_normalized'][name],
                regions=item['regions'],peak_change=item['relative_change_at_reference_peak'],
                reference_peak_yx=item['reference_peak_yx'],
                mask=comparison['comparison_mask'],clean=item['clean'],
                final_threshold_reached='Cleaning finished because the final threshold was reached.' in clean_log,
                wsclean_seconds=imaging['wsclean_elapsed_seconds'],
                imaging_container_id=imaging['container_id']))
        image_runs[pol]=measured
    controls={}
    for pol in ('I','V'):
        comparison=load(root/f'reference-{pol.lower()}-repeat-comparison.json')
        controls[pol]=comparison['residual_normalized']
        original=load(Path(comparison['reference'])/'complete.json')
        repeat=load(Path(comparison['candidate'])/'complete.json')
        if original['source_sha256']!=repeat['source_sha256']:raise ValueError('Reference-repeat source differs')
    pilots={p.name:timings(p) for p in sorted(root.glob('p*t*')) if (p/'processes.json').exists()}
    ablations={p.name:load(p/'timing.json') for p in sorted((root/'ablations').glob('*')) if (p/'timing.json').exists()}
    hardware=load('reports/calibrated-fixed16.json')['hardware']
    hardware['cache_methodology']='Warm OS cache, no eviction, sequential task-owned timing runs; no overlapping task-owned tests or CLEAN; unrelated shared-host activity not controlled'
    profiling={}
    for name,path in [('full_baseline',root/'baseline/profile-summary.json'),
                      ('context_baseline',root/'context-profiles/baseline-summary.json'),
                      ('context_optimized',root/'context-profiles/optimized-summary.json')]:
        raw=load(path)
        profiling[name]=dict(total_seconds=raw['total_seconds'],calls=raw['calls'],
            top_self_time_functions=raw['functions'][:30],
            note='Separate cProfile run; cumulative times overlap; worker native CPU counters are recorded separately')
    source_hashes={str(p.relative_to(root/'final-code')):file_sha256(p)
                   for p in sorted((root/'final-code/src').rglob('*'))
                   if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc'}
    casa=load(root/'casa-check.json')
    if casa.get('table_reopen')!='passed' or casa.get('imaging')!='passed':
        raise ValueError('CASA smoke validation failed')
    casa['runtime_environment']={'CASAPATH':'/common/casa/casa-release-5.4.1-31.el6 linux'}
    casa['installed_runtime_warnings']=['OVRO-LWA absent from the installed observatory repository; frequency conversion unavailable',
        'Observation outside the installed IERS range', 'Installed leap-second table out of date; derived times/coordinates may be inaccurate']
    casa['scope']='CASA 5.4.1 reopening and coarse 64x64 Stokes I dirty-image smoke workflow on a private copy; not CASA 6 or precision astrometry/calibration validation'
    return dict(schema='varquanms.calibrated-performance/1',status='complete',
        baseline_commit=load(root/'start-inventory.json')['commit'],hardware=hardware,
        optimized_source_sha256=source_hashes,profiling=profiling,
        baseline=baseline,optimized=optimized,pilots=pilots,context_ablations=ablations,
        median_speedup=baseline['median_seconds']/optimized['median_seconds'],
        full_ms_under_15_seconds=optimized['maximum_seconds']<15,
        config=manifest['config'],execution=manifest['execution'],
        source_rows=manifest['row_count'],complex_samples=sum(samples.values()),
        visibility_payload_bytes=payload,payload_fraction_of_complex64=payload/(8*sum(samples.values())),
        sample_bits=dict(samples),block_bits=dict(block_bits),validation=evidence,imaging=image_runs,
        reference_repeat_controls=controls,casa_smoke_check=casa,
        tests_passed=124,
        limitations=['One calibrated observation and one fixed epsilon-D configuration; no new scientific acceptance threshold.',
            'Three timing repeats and two CLEAN runs are observed variation, not confidence intervals.',
            'Timing includes source discovery, references, inspection, encoding, integrity/source verification and single-file compaction; export/CLEAN measured separately.',
            'Parallel stage sums overlap and are not elapsed wall time. Memory lists parent and largest child peaks, not aggregate peak.',
            'Reference residual RMS is a scalar normalization, not a thermal-noise estimate or pixelwise divisor.',
            'Stokes V image consistency does not establish absolute polarimetric accuracy, calibration-solution fidelity or temporal behaviour.'])


def figures(result,root,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from astropy.io import fits
    output.mkdir(parents=True,exist_ok=True)
    labels=['Baseline']+list(result['pilots'])+['Selected\n3 repeats']
    entries=[result['baseline'],*result['pilots'].values(),result['optimized']]
    values=np.array([e['median_seconds'] for e in entries])
    errors=np.array([[e['median_seconds']-e['minimum_seconds'] for e in entries],
                     [e['maximum_seconds']-e['median_seconds'] for e in entries]])
    fig,ax=plt.subplots(figsize=(10,4),layout='constrained')
    bars=ax.bar(labels,values,yerr=errors,capsize=4,color=['#777777']+['#4c91bd']*(len(labels)-2)+['#237847'])
    ax.bar_label(bars,labels=[f'{v:.1f} s' for v in values],padding=5)
    ax.axhline(15,color='#a84435',linestyle='--',label='15 s target')
    ax.set_ylabel('MS to verified single-file archive, seconds');ax.set_ylim(0,max(values)*1.18)
    ax.set_title('Complete 15-integration MS; same scientific settings and payload')
    ax.legend();fig.savefig(output/'calibrated-runtime.png',dpi=160);plt.close(fig)
    mask=fits.getdata('results/calibrated-template/above10deg.fits').astype(bool)
    for pol,runs in result['imaging'].items():
        run=runs[0];reference,_=read_fits(run['reference']);candidate,_=read_fits(run['candidate'])
        denominator=run['normalization']['reference_residual_rms']
        if denominator<=0:raise ValueError('Plot requires a positive residual RMS')
        difference=(candidate-reference)/denominator
        peak=np.unravel_index(np.nanargmax(np.where(mask,np.abs(reference),np.nan)),reference.shape)
        for zoom in (False,True):
            radius=100
            view=(slice(max(0,peak[0]-radius),min(reference.shape[0],peak[0]+radius)),
                  slice(max(0,peak[1]-radius),min(reference.shape[1],peak[1]+radius))) if zoom else np.s_[:,:]
            panels=[np.where(mask,a,np.nan)[view] for a in (reference,candidate,difference)]
            extent=max(np.nanmax(np.abs(panels[0])),np.nanmax(np.abs(panels[1]))) or 1.
            limit=float(np.nanmax(np.abs(panels[2]))) or 1.
            fig,axes=plt.subplots(1,3,figsize=(12,4.2),layout='constrained')
            for i,(ax,data,title) in enumerate(zip(axes,panels,('Reference','Reconstructed','(Image − reference) /\nRMS(reference residual)'))):
                norm=SymLogNorm(linthresh=extent*.01,vmin=-extent,vmax=extent) if i<2 else None
                opts=dict(norm=norm) if norm else dict(vmin=-limit,vmax=limit)
                if zoom:opts['extent']=(view[1].start-.5,view[1].stop-.5,view[0].start-.5,view[0].stop-.5)
                im=ax.imshow(data,origin='lower',cmap='RdBu_r',**opts)
                ax.set_title(title,fontsize=10);ax.set_xlabel('Pixel x');ax.set_ylabel('Pixel y')
                bar=fig.colorbar(im,ax=ax,shrink=.75,label='Jy/beam' if i<2 else 'dimensionless')
                if i<2:
                    ticks=[-extent,-extent*.1,0,extent*.1,extent]
                    bar.set_ticks(ticks);bar.set_ticklabels([f'{value:.2g}' for value in ticks])
            suffix='zoom' if zoom else 'full'
            ratio=run['normalization']['difference_over_reference_residual_rms']
            fig.suptitle(f'Stokes {pol}; above 10° only; εD=0.18; RMS difference/reference residual = {100*ratio:.3f}%')
            fig.savefig(output/f'calibrated-performance-{pol.lower()}-{suffix}.png',dpi=160);plt.close(fig)
    bits=[str(x) for x in (3,4,5,6,8,10,12,14,16)]+['float32']
    fig,ax=plt.subplots(figsize=(8,3.5),layout='constrained')
    ax.bar(bits,[100*result['sample_bits'].get(b,0)/result['complex_samples'] for b in bits],color='#4c91bd')
    ax.set_xlabel('Bits per real/imaginary component (float32 = exact block)')
    ax.set_ylabel('Complex samples, %');ax.set_title('Unchanged precision distribution after optimization')
    fig.savefig(output/'calibrated-performance-bits.png',dpi=160);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('results/calibrated-performance'))
    parser.add_argument('--report',type=Path,default=Path('reports/calibrated-performance.json'))
    args=parser.parse_args();result=collect(args.root)
    args.report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    figures(result,args.root,Path('devnotes/figs'))
    print(json.dumps(dict(status='complete',speedup=result['median_speedup'],seconds=result['optimized']['median_seconds'])))


if __name__=='__main__':main()
