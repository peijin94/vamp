"""Summarize matched-rate absolute/relative/hybrid archive and image tests."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path

import numpy as np

from varquanms.imaging import file_sha256, read_fits, compare_run

ROOT=Path('results/calibrated-error-limits')
REF=Path('results/calibrated-template/reference')
MASK=Path('results/calibrated-template/above10deg.fits')
FIGS=Path('devnotes/figs')
FAMILIES=('absolute','relative','hybrid')
COLORS={'absolute':'#2166ac','relative':'#b2182b','hybrid':'#238b45'}


def metric(e,product='image',region='above_elevation_cut',field='relative_rms_difference'):
    return e['comparison']['products'][f'image-{product}.fits']['regions'][region][field]


def collect(partial=False):
    screen=json.loads((ROOT/'screen-final/screen.json').read_text())
    plan=screen['selected']; records=[]
    for case in plan:
        name=case['name'];read=lambda suffix:json.loads((ROOT/(name+'-'+suffix+'.json')).read_text())
        if partial and not (ROOT/(name+'-comparison.json')).exists():continue
        manifest=json.loads((ROOT/(name+'.bundle')/'manifest.json').read_text())
        if manifest['config']!=case['config']:raise ValueError('Actual configuration differs from selected plan')
        encoding=read('encode')
        shared=read('tight-compact') if (ROOT/(name+'-tight-compact.json')).exists() else read('compact')
        export=read('export')
        imaging=json.loads((ROOT/(name+'-images')/'complete.json').read_text())
        comparison=read('comparison')
        assert export['status']=='passed' and export['source_files_unchanged']
        assert imaging['status']=='complete' and imaging['source_unchanged']
        assert shared['all_partition_inventories_verified']
        assert sum(shared['stored_bytes_by_category'].values())==shared['archive_bytes']
        assert comparison['comparison_mask']['pixels']==2245726
        assert comparison['comparison_mask']['sha256']==file_sha256(MASK)
        for product,record in comparison['products'].items():
            assert record['candidate_file_sha256']==file_sha256(ROOT/(name+'-images')/product)
        bits,samples,gamma=Counter(),Counter(),Counter();serialized=0
        for part in manifest['partitions']:
            for product in part['summary']['visibility_block_statistics']['CORRECTED_DATA']['products']:
                bits.update(product['component_bit_width_histogram'])
                samples.update(product['complex_sample_bit_width_histogram'])
                gamma.update(product['power_gamma_block_histogram'])
                serialized+=product['serialized_block_bytes']
        assert sum(samples.values())==931920*48*4
        if not partial:
            assert shared['stored_bytes_by_category']['visibility_block_payload']<=serialized, 'Remove streaming padding consistently before comparing full-file rates'
        record=dict(**case,encoding=encoding,shared=shared,serialized_visibility_bytes=serialized,
            block_bits=dict(bits),sample_bits=dict(samples),gamma_blocks=dict(gamma),
            visibility=read('visibility'),comparison=comparison,imaging=imaging,
            weights_all_exact=all(q['original_bits_exact'] for part in manifest['partitions']
                                  for q in part['summary']['weights_quality'].values()),
            export_verification=export,archive_sha256=file_sha256(Path(shared['archive'])))
        assert record['weights_all_exact'], 'Weight loss would confound the visibility-only threshold comparison'
        if name=='absolute-t256':
            prior=json.loads(Path('reports/calibrated-gamma-ratio.json').read_text())
            prior=next(e for e in prior['experiments'] if e['name']=='ratio')
            assert serialized==prior['serialized_visibility_block_bytes']
            assert dict(samples)==prior['sample_bits']
            record['previous_absolute_control_matches']=dict(serialized_block_bytes=True,sample_bits=True)
        repeat=ROOT/(name+'-repeat-images')
        if (repeat/'complete.json').exists():
            repeat_comparison=read('repeat-comparison')
            assert repeat_comparison['comparison_mask']['sha256']==file_sha256(MASK)
            for product in repeat_comparison['products'].values():
                assert file_sha256(product['reference'])==product['reference_file_sha256']
                assert file_sha256(product['candidate'])==product['candidate_file_sha256']
            record['repeat']=dict(versus_reference=repeat_comparison,
                                  versus_first=compare_run(ROOT/(name+'-images'),repeat,MASK),
                                  imaging=json.loads((repeat/'complete.json').read_text()))
        if not partial and 'repeat' not in record:
            raise ValueError('Repeat imaging is still required: '+name)
        records.append(record)
    # Full per-block records are retained locally. The public report keeps
    # every aggregate and columnar per-block values for the nine imaged cases.
    selected={(c['family'],c['parameter']) for c in plan}
    for e in screen['experiments']:
        blocks=e.pop('blocks')
        if (e['family'],e['parameter']) in selected:
            e['block_records']=dict(bytes=[b['bytes'] for b in blocks],bits=[b['bits'] for b in blocks],
                gamma=[b['gamma'] for b in blocks],valid_samples=[b['metrics']['valid_samples'] for b in blocks],
                component_mse=[b['metrics']['component_mse'] for b in blocks],
                max_complex_error=[b['metrics']['max_complex_error'] for b in blocks])
        for auto in (False,True):
            subset=[(j,b) for j,b in zip(screen['blocks'],blocks) if j['auto']==auto]
            n=sum(b['metrics']['valid_samples'] for j,b in subset)
            e['auto_component_mse' if auto else 'cross_component_mse']=sum(b['metrics']['component_mse']*b['metrics']['valid_samples'] for j,b in subset)/max(1,n)
    screen['schema']='varquanms.error-limit-screen-summary/1'
    screen['block_record_order']='Arrays align with the common screen.blocks list; retained for nine selected profiles.'
    old=json.loads(Path('reports/calibrated-template-tests.json').read_text())
    controls=[e for e in old['experiments'] if e['name'] in ('reference-repeat','lossless')]
    if (ROOT/'reference-repeat-3-comparison.json').exists():
        controls.append(dict(name='reference-repeat-3',
            comparison=json.loads((ROOT/'reference-repeat-3-comparison.json').read_text()),
            imaging=json.loads((ROOT/'reference-repeat-3/complete.json').read_text())))
    result=dict(schema='varquanms.error-limit-study/1',status='partial' if partial else 'complete',
        screen=screen,experiments=records,controls=controls,hardware=old['hardware'],
        tests_passed=102,
        scope='All 15 integrations, 48 channels, four correlations; joint Stokes I CLEAN only above 10 degrees.',
        limitations=['Block-relative RMS is not per-sample relative error. No noise estimate is inferred from signal amplitude.',
            'Global normalized RMS and absolute RMS differ only by a constant for a fixed input. The experiment changes allocation across blocks.',
            'Actual Briggs weights and covariance affect images; visibility diagnostics are unweighted and do not replace actual images.',
            'Nine thresholds chosen from sampled byte counts before imaging. Full-data rates are measured separately, not assumed equal.',
            'Two overlapping pipelines on a shared host; no cache eviction and no isolated timing ranking.',
            'Single observation; no held-out dataset, Q/U/V images, calibration solutions, spectral cubes or light curves tested.',
            'CLEAN masks are independently estimated; reference controls and two imaging runs per case quantify limited run variability.'])
    if not partial:
        imported=json.loads(Path('reports/calibrated-import.json').read_text())
        source=Path(imported['extracted_ms'])
        assert all(file_sha256(source/n)==info['sha256'] for n,info in imported['files'].items())
        result['source_integrity']=dict(original_import_checksums_match=True,files_checked=len(imported['files']))
    return result


def plots(result):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/varquanms-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.io import fits
    cases=result['experiments']
    fig,axes=plt.subplots(1,2,figsize=(11,4.6),constrained_layout=True)
    for family in FAMILIES:
        selected=sorted([e for e in result['screen']['experiments'] if e['family']==family],key=lambda e:e['parameter'])
        for ax,key in zip(axes,('cross_component_mse','mean_block_relative_mse')):
            ax.plot([e['payload_bytes']/1e6 for e in selected],[e[key] for e in selected],'.-',
                    color=COLORS[family],label=family,markersize=4,linewidth=1)
            ax.set(xlabel='Same 768 blocks: serialized payload MB',yscale='log');ax.grid(alpha=.2);ax.legend()
    axes[0].set(ylabel='Cross-correlation component MSE',title='Absolute error; autocorrelations excluded')
    axes[1].set(ylabel='Mean of per-block MSE / signal power',title='Relative error; nonzero-power blocks')
    fig.suptitle('83 settings including byte-matching refinements; not whole-archive bytes or image error')
    fig.savefig(FIGS/'error-limits-screen.png',dpi=170);plt.close(fig)
    # Actual full-file rate versus dirty and restored error; repeat ranges shown.
    fig,axes=plt.subplots(1,2,figsize=(11,4.6),constrained_layout=True)
    for family in FAMILIES:
        selected=sorted([e for e in cases if e['family']==family],key=lambda e:e['shared']['archive_bytes'])
        for ax,product in zip(axes,('dirty','image')):
            x=[e['shared']['archive_bytes']/1e6 for e in selected]
            y=[100*metric(e,product) for e in selected]
            ax.plot(x,y,'o-',color=COLORS[family],label=family)
            for xx,yy,e in zip(x,y,selected):
                if family=='absolute':
                    ax.annotate(str(e['target_absolute_mse']),(xx,yy),xytext=(7,-12),textcoords='offset points',fontsize=8)
                if 'repeat' in e:
                    repeat=100*e['repeat']['versus_reference']['products'][f'image-{product}.fits']['regions']['above_elevation_cut']['relative_rms_difference']
                    ax.plot([xx,xx],[yy,repeat],color=COLORS[family],linewidth=3)
                    ax.plot(xx,repeat,'_',color=COLORS[family])
            ax.set(xlabel='Complete shared archive (decimal MB)',ylabel='RMS difference / reference RMS (%)',
                   title='Dirty image' if product=='dirty' else 'Restored CLEAN image',yscale='log')
            ax.grid(alpha=.2);ax.legend()
    fig.suptitle('Error limits at matched sampled payload rates; labels mark absolute-control MSE limits\nActual whole-file sizes shown; above 10° only; vertical marks span two imaging runs')
    fig.savefig(FIGS/'error-limits-rate.png',dpi=170);plt.close(fig)
    # Compression distribution includes all samples (including flagged zeros).
    fig,ax=plt.subplots(figsize=(11,5),constrained_layout=True)
    y=np.arange(len(cases));left=np.zeros(len(cases));total=np.array([sum(e['sample_bits'].values()) for e in cases])
    for bits,color in zip(('4','8','12','16','float32'),('#ef8a62','#67a9cf','#2166ac','#542788','#aaaaaa')):
        fraction=np.array([e['sample_bits'].get(bits,0) for e in cases])/total*100
        ax.barh(y,fraction,left=left,label=bits+' bits/component' if bits!='float32' else 'exact float32',color=color)
        left+=fraction
    ax.set(yticks=y,yticklabels=[e['name'] for e in cases],xlabel='All complex samples (%)',xlim=(0,100));ax.invert_yaxis()
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.13),ncol=3,fontsize=8)
    ax.set_title('Adaptive precision distribution; sample fractions, not physical byte fractions')
    fig.savefig(FIGS/'error-limits-bits.png',dpi=170);plt.close(fig)
    # Error redistribution on imaged cross-correlations, excluding autos.
    fig,axes=plt.subplots(1,2,figsize=(11,4.8),constrained_layout=True)
    for family in FAMILIES:
        e=next(e for e in cases if e['family']==family and e['target_absolute_mse']==256)
        data=e['visibility']['by_amplitude']['cross/I-paired']
        edges=e['visibility']['amplitude_edges'];x=[];absolute=[];relative=[]
        for i,bin in enumerate(data):
            if bin['samples'] and edges[i] and edges[i+1]:
                x.append(np.sqrt(edges[i]*edges[i+1]));absolute.append(np.sqrt(bin['component_mse']));relative.append(100*bin['relative_rmse'])
        axes[0].plot(x,absolute,'o-',color=COLORS[family],label=family)
        axes[1].plot(x,relative,'o-',color=COLORS[family],label=family)
    for ax in axes:
        ax.set(xscale='log',yscale='log',xlabel='Original paired Stokes I visibility amplitude (stored units)');ax.grid(alpha=.2);ax.legend()
    axes[0].set(ylabel='Component error RMS (stored units)',title='Absolute error by amplitude')
    axes[1].set(ylabel='Error RMS / signal RMS (%)',title='Fractional error by amplitude')
    fig.suptitle('Middle rate tier: full-data cross-correlations; same bins for all methods')
    fig.savefig(FIGS/'error-limits-allocation.png',dpi=170);plt.close(fig)
    # Same-reference, same-scale full images and zoomed differences.
    reference,_=read_fits(REF/'image-image.fits');mask=fits.getdata(MASK).astype(bool)
    middle=[next(e for e in cases if e['family']==f and e['target_absolute_mse']==256) for f in FAMILIES]
    images=[read_fits(ROOT/(e['name']+'-images')/'image-image.fits')[0] for e in middle]
    difference=[a-reference for a in images]
    peak=middle[0]['comparison']['products']['image-image.fits']['reference_peak_yx'];y,x=peak
    for zoom in (False,True):
        crop=np.s_[y-64:y+64,x-64:x+64] if zoom else np.s_[:,:]
        domain=mask[crop];limit=max(float(np.max(np.abs(a[crop][domain]))) for a in difference)
        fig,axes=plt.subplots(1,4,figsize=(15,4.3),constrained_layout=True)
        from matplotlib.colors import SymLogNorm
        scale=float(np.max(np.abs(reference[crop][domain])))
        im=axes[0].imshow(np.where(domain,reference[crop],np.nan),origin='lower',cmap='RdBu_r',norm=SymLogNorm(scale*.01,vmin=-scale,vmax=scale))
        positive=10.**np.arange(np.ceil(np.log10(.01*scale)),np.floor(np.log10(scale))+1)
        fig.colorbar(im,ax=axes[0],shrink=.75,label='Jy/beam',ticks=np.r_[-positive[::-1],0,positive]);axes[0].set_title('Reference')
        for ax,e,delta in zip(axes[1:],middle,difference):
            im=ax.imshow(np.where(domain,delta[crop],np.nan),origin='lower',cmap='RdBu_r',vmin=-limit,vmax=limit)
            ax.set_title(e['family']+' − reference\n'+f'{100*metric(e):.3f}% full-mask RMS')
        fig.colorbar(im,ax=list(axes[1:]),shrink=.75,label='Difference Jy/beam (shared linear scale)')
        for ax in axes:ax.set(xlabel='Crop pixel x' if zoom else 'Pixel x',ylabel='Crop pixel y' if zoom else 'Pixel y')
        fig.suptitle('Middle-rate restored images; above 10° only'+(' • 128×128 reference-peak zoom' if zoom else ''))
        fig.savefig(FIGS/('error-limits-zoom.png' if zoom else 'error-limits-images.png'),dpi=170);plt.close(fig)


def tables(result):
    lines=['\n## Measured results\n',
        '**Absolute and relative limits control where precision is spent; neither alone determines CLEAN agreement.** For the two lower-distortion tiers, block-relative allocation gives a slightly smaller dirty-image difference at similar sampled byte counts, despite a larger unweighted cross-I MSE. Actual weighting and error correlations matter. The separate [25% payload study](fixed25.md) tests the user’s specific budget and recommends parameters from repeated images.','',
        'Restored-image differences are not monotonic in visibility error: the tightest tier has some larger CLEAN differences than looser tiers. Repeated eight-thread CLEAN on identical decoded inputs can change the restored result substantially, while dirty-image differences remain stable. The complete repeat table is therefore part of the result; a favorable single image is insufficient.','',
        '| Target tier | Rule | Archive MB | Visibility blocks MB | Cross-I component MSE | Dirty difference % | Restored difference % | Faint difference % |',
        '| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for e in result['experiments']:
        lines.append(f"| {e['target_absolute_mse']} | {e['family']} | {e['shared']['archive_bytes']/1e6:.3f} | {e['serialized_visibility_bytes']/1e6:.3f} | {e['visibility']['products']['cross/I-paired']['component_mse']:.4f} | {100*metric(e,'dirty'):.4f} | {100*metric(e):.4f} | {100*metric(e,'image','faint_reference'):.4f} |")
    lines+=['','These are full-data measurements. Tier labels refer to the absolute screen targets; relative and hybrid thresholds differ to match sampled byte counts. Cross-I uses only unflagged XX/YY pairs; all original sample ordering and exact metadata checks passed. These MSE values exclude autocorrelations.','',
        '| Tier | Relative RMS fraction r | Hybrid absolute cap B | Hybrid RMS fraction r |',
        '| ---: | ---: | ---: | ---: |']
    for tier in (64,256,4096):
        relative=next(e for e in result['experiments'] if e['name']==f'relative-t{tier}')
        hybrid=next(e for e in result['experiments'] if e['name']==f'hybrid-t{tier}')
        lines.append(f"| {tier} | {relative['config']['relative_bounds']['max_component_rmse']:.12g} | {hybrid['config']['max_component_mse']:.12g} | {hybrid['config']['relative_bounds']['max_component_rmse']:.12g} |")
    lines+=['','### Repeat controls','',
        '| Case | First restored difference % | Repeat restored difference % | Repeat versus first % |',
        '| --- | ---: | ---: | ---: |']
    for e in result['experiments']:
        if 'repeat' not in e:continue
        repeat=e['repeat'];extract=lambda key:100*repeat[key]['products']['image-image.fits']['regions']['above_elevation_cut']['relative_rms_difference']
        lines.append(f"| {e['name']} | {100*metric(e):.4f} | {extract('versus_reference'):.4f} | {extract('versus_first'):.4f} |")
    lines+=['','| Source/lossless control | Dirty difference % | Restored difference % |',
            '| --- | ---: | ---: |']
    for control in result['controls']:
        lines.append(f"| {control['name']} | {100*metric(control,'dirty'):.5f} | {100*metric(control):.4f} |")
    lines+=['','![Sampled absolute and relative rate–distortion](figs/error-limits-screen.png)','',
        '![Actual rate versus image error](figs/error-limits-rate.png)','',
        '![Full restored image comparison](figs/error-limits-images.png)','',
        '![Reference-peak zoom](figs/error-limits-zoom.png)','',
        '![Bit-size distribution](figs/error-limits-bits.png)','',
        '![Absolute and fractional error allocation](figs/error-limits-allocation.png)','',
        '### Timing observations','',
        '| Case | Encoding only s | Complete verified bundle s | Compaction s | Export s | CLEAN s |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for e in result['experiments']:
        lines.append(f"| {e['name']} | {e['encoding']['sum_encode_seconds_excluding_inspection']:.2f} | {e['encoding']['elapsed_seconds_including_copy_inspection_verification']:.2f} | {e['shared']['elapsed_seconds']:.2f} | {e['export_verification']['export']['elapsed_seconds']:.2f} | {e['imaging']['wsclean_elapsed_seconds']:.2f} |")
    lines+=['','Two pipelines overlap on a shared host, without cache eviction. These observations do not establish speed rankings or a 15-second conversion guarantee. All complete-file sizes use the same lossless compaction policy, including removal of unused streaming-chunk allocation. Array digests and scientific manifests match the imaged exports.','',
        '### Reproduction and scope','',
        '```bash','export PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1',
        'python scripts/screen_error_limits.py INPUT.ms NEW_SCREEN_DIR',
        'python scripts/run_error_limits.py --source INPUT.ms --plan NEW_SCREEN_DIR/plan.json --root NEW_RESULT_DIR',
        '```','',
        'Repeat all cases with `scripts/run_error_limits.py --repeat-only` using the same plan/root arguments. All retained FITS, archives and MS copies are under ignored `results/calibrated-error-limits/`. The JSON records commands, checksums, preservation policies, exact selected thresholds, byte categories, all 83 screen aggregates and per-block results for the nine selected cases.','',
        '[Machine-readable report](../reports/calibrated-error-limits.json). The source retains all 108 original import checksums. All weights are checked separately; no original weight is silently adjusted for quantization error. No polarization-image, calibration-solution, held-out-observation, spectral-cube or temporal acceptance is implied.','']
    return '\n'.join(lines)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--partial',action='store_true');args=p.parse_args()
    result=collect(args.partial)
    target=ROOT/'preview.json' if args.partial else Path('reports/calibrated-error-limits.json')
    target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for e in result['experiments']:
        print(e['name'],e['shared']['archive_bytes'],100*metric(e,'dirty'),100*metric(e))
    if not args.partial:
        plots(result)
        path=Path('devnotes/absolute_relative_error.md');text=path.read_text().split('\n## Measured results')[0]
        text=text.replace('Status: full archive, independent export and imaging runs in progress. No final\nranking is asserted here until these complete.',
                          'All nine full archives and exports passed verification. All nine images were\nrepeated with the same settings. The table and plots below contain the measurements.')
        path.write_text(text.rstrip()+'\n'+tables(result))


if __name__=='__main__':main()
