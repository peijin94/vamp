"""Collect the max/std gamma test, matched controls and elevation-masked images."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import re

import numpy as np

from vamp.imaging import compare_run, file_sha256, plot_comparison, read_fits

ROOT = Path('results/calibrated-gamma-ratio')
REFERENCE = Path('results/calibrated-template/reference')
MASK = Path('results/calibrated-template/above10deg.fits')
FIGS = Path('devnotes/figs')


def collect():
    screen = json.loads((ROOT/'screen-final/screen.json').read_text())
    # Keep per-block rates/errors and the common sampled-block identities in
    # the public report; the large per-block timing trace stays local.
    screen['schema'] = 'vamp.gamma-ratio-screen-summary/1'
    for experiment in screen['experiments']:
        experiment['blocks'] = [dict(bytes=b['bytes'], bits=b['bits'], gamma=b['gamma'],
                                      samples=b['samples'],valid_samples=b['metrics']['valid_samples'],
                                      component_mse=b['metrics']['component_mse'],
                                      max_complex_error=b['metrics']['max_complex_error'],
                                      exceptions=b['exceptions']) for b in experiment['blocks']]
        experiment['gamma_blocks'].pop('None',None)
    screen['block_record_order'] = 'Each experiment.blocks list aligns exactly with screen.blocks'
    old = json.loads(Path('reports/calibrated-shared-storage.json').read_text())
    baseline = next(e for e in old['experiments'] if e['name']=='amplitude256')
    records=[]
    for label in ('search','ratio','inlier'):
        manifest=json.loads((ROOT/(label+'.bundle')/'manifest.json').read_text())
        encode=json.loads((ROOT/(label+'-encode.json')).read_text())
        compact=json.loads((ROOT/(label+'-compact.json')).read_text())
        export=json.loads((ROOT/(label+'-export.json')).read_text())
        imaging=json.loads((ROOT/(label+'-images')/'complete.json').read_text())
        assert export['status']=='passed' and export['source_files_unchanged']
        comparison=compare_run(REFERENCE,ROOT/(label+'-images'),MASK)
        (ROOT/(label+'-comparison.json')).write_text(json.dumps(comparison,indent=2,allow_nan=False)+'\n')
        gammas,bits,samples,worker=Counter(),Counter(),Counter(),Counter()
        serialized=0; quality={}
        for part in manifest['partitions']:
            s=part['summary']; spec=s['visibility_block_statistics']['CORRECTED_DATA']
            worker.update(spec['worker_timings_seconds'])
            for product in spec['products']:
                gammas.update(product['power_gamma_block_histogram'])
                bits.update(product['component_bit_width_histogram'])
                samples.update(product['complex_sample_bit_width_histogram'])
                serialized+=product['serialized_block_bytes']
            for name,q in s['visibility_quality']['CORRECTED_DATA']['by_product'].items():
                acc=quality.setdefault(name,dict(valid_samples=0,squared_error=0.,max_complex_error=0.,bias_real_sum=0.,bias_imag_sum=0.))
                n=q['valid_samples']; acc['valid_samples']+=n;acc['squared_error']+=q['component_mse']*n
                acc['max_complex_error']=max(acc['max_complex_error'],q['max_complex_error'])
                acc['bias_real_sum']+=q['bias_real']*n;acc['bias_imag_sum']+=q['bias_imag']*n
        products={name:dict(valid_samples=q['valid_samples'],component_mse=q['squared_error']/q['valid_samples'],
                           max_complex_error=q['max_complex_error'],bias_real=q['bias_real_sum']/q['valid_samples'],
                           bias_imag=q['bias_imag_sum']/q['valid_samples']) for name,q in quality.items()}
        count=sum(q['valid_samples'] for q in products.values())
        matches=re.findall(r'Iteration\s+(\d+)',(ROOT/(label+'-images')/'wsclean.log').read_text())
        records.append(dict(name=label,config=manifest['config'],encoding=encode,shared=compact,
            serialized_visibility_block_bytes=serialized,gamma_blocks=dict(gammas),block_bits=dict(bits),sample_bits=dict(samples),
            worker_timings_seconds=dict(worker),visibility_by_product=products,
            component_mse=sum(q['component_mse']*q['valid_samples'] for q in products.values())/count,
            export_verification=export,imaging=imaging,comparison=comparison,
            maximum_logged_iteration=max(int(v) for v in matches) if matches else None,
            shared_file_sha256=file_sha256(compact['archive'])))
    imported=json.loads(Path('reports/calibrated-import.json').read_text())
    source=Path(imported['extracted_ms'])
    assert all(file_sha256(source/name)==info['sha256'] for name,info in imported['files'].items())
    repeated=dict(versus_reference=compare_run(REFERENCE,ROOT/'ratio-repeat-images',MASK),
                  versus_first=compare_run(ROOT/'ratio-images',ROOT/'ratio-repeat-images',MASK),
                  imaging=json.loads((ROOT/'ratio-repeat-images/complete.json').read_text()))
    return dict(schema='vamp.gamma-ratio-experiment/1',screen=screen,experiments=records,ratio_repeat=repeated,
        baseline=dict(name='symmetric-amplitude256',shared_bytes=baseline['archive_bytes'],
                      comparison=baseline['image_comparison'],previous_image_repeat_percent=[.68605007,.54955293]),
        source_integrity=dict(original_import_checksums_match=True,files_checked=len(imported['files'])),
        hardware=old['hardware'],software=old['software'],tests_passed=97,
        decision='Keep both ratio selectors opt-in. This discrete-gamma trial shows reduced native search work, but no practical whole-archive size improvement or demonstrated full-encoding latency improvement; defaults unchanged.',
        limitations=['Only this calibrated dataset and Stokes I; no polarization, calibration-solution or held-out validation.',
                     'Ratio and inlier variants use R0=4 selected before full-data tests; screen also reports R0=3 and 6.',
                     'Three full encodes used 24 workers on a shared host and overlapped; wall times are observations, not isolated speed rankings.',
                     'Gamma selection narrows the candidate set; all stored-value error checks and exact fallback still apply.',
                     'A smaller max/std ratio does not establish absence of RFI or outliers; std is signal spread and depends on sample count.'])


def screen_plots(screen):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/vamp-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    entries={e['name']:e for e in screen['experiments']}
    ratios=np.array([b['statistics']['max_std_ratio'] if b['statistics']['max_std_ratio'] is not None else np.nan for b in screen['blocks']])
    search=entries['power-search3-b256']
    gamma=np.array([b['gamma'] if b['gamma'] is not None else np.nan for b in search['blocks']])
    valid=np.isfinite(ratios)&np.isfinite(gamma)
    bins=np.geomspace(2,40,35)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for product in ('XX','XY','YX','YY'):
        selected=np.array([b['product']==product for b in screen['blocks']])&np.isfinite(ratios)
        axes[0].hist(ratios[selected],bins=bins,histtype='step',label=product,linewidth=1.5)
    for boundary in (4,8):axes[0].axvline(boundary,color='black',alpha=.4,linestyle='--')
    axes[0].set(xscale='log',xlabel='Maximum absolute residual / pooled component std',ylabel='Sampled blocks',title='Residual tail ratio after stored-median subtraction')
    axes[0].legend()
    count=axes[1].hist2d(np.log2(ratios[valid]),-np.log2(gamma[valid]),bins=[np.log2(bins),[-.5,.5,1.5,2.5]],norm=LogNorm(),cmap='Blues')
    axes[1].plot(np.log2([2,4,4,8,8,40]),[0,0,1,1,2,2],color='#d62728',label='Ratio rule, R0=4')
    axes[1].set(xticks=np.log2([2,4,8,16,32]),xticklabels=['2','4','8','16','32'],
                yticks=[0,1,2],yticklabels=['1','1/2','1/4'],xlabel='Maximum / std',ylabel='Selected gamma',title='Full search choices versus the ratio rule')
    axes[1].legend(loc='upper left');fig.colorbar(count[3],ax=axes[1],label='Blocks')
    fig.suptitle('Matched blocks from integrations 0, 7, 14; undefined ratios excluded; MSE ceiling 256')
    fig.savefig(FIGS/'gamma-ratio-selection.png',dpi=170);plt.close(fig)
    labels=['symmetric','median-uniform','power-search2','power-search3',
            'max-std-v1-r3','max-std-inlier-v1-r3','max-std-v1-r4','max-std-inlier-v1-r4','max-std-v1-r6','max-std-inlier-v1-r6']
    titles=['Symmetric','Median, gamma 1','Search 1 / 1/2','Search 1 / 1/2 / 1/4',
            'Ratio R0=3','Inlier R0=3','Ratio R0=4','Inlier R0=4','Ratio R0=6','Inlier R0=6']
    data=[entries[label+'-b256'] for label in labels]
    fig,axes=plt.subplots(1,2,figsize=(11,6),sharey=True,constrained_layout=True)
    y=np.arange(len(data));base=entries['power-search3-b256']['payload_bytes']
    for ax,values,title in [(axes[0],[100*(e['payload_bytes']/base-1) for e in data],'Payload increase over full search (%)'),
                            (axes[1],[e['elapsed_seconds'] for e in data],'Sampled encoding wall time (s), 8 workers')]:
        bars=ax.barh(y,values,color=['#888888']*4+['#1f77b4','#2ca02c']*3)
        ax.bar_label(bars,fmt='%.2f',padding=3,fontsize=8)
        ax.set(yticks=y,yticklabels=titles,xlabel=title,xlim=(0,max(values)*1.2));ax.grid(axis='x',alpha=.2)
    axes[0].invert_yaxis()
    fig.suptitle('Same 768 blocks, bounds, scales and exception policy; budget 256\nPayload includes block headers/exceptions; excludes archive metadata. Shared-host timing.')
    fig.savefig(FIGS/'gamma-ratio-screen.png',dpi=170);plt.close(fig)


def image_plots(result):
    import matplotlib.pyplot as plt
    from astropy.io import fits
    mask=fits.getdata(MASK).astype(bool)
    plot_comparison(REFERENCE/'image-image.fits',ROOT/'ratio-images/image-image.fits',FIGS/'gamma-ratio-images.png',comparison_mask=mask)
    reference,_=read_fits(REFERENCE/'image-image.fits')
    peak=result['experiments'][0]['comparison']['products']['image-image.fits']['reference_peak_yx']
    y,x=peak;crop=np.s_[y-64:y+64,x-64:x+64]
    deltas=[]
    for e in result['experiments']:
        a,_=read_fits(ROOT/(e['name']+'-images')/'image-image.fits');deltas.append((a-reference)[crop])
    limit=max(float(np.max(np.abs(a[mask[crop]]))) for a in deltas)
    fig,axes=plt.subplots(1,3,figsize=(12,4),constrained_layout=True)
    for ax,e,delta in zip(axes,result['experiments'],deltas):
        im=ax.imshow(np.where(mask[crop],delta,np.nan),origin='lower',cmap='RdBu_r',vmin=-limit,vmax=limit,
                     extent=(x-64.5,x+63.5,y-64.5,y+63.5))
        ax.set(title=e['name']+' minus reference',xlabel='Pixel x',ylabel='Pixel y')
    fig.colorbar(im,ax=axes,label='Jy/beam',shrink=.8)
    fig.suptitle('128-pixel zoom at reference peak; elevation >10°; common full linear difference scale')
    fig.savefig(FIGS/'gamma-ratio-zoom.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    y=np.arange(len(result['experiments']))
    for ax,field,keys,colors,title in [
        (axes[0],'gamma_blocks',['1.0','0.5','0.25'],['#888888','#1f77b4','#d62728'],'Gamma: fraction of quantized blocks'),
        (axes[1],'sample_bits',['4','8','12','16','float32'],['#9467bd','#1f77b4','#2ca02c','#ff7f0e','#888888'],'Precision: fraction of all complex samples')]:
        left=np.zeros(len(y))
        totals=np.array([sum(e[field].values()) for e in result['experiments']])
        for key,color in zip(keys,colors):
            values=np.array([e[field].get(key,0) for e in result['experiments']])/totals*100
            ax.barh(y,values,left=left,color=color,label=key);left+=values
        ax.set(xlim=(0,100),yticks=y,yticklabels=[e['name'] for e in result['experiments']],xlabel='Percent',title=title)
        ax.invert_yaxis();ax.legend(loc='upper center',bbox_to_anchor=(.5,-.17),ncol=len(keys),fontsize=9)
    fig.suptitle('Whole dataset; bits apply to each real/imaginary component\nFractions are not byte shares; float32 includes all-flagged zero blocks')
    fig.savefig(FIGS/'gamma-ratio-distributions.png',dpi=170);plt.close(fig)


def write_note(result):
    full={e['name']:e for e in result['experiments']}
    percentage=lambda e:100*e['comparison']['products']['image-image.fits']['regions']['above_elevation_cut']['relative_rms_difference']
    repeat_percentage=100*result['ratio_repeat']['versus_reference']['products']['image-image.fits']['regions']['above_elevation_cut']['relative_rms_difference']
    lines=['# Gamma from maximum / standard deviation','',
        '24 October 2026. Follow-up to [calibrated dataset tests](test_calibrated_dataset.md) and [issue #1](https://github.com/peijin94/vamp/issues/1). The paper is unchanged.','',
        f"**Result:** the direct rule produces {full['ratio']['shared']['archive_bytes']/1e6:.3f} MB and {min(percentage(full['ratio']),repeat_percentage):.3f}–{max(percentage(full['ratio']),repeat_percentage):.3f}% restored-image differences in two runs. "
        f"Full gamma search gives {full['search']['shared']['archive_bytes']/1e6:.3f} MB / {percentage(full['search']):.3f}%; the inlier rule gives {full['inlier']['shared']['archive_bytes']/1e6:.3f} MB / {percentage(full['inlier']):.3f}% in its single run. "
        'The ratio rules reduce measured native search work, but show no practical whole-archive compression or full-encoding latency advantage here. Both remain optional; archive defaults are unchanged. '
        'This tests the discrete gamma set {1, 1/2, 1/4}, not every possible smoother max/std-to-gamma mapping.','',
        'The proposed direction is sensible as a tail-concentration heuristic: a larger maximum relative to the spread can justify more code levels near the center. It is not a sufficient distribution summary or an outlier detector. Outliers also increase standard deviation, and the ratio depends on the number of eligible samples and real/imaginary component structure.','',
        '## Tested rule','',
        'Subtract the separately stored float32 real/imaginary medians. From eligible float64 residual pairs r, compute population variance within each component, then pool:', '',
        r'$$\sigma^2=\frac{\operatorname{Var}(r_\mathrm{Re})+\operatorname{Var}(r_\mathrm{Im})}{2},\qquad R=\frac{\max_{j,c}|r_{j,c}|}{\sigma}.$$','',
        r'$$\gamma(R)=\begin{cases}1 & R\le R_0\\1/2 & R_0<R\le2R_0\\1/4 & R>2R_0.\end{cases}$$','',
        'Use R0=4 for the full tests; the screen also tries 3 and 6. Smaller gamma means denser levels near zero residual. Restricting the first test to the existing exact powers leaves the decoder unchanged. This is a stepwise approximation to a decreasing gamma-versus-tail-ratio rule, not a continuous-power codec.','',
        '`max-std-v1` excludes flagged, nonfinite and signed-zero-protected samples. `max-std-inlier-v1` additionally removes each range candidate’s exact outliers before computing its ratio. Thus an outlier already stored exactly need not determine the remaining samples’ gamma. Empty or zero-variance sets choose gamma 1; a restricted allowed-gamma set uses the nearest allowed log2 value.','',
        'One gamma is chosen per range candidate and shared by both components. The stored median, range, gamma and unsigned codes still reconstruct with the existing `power-midrise-v1` decoder. Selection across 4/8/12/16 bits, exceptions and exact float32 fallback still evaluates every configured error bound. Sigma here is signal spread, **not thermal noise**. Defaults remain the existing full gamma search and original uniform codebook.','',
        '## Matched-block screen','',
        'Use exactly the same 768 amplitude-ordered 128×48×1 blocks from the first, middle and last integrations for all 20 profiles. Four products are read from MS metadata. Budgets are 64 and 256 with baseline MSE ≤2B, max complex error ≤10√B and bias ≤√B. Both maximum/.995 ranges compete. The table below shows budget 256; all blocks and both budgets are retained in the JSON.','',
        '| Policy | Block payload MB | Component MSE | Screen time (s) |',
        '| --- | ---: | ---: | ---: |']
    for e in result['screen']['experiments']:
        if e['name'].endswith('b256'):
            lines.append(f"| {e['name'][:-5]} | {e['payload_bytes']/1e6:.4f} | {e['component_mse']:.3f} | {e['elapsed_seconds']:.3f} |")
    lines+=['','Payload bytes include headers and exceptions, but exclude archive metadata. Screen times use eight workers on the shared host. They are not whole-MS conversion times. Equal configured bounds do not imply equal achieved MSE.','',
        '![Ratio distribution and gamma choices](figs/gamma-ratio-selection.png)','',
        '![Matched-block payload and timing comparison](figs/gamma-ratio-screen.png)','',
        'The median observed ratio is 4.29. At R0=4, the direct rule agrees with the full search’s selected gamma in 228/588 quantized blocks; the inlier variant agrees in 380/588. The full search contains every heuristic candidate, so it cannot have a larger minimum encoded block size at identical bounds. The heuristic’s possible advantage is reduced search work, not a guaranteed better rate–distortion optimum.','',
        '## Whole-data archives and CLEAN','',
        'All three new runs use all 15 integrations, amplitude groups 128×48×1, budget 256, four-bit weights and the same calibrated-only preservation policy. Each has passed source/export verification and is compacted with the same exact HDF5 sharing. Images use the requested 2048², 4 arcmin, eight-thread CLEAN template and reference beam. All metrics and plots use the same 2,245,726 pixels above 10° at every integration.','',
        '| Policy | Shared archive MB | Serialized visibility blocks MB | Component MSE | Dirty difference % | Restored difference % | Faint-region difference % |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for e in result['experiments']:
        products=e['comparison']['products'];restored=products['image-image.fits']['regions'];dirty=products['image-dirty.fits']['regions']['above_elevation_cut']
        lines.append(f"| {e['name']} | {e['shared']['archive_bytes']/1e6:.3f} | {e['serialized_visibility_block_bytes']/1e6:.3f} | {e['component_mse']:.3f} | {100*dirty['relative_rms_difference']:.4f} | {100*restored['above_elevation_cut']['relative_rms_difference']:.4f} | {100*restored['faint_reference']['relative_rms_difference']:.4f} |")
    repeat=result['ratio_repeat']
    rr=100*repeat['versus_reference']['products']['image-image.fits']['regions']['above_elevation_cut']['relative_rms_difference']
    rf=100*repeat['versus_first']['products']['image-image.fits']['regions']['above_elevation_cut']['relative_rms_difference']
    lines+=['',f'The direct-rule repeat gives {rr:.4f}% restored-image difference against the reference and {rf:.4f}% against its first run. '
        'The previous symmetric-amplitude profile is 414.599 MB with 0.55–0.69% restored-image difference in two runs. '
        'Repeated-source/lossless controls gave 0.168–0.229%; these few repeats do not define a confidence interval. '
        'Search and inlier variants have one image run each, so small differences are not evidence of a stable scientific improvement.','',
        '![Direct max/std rule: reference, decoded and difference](figs/gamma-ratio-images.png)','',
        '![Zoomed image differences for all three gamma policies](figs/gamma-ratio-zoom.png)','',
        '![Gamma and component precision distributions](figs/gamma-ratio-distributions.png)','',
        '| Policy | Encoding only (s) | Complete verified bundle (s) | Additional compacting (s) | MS export (s) | WSClean (s) |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for e in result['experiments']:
        a=e['encoding']
        lines.append(f"| {e['name']} | {a['sum_encode_seconds_excluding_inspection']:.2f} | {a['elapsed_seconds_including_copy_inspection_verification']:.2f} | {e['shared']['elapsed_seconds']:.2f} | {e['export_verification']['export']['elapsed_seconds']:.2f} | {e['imaging']['wsclean_elapsed_seconds']:.2f} |")
    cpu={e['name']:sum(e['worker_timings_seconds'].get(k,0) for k in ('native_search_quantization_cpu','native_search_compression_cpu')) for e in result['experiments']}
    lines+=['',f"Native candidate-search and compression thread CPU summed to {cpu['search']:.1f} s for search, {cpu['ratio']:.1f} s for ratio and {cpu['inlier']:.1f} s for inlier. "
        'These are accumulated worker CPU times, not wall latency. The selector reduces this measured work, but the full encoding wall times do not demonstrate an end-to-end speedup.','',
        'Full encodes use 24 workers and overlapped other work; caches were not evicted. Screen speed does not establish end-to-end speed. Stored block streams have chunk-allocation overhead, so serialized-byte improvements need not reduce whole-file size by the same amount. The JSON includes both byte measures, actual precision/gamma distributions, error moments and worker timings.','',
        '## Reproduction and scope','',
        '```bash','export PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1',
        'python scripts/screen_gamma_ratio.py INPUT.ms NEW_SCREEN_DIR',
        'python scripts/timeslot_archive.py encode INPUT.ms NEW.bundle --config configs/calibrated-gamma-ratio.json --report encode.json',
        'python scripts/timeslot_archive.py compact NEW.bundle SHARED.h5 --report compact.json',
        'python scripts/validate_partition_export.py SHARED.h5 INPUT.ms NEW.ms --report export.json',
        'python scripts/run_wsclean.py NEW.ms --output IMAGES --data-column CORRECTED_DATA --clean-profile calibrated --threads 8 --beam-from REFERENCE/image-image.fits',
        'python scripts/compare_wsclean.py REFERENCE IMAGES --mask above10deg.fits --report comparison.json','```','',
        'Use `calibrated-gamma-inlier.json` and `calibrated-gamma-search.json` for the matched variants. New configurations use `power_gamma_policy` and `power_ratio_reference`; the decoder still reads only stored block parameters.','',
        'All 97 tests pass, including analytic ratios, threshold boundaries, protected values, range-exception exclusion, constant blocks, native/portable equivalence, independent decoding and exact fallback. All 108 original source-file checksums still match. No Q/U/V, calibration-solution, temporal or held-out-data fidelity claim is made.','',
        '[Machine-readable results](../reports/calibrated-gamma-ratio.json). Large archives, exports and FITS remain ignored under `results/calibrated-gamma-ratio/`.','']
    Path('devnotes/gamma_max_std.md').write_text('\n'.join(lines))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--screen-only',action='store_true');args=p.parse_args()
    if args.screen_only:
        screen_plots(json.loads((ROOT/'screen-final/screen.json').read_text()))
    else:
        result=collect()
        Path('reports/calibrated-gamma-ratio.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        screen_plots(result['screen']);image_plots(result);write_note(result)
