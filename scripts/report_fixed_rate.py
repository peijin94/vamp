"""Report actual <=25%-payload finalists and repeated calibrated CLEAN images."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path

import numpy as np
import h5py

from varquanms.imaging import compare_run, file_sha256, read_fits
from varquanms.codec import unpack_header, decode_block

ROOT=Path('results/calibrated-fixed25')
REF=Path('results/calibrated-template/reference')
MASK=Path('results/calibrated-template/above10deg.fits')
FIGS=Path('devnotes/figs')


def payload_distribution(path, *, by_class=False):
    sizes=Counter();samples=Counter();zero_samples=zero_bytes=0
    class_samples={'auto':Counter(),'cross':Counter()}
    with h5py.File(path,'r') as h5:
        manifest=json.loads(h5['manifest'][...].tobytes())
        for part in manifest['partitions']:
            group=h5[part['group']]
            spec=json.loads(group['manifest'][...].tobytes())
            classes={}
            if by_class:
                from types import SimpleNamespace
                from varquanms.archive import read_column
                reader=SimpleNamespace(h5=group)
                for table in spec['tables']:
                    columns=table['columns']
                    for visibility in columns.values():
                        if visibility['encoding']!='quantized':continue
                        auto=(read_column(reader,columns['ANTENNA1'],0,table['row_count'])==
                              read_column(reader,columns['ANTENNA2'],0,table['row_count']))
                        segments=(visibility['cell_shape'][0]+visibility['channel_segment']-1)//visibility['channel_segment']
                        for product in visibility['products']:
                            perm=(np.arange(table['row_count']) if product['permutation'] is None
                                  else group[product['permutation']][...])
                            boundaries=group[product['group_boundaries']][...]
                            ordered=auto[perm]
                            sums=np.r_[np.int64(0),np.cumsum(ordered,dtype='i8')]
                            counts=sums[boundaries[1:]]-sums[boundaries[:-1]]
                            sizes_by_group=np.diff(boundaries)
                            assert np.all(sizes_by_group>0) and np.all((counts==0)|(counts==sizes_by_group))
                            classes[product['payload']]=np.repeat(counts==sizes_by_group,segments)
            for entry in spec['inventory']:
                if entry['role']!='visibility_block_payload':continue
                ds=group[entry['path']]
                parent=ds.parent;offsets=parent['offsets'][...]
                # Millions of tiny per-baseline blocks must not become millions
                # of HDF5 reads. Keep the report's cache bounded to 64 MiB.
                payload=memoryview(ds[...]) if ds.size <= 64*1024**2 else None
                for index,(lo,hi) in enumerate(zip(offsets[:-1],offsets[1:])):
                    blob=payload[int(lo):int(hi)] if payload is not None else ds[int(lo):int(hi)].tobytes()
                    bits,_,_,count,_,_=unpack_header(blob)
                    label=str(bits) if bits else 'float32'
                    sizes[label]+=len(blob);samples[label]+=count
                    if by_class:class_samples['auto' if classes[entry['path']][index] else 'cross'][label]+=count
                    if bits==0 and np.all(decode_block(blob,(count,)).view('<u4')==0):
                        zero_samples+=count;zero_bytes+=len(blob)
    result=dict(serialized_bytes_by_bits=dict(sizes),samples_by_bits=dict(samples),
                exact_positive_zero_samples=zero_samples,exact_positive_zero_bytes=zero_bytes)
    if by_class:
        assert class_samples['auto']+class_samples['cross']==samples
        result['samples_by_baseline_class_and_bits']={key:dict(value) for key,value in class_samples.items()}
    return result


def fraction(comparison,product='image',region='above_elevation_cut'):
    return comparison['products'][f'image-{product}.fits']['regions'][region]['relative_rms_difference']


def symmetric_ticks(scale):
    positive=10.**np.arange(np.ceil(np.log10(.01*scale)),np.floor(np.log10(scale))+1)
    return np.r_[-positive[::-1],0,positive]


def collect(partial=False,reuse_preview=False):
    screen=json.loads((ROOT/'screen/screen.json').read_text())
    plan=json.loads((ROOT/'screen/plan.json').read_text())
    cached={e['name']:e for e in json.loads((ROOT/'preview.json').read_text())['experiments']} if reuse_preview else {}
    cases=[]
    for selected in plan:
        name=selected['name'];read=lambda suffix:json.loads((ROOT/(name+'-'+suffix+'.json')).read_text())
        if partial and not (ROOT/(name+'-comparison.json')).exists():continue
        rate=read('rate');assert rate['within_target'] and rate['fraction']<=.25
        manifest=json.loads((ROOT/(name+'.bundle')/'manifest.json').read_text())
        assert manifest['config']==selected['config']
        shared=read('tight-compact') if (ROOT/(name+'-tight-compact.json')).exists() else read('compact')
        export=read('export');comparison=read('comparison')
        assert shared['all_partition_inventories_verified']
        assert shared['archive_bytes']==sum(shared['stored_bytes_by_category'].values())
        assert export['status']=='passed' and export['source_files_unchanged']
        assert comparison['comparison_mask']['sha256']==file_sha256(MASK)
        categories=shared['stored_bytes_by_category']
        assert categories['visibility_block_payload'] <= rate['payload_bytes']
        rate['physical_payload_bytes']=categories['visibility_block_payload']
        rate['physical_payload_fraction']=rate['physical_payload_bytes']/rate['original_complex64_bytes']
        rate['payload_plus_permutations_and_indexes_bytes']=sum(categories.get(k,0) for k in
            ('visibility_block_payload','permutations','visibility_indexes'))
        if name in cached:
            e=cached[name]
            assert all(e[k]==v for k,v in selected.items())
            assert e['shared']==shared and e['rate']==rate and e['export_verification']==export
            assert e['comparison']==comparison and e['encoding']==read('encode')
            assert e['visibility']==read('visibility')
            assert e['archive_sha256']==file_sha256(Path(shared['archive']))
            assert e['imaging']==json.loads((ROOT/(name+'-images')/'complete.json').read_text())
            assert e['repeat']['imaging']==json.loads((ROOT/(name+'-repeat-images')/'complete.json').read_text())
            for c in (e['comparison'],e['repeat']['versus_reference'],e['repeat']['versus_first']):
                assert c['comparison_mask']['sha256']==file_sha256(MASK)
                for product in c['products'].values():
                    assert file_sha256(product['reference'])==product['reference_file_sha256']
                    assert file_sha256(product['candidate'])==product['candidate_file_sha256']
            cases.append(e)
            continue
        samples=Counter();bits=Counter();gammas=Counter()
        for part in manifest['partitions']:
            assert all(q['original_bits_exact'] for q in part['summary']['weights_quality'].values())
            for product in part['summary']['visibility_block_statistics']['CORRECTED_DATA']['products']:
                samples.update(product['complex_sample_bit_width_histogram']);bits.update(product['component_bit_width_histogram'])
                gammas.update(product['power_gamma_block_histogram'])
        assert sum(samples.values())==931920*48*4
        e=dict(**selected,rate=rate,shared=shared,encoding=read('encode'),export_verification=export,
            comparison=comparison,visibility=read('visibility'),sample_bits=dict(samples),block_bits=dict(bits),gamma_blocks=dict(gammas),
            imaging=json.loads((ROOT/(name+'-images')/'complete.json').read_text()),
            archive_sha256=file_sha256(Path(shared['archive'])))
        e['payload_distribution']=payload_distribution(shared['archive'])
        assert e['payload_distribution']['samples_by_bits']==dict(samples)
        assert sum(e['payload_distribution']['serialized_bytes_by_bits'].values())==rate['payload_bytes']
        repeat=ROOT/(name+'-repeat-images')
        if (repeat/'complete.json').exists():
            e['repeat']=dict(versus_reference=compare_run(REF,repeat,MASK),versus_first=compare_run(ROOT/(name+'-images'),repeat,MASK),
                             imaging=json.loads((repeat/'complete.json').read_text()))
            e['worst_restored_fraction_two_runs']=max(fraction(comparison),fraction(e['repeat']['versus_reference']))
        elif not partial:raise ValueError('Repeat required: '+name)
        cases.append(e)
    result=dict(schema='varquanms.fixed25-study/1',status='partial' if partial else 'complete',
        original_complex64_bytes=1431429120,payload_budget_bytes=357857280,screen=screen,experiments=cases,
        selection_objective='Primary empirical score: larger restored-image RMS difference of two CLEAN runs, all above 10 degrees. Also report dirty, faint, bright, peak and aperture metrics; two runs do not establish a confidence interval.',
        scope='Payload fraction includes serialized block factors/headers/exceptions and lossless fallback. Full archive bytes are separate. All products are archived; images are joint Stokes I only.',
        limitations=['Best means best among the fully imaged finalists on this observation, not a global parameter optimum or a held-out guarantee.',
                     'Pilot uses grouping-specific rank samples; full archives validate actual byte counts and scientific metrics.',
                     'Independent eight-thread CLEAN may amplify tiny input or numerical changes; repeated source controls are reported in the companion error-limit study.',
                     'Relative means block-normalized RMS, not a samplewise percentage bound. No inferred thermal-noise guarantee.'])
    if not partial:
        result['empirical_best_worst_of_two']=min(cases,key=lambda e:e['worst_restored_fraction_two_runs'])['name']
        result['best_first_dirty']=min(cases,key=lambda e:fraction(e['comparison'],'dirty'))['name']
        broad=json.loads(Path('reports/calibrated-error-limits.json').read_text())
        result['eligible_lower_rate_controls']=[dict(name=e['name'],config=e['config'],
            payload_bytes=e['serialized_visibility_bytes'],
            physical_payload_bytes=e['shared']['stored_bytes_by_category']['visibility_block_payload'],
            payload_fraction=e['serialized_visibility_bytes']/1431429120,
            first_restored_fraction=fraction(e['comparison']),
            repeat_restored_fraction=fraction(e['repeat']['versus_reference']),
            worst_restored_fraction_two_runs=max(fraction(e['comparison']),fraction(e['repeat']['versus_reference'])),
            dirty_fraction=fraction(e['comparison'],'dirty'))
            for e in broad['experiments'] if e['serialized_visibility_bytes']<=357857280]
        all_repeated=cases+result['eligible_lower_rate_controls']
        result['best_across_repeated_feasible_tests']=min(all_repeated,key=lambda e:e['worst_restored_fraction_two_runs'])['name']
        historical=[]
        for e in json.loads(Path('reports/calibrated-template-tests.json').read_text())['experiments']:
            path=Path('results/calibrated-template')/(e['name']+'.bundle')/'manifest.json'
            if not path.exists():continue
            manifest=json.loads(path.read_text())
            payload=sum(p['serialized_block_bytes'] for part in manifest['partitions']
                for p in part['summary']['visibility_block_statistics'].get('CORRECTED_DATA',{}).get('products',[]))
            if 0<payload<=357857280:
                historical.append(dict(name=e['name'],payload_bytes=payload,
                    payload_fraction=payload/1431429120,first_restored_fraction=fraction(e['comparison']),
                    cannot_beat_selected_worst_of_two=fraction(e['comparison'])>=min(c['worst_restored_fraction_two_runs'] for c in all_repeated)))
        result['historical_lower_rate_controls']=historical
        assert all(e['physical_payload_bytes']<=357857280 for e in result['eligible_lower_rate_controls'])
        assert all(e['cannot_beat_selected_worst_of_two'] for e in historical)
        result['control_rate_definition']='Serialized codec bytes; older files may retain streaming padding. Repacking preserves codec bytes. All five new finalists also pass the physical-payload limit.'
        imported=json.loads(Path('reports/calibrated-import.json').read_text())
        assert all(file_sha256(Path(imported['extracted_ms'])/n)==v['sha256'] for n,v in imported['files'].items())
        result['source_integrity']=dict(original_import_checksums_match=True,files_checked=len(imported['files']))
    return result


def plots(result, *, root=ROOT, prefix='fixed25', output=FIGS, include_best=True):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/varquanms-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from matplotlib.ticker import StrMethodFormatter
    from astropy.io import fits
    cases=result['experiments']
    ceiling=(100*result['payload_budget_bytes']/result['original_complex64_bytes']
             if result.get('payload_budget_bytes') is not None else None)
    residual_rms=result.get('reference_residual_statistics',{}).get('rms')
    image_fraction=lambda c,product='image': (c['products'][f'image-{product}.fits']['regions']['above_elevation_cut']['difference_rms']/residual_rms
        if residual_rms else fraction(c,product))
    labels=[e.get('plot_label',f"{e['config']['baseline_group']}×{e['config']['channel_segment']}\n"+
            ('search' if e['config']['power_gamma_policy']=='search' else 'max/std')+'\n'+e['family']) for e in cases]
    fig,axes=plt.subplots(1,2,figsize=(max(11,1.65*len(cases)),5),constrained_layout=True)
    for x,e in enumerate(cases):
        axes[0].bar(x,100*e['rate']['fraction'],color='#4393c3')
        a=100*image_fraction(e['comparison']);b=100*image_fraction(e['repeat']['versus_reference'])
        axes[1].plot([x,x],[a,b],color='#b2182b',linewidth=2)
        axes[1].scatter([x,x],[a,b],color='#b2182b',label='Restored: two runs' if x==0 else None)
        axes[1].scatter(x,100*image_fraction(e['comparison'],'dirty'),marker='s',color='#2166ac',label='Dirty: first run' if x==0 else None)
    if ceiling is not None:axes[0].axhline(ceiling,color='black',linestyle='--')
    upper=ceiling if ceiling is not None else max(100*e['rate']['fraction'] for e in cases)
    axes[0].set(ylabel='Serialized visibility bytes / complex64 (%)',ylim=(0,1.12*upper))
    axes[1].set(ylabel='Difference RMS / reference residual RMS (%)' if residual_rms else 'RMS difference / reference RMS (%)',yscale='log');axes[1].legend(loc='upper left')
    axes[1].yaxis.set_major_formatter(StrMethodFormatter('{x:g}'))
    axes[1].yaxis.set_minor_formatter(StrMethodFormatter('{x:g}'))
    for ax in axes:ax.set(xticks=range(len(cases)),xticklabels=labels);ax.tick_params(axis='x',labelsize=8);ax.grid(axis='y',alpha=.2)
    fig.suptitle(result.get('plot_title','Fixed visibility-payload budget; same 15-integration MS and above-10° image mask'))
    fig.savefig(output/(prefix+'-rate.png'),dpi=170);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(13,max(5,.8*len(cases))),constrained_layout=True);y=np.arange(len(cases))
    for ax,key,title in zip(axes,('samples_by_bits','serialized_bytes_by_bits'),('All complex samples (%)','Serialized payload bytes (%)')):
        left=np.zeros(len(cases))
        precision=result.get('precision_labels',('3','4','6','8','12','16','float32'))
        colors=(list(plt.get_cmap('viridis')(np.linspace(.1,.9,len(precision)-1)))+['#aaa']
                if 'precision_labels' in result else ['#fddbc7','#ef8a62','#92c5de','#4393c3','#2166ac','#542788','#aaa'])
        for bits,color in zip(precision,colors):
            values=np.array([e['payload_distribution'][key].get(bits,0)/sum(e['payload_distribution'][key].values())*100 for e in cases])
            ax.barh(y,values,left=left,label=bits,color=color);left+=values
        ax.set(yticks=y,yticklabels=labels,xlabel=title,xlim=(0,100));ax.invert_yaxis()
    axes[1].legend(title='Bits/component; float32 = exact',ncol=min(5,len(precision)) if len(precision)>7 else 7,loc='upper center',bbox_to_anchor=(.5,-.15),fontsize=8)
    fig.suptitle('Bit-size distribution: exact zero blocks can contain many samples but few bytes')
    fig.savefig(output/(prefix+'-bits.png'),dpi=170);plt.close(fig)
    reference,_=read_fits(REF/'image-image.fits');mask=fits.getdata(MASK).astype(bool)
    deltas=[read_fits(root/(e['name']+'-images')/'image-image.fits')[0]-reference for e in cases]
    y,x=cases[0]['comparison']['products']['image-image.fits']['reference_peak_yx']
    for zoom in (False,True):
        crop=np.s_[y-64:y+64,x-64:x+64] if zoom else np.s_[:,:]
        domain=mask[crop];lim=max(float(np.max(np.abs(d[crop][domain]))) for d in deltas)
        if len(cases)>=4:
            rows=(len(cases)+3)//3
            fig,axes=plt.subplots(rows,3,figsize=(12,4*rows),constrained_layout=True)
            flat=axes.ravel()
            for unused in flat[len(cases)+1:]:fig.delaxes(unused)
            axes=flat[:len(cases)+1]
        else:
            fig,axes=plt.subplots(1,len(cases)+1,figsize=(3.2*(len(cases)+1),4.3),constrained_layout=True)
        scale=float(np.max(np.abs(reference[crop][domain])))
        im=axes[0].imshow(np.where(domain,reference[crop],np.nan),origin='lower',cmap='viridis',norm=SymLogNorm(.01*scale,vmin=-scale,vmax=scale))
        axes[0].set_title('Reference');fig.colorbar(im,ax=axes[0],shrink=.7,label='Jy/beam',ticks=symmetric_ticks(scale))
        for ax,label,delta in zip(axes[1:],labels,deltas):
            im=ax.imshow(np.where(domain,delta[crop],np.nan),origin='lower',cmap='bwr',vmin=-lim,vmax=lim);ax.set_title(label+' − reference',fontsize=9)
        fig.colorbar(im,ax=list(axes[1:]),shrink=.7,label='Difference Jy/beam; common linear scale')
        for ax in axes:ax.set(xlabel='Crop pixel x' if zoom else 'Pixel x',ylabel='Crop pixel y' if zoom else 'Pixel y')
        fig.suptitle('First CLEAN run for every finalist; above 10° only'+(' • reference-peak zoom' if zoom else ''))
        fig.savefig(output/(prefix+('-zoom.png' if zoom else '-images.png')),dpi=170);plt.close(fig)
    if not include_best:
        return
    best=next(e for e in cases if e['name']==result['empirical_best_worst_of_two'])
    restored,_=read_fits(root/(best['name']+'-images')/'image-image.fits')
    for zoom in (False,True):
        crop=np.s_[y-64:y+64,x-64:x+64] if zoom else np.s_[:,:]
        domain=mask[crop];a=reference[crop];b=restored[crop];delta=b-a
        scale=max(float(np.max(np.abs(a[domain]))),float(np.max(np.abs(b[domain]))))
        limit=float(np.max(np.abs(delta[domain])))
        fig,axes=plt.subplots(1,3,figsize=(12,4.5),constrained_layout=True)
        for ax,values,title in zip(axes[:2],(a,b),('Uncompressed reference','Selected quantized archive')):
            im=ax.imshow(np.where(domain,values,np.nan),origin='lower',cmap='viridis',
                         norm=SymLogNorm(.01*scale,vmin=-scale,vmax=scale))
            ax.set_title(title)
        fig.colorbar(im,ax=list(axes[:2]),shrink=.7,label='Jy/beam; shared scale',ticks=symmetric_ticks(scale))
        im=axes[2].imshow(np.where(domain,delta,np.nan),origin='lower',cmap='bwr',vmin=-limit,vmax=limit)
        axes[2].set_title('Quantized − reference')
        fig.colorbar(im,ax=axes[2],shrink=.7,label='Difference Jy/beam')
        for ax in axes:ax.set(xlabel='Crop pixel x' if zoom else 'Pixel x',ylabel='Crop pixel y' if zoom else 'Pixel y')
        fig.suptitle(best['layout']+' / '+best['family']+'; first CLEAN run; above 10°'+(' • peak zoom' if zoom else ''))
        fig.savefig(output/(prefix+('-best-zoom.png' if zoom else '-best-images.png')),dpi=170);plt.close(fig)


def note(result):
    best=next(e for e in result['experiments'] if e['name']==result['empirical_best_worst_of_two'])
    config=best['config']
    Path('configs/calibrated-fixed25-best.json').write_text(json.dumps(config,indent=2,allow_nan=False)+'\n')
    lines=['# Fixed 25% visibility-payload budget','',
        f"Among the fully imaged finalists, **{best['layout']} with {best['family']} limits** gives the smallest worse-of-two restored-image difference. Its physical visibility payload is **{best['rate']['physical_payload_bytes']/1e6:.3f} MB ({100*best['rate']['physical_payload_fraction']:.3f}%)**; the complete file is **{best['shared']['archive_bytes']/1e6:.3f} MB**. The two restored-image differences are **{100*fraction(best['comparison']):.4f}% and {100*fraction(best['repeat']['versus_reference']):.4f}%**. This measures agreement with the uncompressed Stokes I reference, not accuracy against sky truth.",'',
        '[Exact opt-in configuration](../configs/calibrated-fixed25-best.json). Defaults remain unchanged. This calibrated-only profile omits original DATA/MODEL_DATA and zeros flagged visibilities while preserving their flags; original CORRECTED_DATA values are used independently for every experiment. Four-bit weights are exact for this dataset.','',
        'Follow-up to [absolute versus relative error limits](absolute_relative_error.md). The target is **25% of original complex64 visibility bytes**, equivalent to 25% of the two original float32 component streams. It means 4× smaller, not a 25% saving.','',
        'This MS contains 178,928,640 complex samples: 1,431,429,120 original bytes. The payload ceiling is **357,857,280 bytes (357.857 MB)**. Payload includes codes, medians, scales, gamma, block headers, exact exceptions and float32 fallback after conventional lossless compression. Whole-archive coordinates, flags, weights, maps, metadata and quality arrays are counted separately.','',
        '## Tested configurations','',
        'The pilot compares five layout/codec choices and three error rules, tuning thresholds with 14 logarithmic byte-count refinements. It targets 24.8% to allow for sampling differences. Three 128×48 candidates test absolute, relative and hybrid rules directly; the lowest-MSE alternative and a smaller-group control are also imaged. Every full finalist’s serialized block stream had to fit 25% before imaging. The final compacted physical streams are independently checked against the same ceiling. This is a finite experimental comparison, not a proof of the global optimum.','',
        '| Finalist | Baselines × channels | Ordering | Gamma rule | Absolute MSE cap | Relative component RMS fraction |',
        '| --- | --- | --- | --- | ---: | ---: |']
    for e in result['experiments']:
        c=e['config'];relative=(c['relative_bounds'] or {}).get('max_component_rmse')
        lines.append(f"| {e['layout']} / {e['family']} | {c['baseline_group']} × {c['channel_segment']} | {c['ordering']} | {c['power_gamma_policy']} | {c['max_component_mse'] if c['max_component_mse'] is not None else '—'} | {relative if relative is not None else '—'} |")
    lines+=['','All groups contain one correlation product; all spectra retain their original channels. Power levels use stored real/imaginary medians, float32 range, gamma in {1, 1/2, 1/4}, exact exceptions and 4/8/12/16-bit choices. Ratio selection uses R0=4. The error-rule definitions and maximum/bias/baseline guards are specified in the companion note; exact configurations are retained in the JSON.','',
        '## Full-data images and bytes','',
        '| Finalist | Payload MB (%) | Whole archive MB | Dirty difference % | Restored run 1 % | Restored run 2 % | Worst faint-region difference % |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for e in result['experiments']:
        comparison=e['comparison'];repeat=e['repeat']['versus_reference']
        lines.append(f"| {e['layout']} / {e['family']} | {e['rate']['payload_bytes']/1e6:.3f} ({100*e['rate']['fraction']:.3f}%) | {e['shared']['archive_bytes']/1e6:.3f} | {100*fraction(comparison,'dirty'):.4f} | {100*fraction(comparison):.4f} | {100*fraction(repeat):.4f} | {100*max(fraction(comparison,'image','faint_reference'),fraction(repeat,'image','faint_reference')):.4f} |")
    lines+=['',f"Lowest worst-of-two restored-image difference: **{result['empirical_best_worst_of_two']}**. Lowest first-run dirty-image difference: **{result['best_first_dirty']}**. These are empirical rankings on this observation; two repeats do not establish a confidence interval or a stable universal winner.",'',
        'For this budget, use the absolute-limit configuration as a starting point: the tested relative rule has a larger dirty-image difference and a larger worse-of-two restored difference. Hybrid limits do not improve this image comparison. Full gamma search reduces dirty-image error but gives a much less consistent CLEAN result. The 32×24 grouping control is also slightly worse than 128×48 in both restored runs. These results do not establish that absolute limits always win, or that larger groups always win.','',
        'The source-only CLEAN repeat differences are about 0.229% and 0.236%; the same-policy lossless control gives 0.168%. They show the limited reproducibility of this eight-thread, independently masked workflow. They are not a confidence interval or a universal noise floor, and are not subtracted from quantized-image differences.','',
        'All images use the requested 2048², 4 arcmin, eight-thread multiscale Stokes I template, identical reference beam and fixed above-10° mask. Each archive is independently decoded and verified. Weights are exact on this dataset, so visibility-error allocation is the changing factor. No approval of a scientific error tolerance is implied by a low image-difference score.','',
        '![Fixed-budget rate and image differences](figs/fixed25-rate.png)','',
        '![Selected archive: reference, decoded and difference images](figs/fixed25-best-images.png)','',
        '![Selected archive: reference-peak zoom](figs/fixed25-best-zoom.png)','',
        '![Full restored-image differences](figs/fixed25-images.png)','',
        '![Reference-peak zoom](figs/fixed25-zoom.png)','',
        '![Bit-size distribution](figs/fixed25-bits.png)','',
        '### Smaller archives and earlier controls','',
        'The target is a ceiling, so the six lower-rate profiles from the companion study also qualify. They were each imaged twice and repacked with the same lossless storage policy; their physical payloads also fit 25%.','',
        '| Lower-rate rule | Payload % | First restored % | Repeat restored % |',
        '| --- | ---: | ---: | ---: |']
    for e in result['eligible_lower_rate_controls']:
        lines.append(f"| {e['name']} | {100*e['payload_fraction']:.3f} | {100*e['first_restored_fraction']:.4f} | {100*e['repeat_restored_fraction']:.4f} |")
    lines+=['',f"Best worst-of-two result across these controls and the five near-budget finalists: **{result['best_across_repeated_feasible_tests']}**.",'',
        'Earlier serialized-payload controls are also retained in the JSON, including the parallel-hand priority and six-bit experiments. Each has a first-run restored difference greater than the selected candidate’s worse repeat, so adding a second run cannot improve its worst-of-two score below that candidate. Their original physical container sizes are not substituted for codec bytes.','',
        '### Physical accounting and timings','',
        'The original extensible HDF5 streams allocated unused tail padding. Lossless compaction writes their known lengths without that padding. All scientific manifests and array checksums match the already imaged archives; the original files remain available. The physical payload meets 25%, including block headers, factors and exceptions. Permutations, indexes and other categories remain explicit rather than disappearing into a claimed payload ratio.','',
        '| Finalist | Payload + permutation/index MB | Encoding only s | Verified bundle s | Compaction s | Export s | First CLEAN s |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for e in result['experiments']:
        lines.append(f"| {e['layout']} / {e['family']} | {e['rate']['payload_plus_permutations_and_indexes_bytes']/1e6:.3f} | {e['encoding']['sum_encode_seconds_excluding_inspection']:.2f} | {e['encoding']['elapsed_seconds_including_copy_inspection_verification']:.2f} | {e['shared']['elapsed_seconds']:.2f} | {e['export_verification']['export']['elapsed_seconds']:.2f} | {e['imaging']['wsclean_elapsed_seconds']:.2f} |")
    lines+=['','Times cover all 15 integrations. Runs overlap on a shared host with no cache eviction; these are observed durations, not isolated speed rankings or a 15-second full-conversion claim. The JSON retains hardware/software and commands through the companion report and per-stage records.','',
        'For the selected configuration, the measured byte breakdown is:','',
        '| Category | Decimal MB |','| --- | ---: |']
    for key,value in best['shared']['stored_bytes_by_category'].items():
        lines.append(f'| {key} | {value/1e6:.4f} |')
    lines+=['','The sample-weighted precision distribution is different from the byte distribution because exact zero blocks compress well:','',
        '| Bits per component | Sample fraction % | Payload byte fraction % |','| --- | ---: | ---: |']
    dist=best['payload_distribution']
    for bits in ('4','8','12','16','float32'):
        lines.append(f"| {bits} | {100*dist['samples_by_bits'].get(bits,0)/178928640:.4f} | {100*dist['serialized_bytes_by_bits'].get(bits,0)/best['rate']['payload_bytes']:.4f} |")
    lines+=['',f"Exact all-positive-zero blocks contain {dist['exact_positive_zero_samples']:,} samples and occupy only {dist['exact_positive_zero_bytes']/1e6:.4f} MB. The calibrated-only flagged-value policy and ordinary lossless compression contribute to the measured rate; the 4× payload result is not a claim that every retained value uses a fixed 8-bit component code.",'',
        '## Reproduction','',
        '```bash','export PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1',
        'python scripts/screen_fixed_rate.py INPUT.ms NEW_SCREEN_DIR',
        'python scripts/run_error_limits.py --source INPUT.ms --plan NEW_SCREEN_DIR/plan.json --root NEW_RESULT_DIR',
        'python scripts/run_error_limits.py --source INPUT.ms --plan NEW_SCREEN_DIR/plan.json --root NEW_RESULT_DIR --repeat-only',
        '```','',
        '[Machine-readable results](../reports/calibrated-fixed25.json) include thresholds, all pilot trials, measured byte categories, encode/export/CLEAN timings, commands, image metrics and checksums. Large files remain ignored under `results/calibrated-fixed25/`. Source files retain the original 108 import checksums. The paper and root README are unchanged.','']
    Path('devnotes/fixed25.md').write_text('\n'.join(lines))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--partial',action='store_true')
    p.add_argument('--reuse-preview',action='store_true',help='Reuse previously measured cases only after rechecking configurations, reports and every archive/FITS hash')
    args=p.parse_args()
    result=collect(args.partial,args.reuse_preview)
    path=ROOT/'preview.json' if args.partial else Path('reports/calibrated-fixed25.json')
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for e in result['experiments']:print(e['name'],e['rate']['fraction'],100*fraction(e['comparison'],'dirty'),100*fraction(e['comparison']))
    if not args.partial:plots(result);note(result)


if __name__=='__main__':main()
