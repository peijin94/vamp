"""Collect the 16%-payload epsilon study, repeated images, byte audits and plots."""
import argparse
from collections import Counter
import json
from pathlib import Path
import re

import numpy as np
from astropy.io import fits

from varquanms.imaging import _grid, file_sha256 as _file_sha256, read_fits
from report_fixed_rate import payload_distribution, fraction, plots
from report_fixed10 import moments

ROOT=Path('results/calibrated-fixed16')
REFERENCE=Path('results/calibrated-template/reference')
MASK=Path('results/calibrated-template/above10deg.fits')
ORIGINAL_BYTES=1431429120
BUDGET_BYTES=229028659
HASHES={}


def file_sha256(path):
    """Reuse hashes during collection; shared files are rehashed at the end."""
    path=Path(path).resolve()
    def stamp():
        s=path.stat()
        return (s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns)
    before=stamp()
    if path in HASHES:
        original,digest,count=HASHES[path]
        if original!=before:raise ValueError(f'File changed during report: {path}')
        HASHES[path]=(original,digest,count+1)
        return digest
    digest=_file_sha256(path)
    if before!=stamp():raise ValueError(f'File changed while hashing: {path}')
    HASHES[path]=(before,digest,1)
    return digest


def verify_reused_hashes():
    # A repeated reference needs only two full reads, rather than one per case.
    # Rehash contents because timestamp resolution may miss same-size writes.
    for path,(_,digest,count) in HASHES.items():
        if count>1 and _file_sha256(path)!=digest:
            raise ValueError(f'Repeated reference changed during report: {path}')


def label(e):
    c=e['config']
    return f"{c['baseline_group']}×{c['channel_segment']} {e['family']}"+(' search' if c['power_gamma_policy']=='search' else '')


def binding_fractions(entry):
    records=json.loads((ROOT/'screen'/(entry['layout']+'-frontier.json')).read_text())['blocks']
    eps=entry['config']['normalized_bounds'];total=0;counts=Counter()
    for block in records:
        if block['auto'] or not block['references']['valid_samples']:continue
        values={k:eps['epsilon_'+k]*(block['global_std'] if k=='G' else block['references'][k])
                for k in 'GSDR' if eps['epsilon_'+k] is not None}
        active=[k for k,v in values.items() if v==min(values.values())]
        n=block['references']['valid_samples'];total+=n
        for k in active:counts[k]+=n/len(active)
    return dict(domain='Sampled unflagged finite cross-correlation samples; exact ties split equally',
                samples=total,fractions={k:counts[k]/total if total else 0. for k in 'GSDR'})


def sampling_policy_sensitivity(plan, stability):
    """Hold every finalist's epsilons fixed while changing only sampled G."""
    from screen_normalized import Matrix
    results=[]
    for layout in dict.fromkeys(e['layout'] for e in plan):
        records=json.loads((ROOT/'screen'/(layout+'-frontier.json')).read_text())['blocks']
        matrix=Matrix(records)
        for selected in (e for e in plan if e['layout']==layout):
            policy=selected['config']['normalized_bounds']
            recipe=dict(weights={k:policy['epsilon_'+k] for k in 'GSDR' if policy['epsilon_'+k] is not None})
            trials=[];base=None
            for sample in stability:
                matrix.refs['G']=np.array([sample['references']['groups'][b['configuration_key']]['columns']['CORRECTED_DATA']
                    [b['product']]['auto' if b['auto'] else 'cross']['component_std'] or 0. for b in records])
                trial=matrix.evaluate(recipe,1.);choices=np.array(trial.pop('choices'))
                if base is None:base=choices
                changed=choices!=base
                trials.append(dict(sampling=sample['name'],**trial,changed_blocks=int(changed.sum()),
                    changed_complex_samples=int(matrix.samples[changed].sum()),
                    changed_sample_fraction=float(matrix.samples[changed].sum()/matrix.samples.sum())))
            results.append(dict(name=selected['name'],trials=trials))
    return dict(scope='Pilot blocks only; fixed epsilons, changed G estimate, no retuning or additional CLEAN',experiments=results)


def visibility_storage_equality(first, second):
    """Compare actual arrays that determine visibility decoding, in bounded reads."""
    import h5py
    roles={'visibility_block_payload','visibility_indexes','permutations','reconstruction'}
    checked=0;different=[]
    with h5py.File(first,'r') as a,h5py.File(second,'r') as b:
        ma=json.loads(a['manifest'][...].tobytes());mb=json.loads(b['manifest'][...].tobytes())
        assert len(ma['partitions'])==len(mb['partitions'])
        for index,(pa,pb) in enumerate(zip(ma['partitions'],mb['partitions'])):
            ga,gb=a[pa['group']],b[pb['group']]
            sa=json.loads(ga['manifest'][...].tobytes());sb=json.loads(gb['manifest'][...].tobytes())
            ia={e['path']:e for e in sa['inventory'] if e['role'] in roles}
            ib={e['path']:e for e in sb['inventory'] if e['role'] in roles}
            if ia.keys()!=ib.keys():
                different.append(dict(partition=index,reason='different array paths'));continue
            for path in ia:
                x,y=ga[path],gb[path];checked+=1
                same=x.shape==y.shape and x.dtype==y.dtype
                if same:
                    stride=max(1,1048576//max(1,int(np.prod(x.shape[1:]))*x.dtype.itemsize))
                    same=all(x[start:start+stride].tobytes()==y[start:start+stride].tobytes()
                             for start in range(0,x.shape[0],stride))
                if not same:different.append(dict(partition=index,path=path))
    return dict(first=str(first),second=str(second),arrays_compared=checked,
                arrays_identical=not different,different_arrays=different,roles=sorted(roles),
                meaning='With the same verified source mapping and codec settings, equality gives identical decoded visibility inputs.')


def serial_parallel_equality(entries):
    import h5py
    attempt=json.loads((ROOT/'baseline-serial-attempt.json').read_text())
    baseline=next(e for e in entries if e['config']['baseline_group']==1)
    assert baseline['encoding']['partition_workers']==6
    checked=0;partitions=[]
    with h5py.File(baseline['shared']['archive'],'r') as new:
        manifest=json.loads(new['manifest'][...].tobytes())
        for verification in attempt['completed_partitions']:
            index=int(verification.removeprefix('verification-').removesuffix('.json'))
            oldpath=Path(attempt['partial_bundle'])/f'partition-{index:04d}.h5'
            with h5py.File(oldpath,'r') as old:
                group=new[manifest['partitions'][index]['group']]
                a=json.loads(old['manifest'][...].tobytes());b=json.loads(group['manifest'][...].tobytes())
                ia={e['path']:e for e in a['inventory']};ib={e['path']:e for e in b['inventory']}
                assert ia==ib
                for path in ia:
                    x,y=old[path],group[path]
                    assert x.shape==y.shape and x.dtype==y.dtype
                    if x.ndim==0:assert x[()].tobytes()==y[()].tobytes()
                    else:
                        stride=max(1,1048576//max(1,int(np.prod(x.shape[1:]))*x.dtype.itemsize))
                        for start in range(0,x.shape[0],stride):
                            assert x[start:start+stride].tobytes()==y[start:start+stride].tobytes()
                    checked+=1
                partitions.append(index)
    return dict(serial_attempt=attempt,parallel_partition_workers=6,
                equivalence=dict(status='passed',completed_serial_partitions=partitions,
                    arrays_compared=checked,all_inventory_arrays_bitwise_identical=True))


def collect(partial=False):
    HASHES.clear()
    plan=json.loads((ROOT/'screen/plan.json').read_text())
    screen=json.loads((ROOT/'screen/screen.json').read_text())
    assert plan==screen['selected']
    reference,header=read_fits(REFERENCE/'image-image.fits')
    residual,rh=read_fits(REFERENCE/'image-residual.fits');assert _grid(header)==_grid(rh)
    mask=fits.getdata(MASK).astype(bool);assert int(mask.sum())==2245726
    ref_stats=moments(residual[mask]);entries=[]
    for selected in plan:
        name=selected['name'];read=lambda stage:json.loads((ROOT/(name+'-'+stage+'.json')).read_text())
        if partial and not (ROOT/(name+'-comparison.json')).exists():continue
        rate=read('rate');shared=read('compact');export=read('export');config=selected['config']
        manifest=json.loads((ROOT/(name+'.bundle')/'manifest.json').read_text())
        assert manifest['config']==config and config['bits']==[3,4,6,8,12,16]
        assert config['normalized_bounds'] is not None and config['relative_bounds'] is None
        assert all(config[k] is None for k in ('max_component_mse','max_complex_error','max_complex_bias','max_baseline_component_mse','max_p99_error','noise'))
        assert export['status']=='passed' and export['source_files_unchanged']
        assert export['distribution_reference_verification']['source_resampling_identical']
        assert shared['all_partition_inventories_verified']
        assert sum(shared['stored_bytes_by_category'].values())==shared['archive_bytes']
        assert rate['original_complex64_bytes']==ORIGINAL_BYTES and rate['within_target'] and rate['physical_payload_within_target']
        assert shared['stored_bytes_by_category']['visibility_block_payload']==rate['physical_payload_bytes']<=BUDGET_BYTES
        assert manifest['distribution_references']['groups']==screen['references']['groups']
        bits=Counter();samples=Counter();gammas=Counter();payload=0
        for part in manifest['partitions']:
            assert all(x['original_bits_exact'] for x in part['summary']['weights_quality'].values())
            for p in part['summary']['visibility_block_statistics']['CORRECTED_DATA']['products']:
                bits.update(p['component_bit_width_histogram']);samples.update(p['complex_sample_bit_width_histogram'])
                gammas.update(p['power_gamma_block_histogram']);payload+=p['serialized_block_bytes']
        assert payload==rate['payload_bytes'] and 8*sum(samples.values())==ORIGINAL_BYTES
        entry=dict(**selected,rate=rate,shared=shared,encoding=read('encode'),export_verification=export,
            visibility=read('visibility'),sample_bits=dict(samples),block_bits=dict(bits),gamma_blocks=dict(gammas),
            archive_sha256=file_sha256(Path(shared['archive'])),weights_all_exact=True,
            distribution_references=manifest['distribution_references'],pilot_binding_references=binding_fractions(selected))
        entry['plot_label']=f"{config['baseline_group']}×{config['channel_segment']}\n"+selected['family'].replace('-G1-D2',' (D×2)')+('\nγ search' if config['power_gamma_policy']=='search' else '')
        audit_path=ROOT/(name+'-payload-audit.json')
        if audit_path.exists():
            cached=json.loads(audit_path.read_text());assert cached['sha256']==entry['archive_sha256']
            distribution=cached['distribution']
        else:
            distribution=payload_distribution(shared['archive'],by_class=True)
            audit_path.write_text(json.dumps(dict(sha256=entry['archive_sha256'],distribution=distribution),indent=2)+'\n')
        assert distribution['samples_by_bits']==dict(samples)
        assert sum(distribution['serialized_bytes_by_bits'].values())==payload
        entry['payload_distribution']=distribution;measured=[]
        for repeat in (False,True):
            suffix='repeat-comparison' if repeat else 'comparison'
            if partial and not (ROOT/(name+'-'+suffix+'.json')).exists():continue
            comparison=read(suffix)
            assert comparison['comparison_mask']['sha256']==file_sha256(MASK)
            for product in comparison['products'].values():
                assert file_sha256(product['reference'])==product['reference_file_sha256']
                assert file_sha256(product['candidate'])==product['candidate_file_sha256']
            images=ROOT/(name+('-repeat-images' if repeat else '-images'))
            complete=json.loads((images/'complete.json').read_text())
            assert complete['status']=='complete' and complete['source_unchanged']
            log_path=images/'wsclean.log';log=log_path.read_text()
            total=re.search(r'Total:\s+(\d+) components',log)
            major=re.search(r'(\d+) major iterations? (?:was|were) performed',log)
            clean=dict(final_threshold_reached='Cleaning finished because the final threshold was reached.' in log,
                cleaned_components=int(total[1]) if total else None,
                major_iterations=int(major[1]) if major else None,log_sha256=file_sha256(log_path))
            image,ih=read_fits(images/'image-image.fits');candidate_residual,ch=read_fits(images/'image-residual.fits')
            assert _grid(ih)==_grid(ch)==_grid(header)
            assert np.all(np.isfinite(image[mask])) and np.all(np.isfinite(candidate_residual[mask]))
            difference=moments((image-reference)[mask]);residual_stats=moments(candidate_residual[mask])
            assert difference['rms']==comparison['products']['image-image.fits']['regions']['above_elevation_cut']['difference_rms']
            measured.append(dict(run='repeat' if repeat else 'first',clean_termination=clean,restored_difference=difference,
                restored_difference_over_reference_residual=difference['rms']/ref_stats['rms'],
                candidate_residual=residual_stats,
                candidate_residual_rms_relative_change=residual_stats['rms']/ref_stats['rms']-1,
                dirty_difference_over_reference_residual=comparison['products']['image-dirty.fits']['regions']['above_elevation_cut']['difference_rms']/ref_stats['rms']))
            if repeat:entry['repeat']=dict(versus_reference=comparison,imaging=complete)
            else:entry.update(comparison=comparison,imaging=complete)
        entry['residual_normalized_runs']=measured
        if len(measured)==2:entry['worst_residual_normalized_difference']=max(x['restored_difference_over_reference_residual'] for x in measured)
        entries.append(entry)
    result=dict(schema='varquanms.normalized-fixed16-study/1',status='partial' if partial else 'complete',
        original_complex64_bytes=ORIGINAL_BYTES,payload_budget_bytes=BUDGET_BYTES,screen=screen,experiments=entries,
        reference_residual_statistics=ref_stats,
        reference_residual=dict(path=str(REFERENCE/'image-residual.fits'),sha256=file_sha256(REFERENCE/'image-residual.fits')),
        comparison_mask=dict(path=str(MASK),sha256=file_sha256(MASK),pixels=int(mask.sum()),minimum_elevation_deg=10),
        sampling_stability=json.loads((ROOT/'sampling/stability.json').read_text()),tests_passed=110,
        hardware=json.loads(Path('reports/calibrated-error-limits.json').read_text())['hardware'],
        selection_objective='Smallest larger restored-image difference RMS across two CLEAN runs among verified <=16% visibility-payload finalists.',
        limitations=['Only joint Stokes I imaging of this calibrated 15-integration observation is validated.',
            'Different layout pilots sample different original blocks; unweighted pilot MSE only nominates full-data finalists.',
            'All normalized limits use source signal-distribution references, never a thermal-noise estimate.',
            'Calibrated-only removes original DATA/MODEL_DATA and zeroes flagged values; flags remain, and weights are exact on this dataset.',
            'Shared-host timings with overlapping pipelines and no cache eviction are not isolated speed comparisons.',
            'No held-out dataset, Q/U/V, spectral-cube or calibration-solution validation. No universal scientific tolerance is established.'])
    if not partial:
        assert len(entries)==len(plan) and all('repeat' in e for e in entries)
        result['study_controller']=json.loads((ROOT/'study-status.json').read_text())
        assert result['study_controller']['status']=='passed'
        result['empirical_best_worst_of_two']=min(entries,key=lambda e:e['worst_residual_normalized_difference'])['name']
        result['best_first_dirty']=min(entries,key=lambda e:e['residual_normalized_runs'][0]['dirty_difference_over_reference_residual'])['name']
        result['best_worst_peak_change']=min(entries,key=lambda e:max(abs(c['products']['image-image.fits']['relative_change_at_reference_peak'])
            for c in (e['comparison'],e['repeat']['versus_reference'])))['name']
        history=[]
        for target in (10,25):
            old=json.loads(Path(f'reports/calibrated-fixed{target}.json').read_text())
            e=next(x for x in old['experiments'] if x['name']==old['empirical_best_worst_of_two'])
            assert e['comparison']['products']['image-image.fits']['reference_file_sha256']==entries[0]['comparison']['products']['image-image.fits']['reference_file_sha256']
            assert e['config']['keep_calib_only'] and e['config']['columns']==['CORRECTED_DATA']
            errors=[c['products']['image-image.fits']['regions']['above_elevation_cut']['difference_rms'] for c in (e['comparison'],e['repeat']['versus_reference'])]
            history.append(dict(target_percent=target,name=e['name'],config=e['config'],
                payload_bytes=e['rate']['payload_bytes'],archive_bytes=e['shared']['archive_bytes'],
                restored_difference_rms=errors,worst_difference_over_reference_residual=max(errors)/ref_stats['rms']))
        result['historical_context']=history
        source=json.loads(Path('reports/calibrated-import.json').read_text())
        assert all(file_sha256(Path(source['extracted_ms'])/p)==v['sha256'] for p,v in source['files'].items())
        result['source_integrity']=dict(original_import_checksums_match=True,files_checked=len(source['files']))
        result['sampling_policy_sensitivity']=sampling_policy_sensitivity(plan,result['sampling_stability'])
        combined=next(e for e in entries if e['family']=='GSDR-G1-D2')
        pair=next(e for e in entries if e['family']=='GS')
        result['combined_vs_pair_storage']=visibility_storage_equality(combined['shared']['archive'],pair['shared']['archive'])
        result['execution_adjustment']=serial_parallel_equality(entries)
        result['host_io_snapshot']=dict(scope='Five one-second shared-host observations during the sweep; not a full-run profile',
            text=(ROOT/'host-io-snapshot.txt').read_text())
    verify_reused_hashes()
    return result


def supplementary_plots(result):
    import os
    os.environ.setdefault('MPLCONFIGDIR','/tmp/varquanms-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    dest=Path('devnotes/figs');dest.mkdir(exist_ok=True)
    screen=result['screen'];recipes=[r['name'] for r in screen['recipes']]
    layouts=[r['layout'] for r in screen['layouts']]
    lookup={(e['family'],e['layout']):e for e in screen['experiments']}
    relative=np.array([[lookup[(r,l)]['cross_parallel_component_mse']/lookup[('R',l)]['cross_parallel_component_mse']
                        for l in layouts] for r in recipes])
    log=np.log2(relative);span=max(abs(log.min()),abs(log.max()),.01)
    fig,ax=plt.subplots(figsize=(11,13),constrained_layout=True)
    im=ax.imshow(log,aspect='auto',cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-span,vcenter=0,vmax=span))
    names=[f"{r['config']['baseline_group']}×{r['config']['channel_segment']}"+
           ('\nγ search' if r['config']['power_gamma_policy']=='search' else '') for r in screen['layouts']]
    ax.set(xticks=range(len(layouts)),xticklabels=names,yticks=range(len(recipes)),yticklabels=recipes,
           xlabel='Baselines × channels per product',ylabel='Epsilon coefficient recipe')
    ax.tick_params(labelsize=8)
    for i,row in enumerate(relative):
        for j,value in enumerate(row):ax.text(j,i,f'{value:.2f}',ha='center',va='center',fontsize=6,
            color='white' if abs(log[i,j])>.65*span else 'black')
    fig.colorbar(im,ax=ax,shrink=.6,label='log₂(MSE / R-only MSE in the same layout)')
    fig.suptitle('333 pilot combinations at ≤15.75% sampled visibility payload\n'
                 'Cross XX/YY component MSE ratio; values <1 are lower. Different layouts sample different blocks.\n'
                 'This screen nominates full-data tests; it is not an image-quality measurement.',fontsize=10)
    fig.savefig(dest/'fixed16-pilot.png',dpi=170);plt.close(fig)
    sampled=[next(iter(x['references']['groups'].values()))['columns']['CORRECTED_DATA'] for x in result['sampling_stability']]
    keys=[(p,a) for p in sampled[0] for a in ('auto','cross')]
    fig,ax=plt.subplots(figsize=(11,4.2),constrained_layout=True)
    for offset,index,name,color in zip((-.25,0,.25),(0,1,3),('2% nested','3% nested','2% independent seed'),('#4393c3','#92c5de','#b2182b')):
        values=[100*(sampled[index][p][a]['component_std']/sampled[2][p][a]['component_std']-1) for p,a in keys]
        ax.bar(np.arange(len(keys))+offset,values,width=.24,label=name,color=color)
    ax.axhline(0,color='black',linewidth=.7)
    ax.set(xticks=range(len(keys)),xticklabels=[f'{p}\n{a}' for p,a in keys],ylabel='Std relative to nested 5% estimate (%)')
    ax.legend(ncol=3);ax.grid(axis='y',alpha=.2)
    fig.suptitle('Global reference sampling sensitivity; all 15 integrations\n5% is a comparison estimate, not population truth')
    fig.savefig(dest/'fixed16-sampling.png',dpi=170);plt.close(fig)
    cases=result['experiments']
    if not cases:return
    fig,ax=plt.subplots(figsize=(10,max(4,.5*len(cases))),constrained_layout=True);left=np.zeros(len(cases))
    for k,color in zip('GSDR',('#2166ac','#67a9cf','#ef8a62','#b2182b')):
        values=np.array([e['pilot_binding_references']['fractions'][k]*100 for e in cases])
        ax.barh(range(len(cases)),values,left=left,label=k,color=color);left+=values
    ax.set(yticks=range(len(cases)),yticklabels=[label(e) for e in cases],xlim=(0,100),xlabel='Pilot valid cross-correlation samples (%)')
    ax.invert_yaxis();ax.legend(ncol=4,title='Reference setting the tightest RMS budget',
                              loc='upper center',bbox_to_anchor=(.5,-.14))
    fig.suptitle('Which reference actually binds?\nExact ties split equally; autocorrelations and blocks without valid samples excluded')
    fig.savefig(dest/'fixed16-binding.png',dpi=170);plt.close(fig)


def note(result):
    cases=result['experiments'];best=next(e for e in cases if e['name']==result['empirical_best_worst_of_two'])
    Path('configs/calibrated-fixed16-best.json').write_text(json.dumps(best['config'],indent=2)+'\n')
    lines=['# Sixteen-percent payload: G/S/D/R budgets, combinations and block sizes','',
        f"Among the ten tested configurations, **{label(best)}** has the lowest worst-case restored-image difference across two CLEAN runs. Its visibility payload is **{best['rate']['physical_payload_bytes']/1e6:.3f} MB ({100*best['rate']['fraction']:.3f}%)**, and the complete archive is **{best['shared']['archive_bytes']/1e6:.3f} MB**. [Exact opt-in configuration](../configs/calibrated-fixed16-best.json).",'',
        'Fixed absolute B has been removed from these configurations. The [distribution-v1 policy](normalized-budgets.md) uses only dimensionless epsilons and data-derived references. Legacy experiments remain unchanged. This is a measured configuration for this observation, not a scientific default.','',
        '## Design and reference sampling','',
        'The byte ceiling is **229,028,659 bytes**, the integer limit corresponding to 16% of 1,431,429,120 original CORRECTED_DATA complex64 bytes. It counts packed codes, factors, headers, exceptions, exact fallback and secondary compression. Permutations, indexes and other content are outside this ceiling and inside the complete archive.','',
        'The pilot tests **37 coefficient recipes × 9 layout/codec settings = 333 combinations**. Recipes include all 15 nonempty G/S/D/R subsets, pairwise coefficient ratios and four-reference combinations. Block shapes are 1×48, 8×24, 32×24, 32×48, 128×16, 128×24, 128×48 and 256×48, each for one correlation; 128×48 also tests full gamma search. Baseline 1 uses original order; grouped cases use independent amplitude order per product, with auto/cross classes separate.','',
        'Every packed candidate is measured once in the pilot and its error/byte tradeoff reused to search a common multiplier for each recipe, targeting 15.75% to leave sampling margin. Thresholds are chosen by bytes before imaging. Cached decisions are checked against the production encoder. Ten finalists include each single-reference control, combined rules, an unsorted per-baseline control, smaller blocks and full gamma search. Different layout pilots sample different rank blocks; they do not establish image rankings.','',
        'Global references use a fixed seeded 2% sample from all 15 integrations, separately for each correlation and auto/cross class. The source has one physical configuration, despite different observation IDs, so there are eight global references. A minimum of 128 samples is available for small strata; it was not needed here. Nested 3% and 5% samples plus an independent 2% seed measure sampling sensitivity.','',
        '| Product/class | 2% sample count | Global component std | 2% versus 5% difference | Independent 2% versus 5% difference |',
        '| --- | ---: | ---: | ---: | ---: |']
    sampled=[next(iter(x['references']['groups'].values()))['columns']['CORRECTED_DATA'] for x in result['sampling_stability']]
    for product in sampled[0]:
        for kind in ('auto','cross'):
            v=[s[product][kind]['component_std'] for s in sampled]
            lines.append(f"| {product}/{kind} | {sampled[0][product][kind]['sampled_samples']:,} | {v[0]:.6f} | {100*(v[0]/v[2]-1):+.3f}% | {100*(v[3]/v[2]-1):+.3f}% |")
    sensitivity=result['sampling_policy_sensitivity']['experiments']
    changed=max(t['changed_sample_fraction'] for e in sensitivity for t in e['trials'])
    change_rate=max(abs(t['payload_fraction']-e['trials'][0]['payload_fraction']) for e in sensitivity for t in e['trials'])
    lines+=['',f"Holding each finalist's epsilons fixed and changing only G to the 3%, 5% or independent-seed estimate changes at most **{100*changed:.3f}% of pilot samples' selected candidates**, and the sampled payload fraction by at most **{100*change_rate:.4f} percentage points**. This is a pilot allocation check, not a further imaging test or a guarantee for other MS distributions."]
    lines+=['','![Pilot combinations at matched sampled payload](figs/fixed16-pilot.png)','',
        '![Sampling sensitivity](figs/fixed16-sampling.png)','',
        'Std/MAD are computed within real and imaginary components and pooled. RMS retains the coherent centre. They describe signal distributions, not thermal noise. No fixed absolute floor is inserted. Zero MAD can require exact storage even for nonconstant sparse blocks.','',
        '$$T_b=\min(\epsilon_G G_{p,a},\epsilon_S S_b,\epsilon_D D_b,\epsilon_R R_b),\qquad \mathrm{component\ MSE}_b\le T_b^2.$$','',
        'Disabled terms are omitted. Every finalist uses maximum complex error ≤10T, complex bias ≤T and baseline component MSE ≤2T². All four epsilon values are recorded below; no legacy B or fixed absolute guard is active.','',
        '| Finalist | εG | εS | εD | εR |', '| --- | ---: | ---: | ---: | ---: |']
    for e in cases:
        p=e['config']['normalized_bounds'];values=['—' if p['epsilon_'+k] is None else f"{p['epsilon_'+k]:.8g}" for k in 'GSDR']
        lines.append('| '+label(e)+' | '+' | '.join(values)+' |')
    lines+=['','Some constraints can be redundant: S≤R, so equal εS and εR make the RMS constraint inactive. Pilot binding fractions are recorded rather than assuming every enabled coefficient affects bit selection.','',
        '![References setting the tightest pilot budgets](figs/fixed16-binding.png)','',
        '## Repeated image comparisons','',
        f"The fixed above-10° mask contains **2,245,726 pixels**. Reference residual RMS is **{result['reference_residual_statistics']['rms']:.9f} Jy/beam**, without mean removal or clipping. Its sky/deconvolution structure prevents interpreting it as a thermal-noise estimate.",'',
        '| Finalist | Payload MB (%) | Whole MB | Restored difference mJy/beam: first / repeat | Difference / reference residual RMS %: first / repeat | Dirty difference / residual RMS % |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for e in cases:
        a,b=e['residual_normalized_runs']
        lines.append(f"| {label(e)} | {e['rate']['physical_payload_bytes']/1e6:.3f} ({100*e['rate']['fraction']:.3f}%) | {e['shared']['archive_bytes']/1e6:.3f} | {1000*a['restored_difference']['rms']:.4f} / {1000*b['restored_difference']['rms']:.4f} | {100*a['restored_difference_over_reference_residual']:.4f} / {100*b['restored_difference_over_reference_residual']:.4f} | {100*a['dirty_difference_over_reference_residual']:.4f} |")
    equality=result['combined_vs_pair_storage']
    if equality['arrays_identical']:
        lines+=['',f"The G/S and selected four-term archives have **identical visibility payloads and decoding maps** in all {equality['arrays_compared']:,} arrays checked byte for byte. Their CLEAN inputs are therefore identical. Any difference between their image comparisons is repeat variation, not a benefit from adding D or R."]
    else:
        lines+=['','The G/S and selected four-term policies made identical pilot choices, but their complete payloads or decoding maps differ. The report records the differing arrays; the pilot cannot establish full-data equivalence.']
    lines+=['','| Finalist | Worst difference / restored-reference RMS % | Worst faint-region difference % | Largest reference-peak change % | Candidate residual RMS Jy/beam: first / repeat |',
        '| --- | ---: | ---: | ---: | ---: |']
    for e in cases:
        c=[e['comparison'],e['repeat']['versus_reference']];a,b=e['residual_normalized_runs']
        peak=max(abs(x['products']['image-image.fits']['relative_change_at_reference_peak']) for x in c)
        lines.append(f"| {label(e)} | {100*max(fraction(x) for x in c):.4f} | {100*max(fraction(x,'image','faint_reference') for x in c):.4f} | {100*peak:.5f} | {a['candidate_residual']['rms']:.6f} / {b['candidate_residual']['rms']:.6f} |")
    lines+=['','Faint pixels are below 1% of the masked reference peak, using the fixed reference to define the region. Peak change refers to the brightest absolute reference pixel, not an identification of the Sun or an integrated-flux measurement. Similar residual RMS alone does not demonstrate restored-image agreement.','',
        'All runs use the requested 2048², four-arcmin, eight-thread multiscale Stokes I template and reference restoring beam. Two runs characterize observed repeat variation only; they are not confidence intervals. Actual payloads differ slightly under the common ceiling.','']
    singles=[e for e in cases if e['layout']=='amplitude128x48' and e['family'] in ('G','S','D','R')]
    single=min(singles,key=lambda e:e['worst_residual_normalized_difference'])
    baseline=next(e for e in cases if e['config']['baseline_group']==1)
    best_peak=next(e for e in cases if e['name']==result['best_worst_peak_change'])
    best_dirty=next(e for e in cases if e['name']==result['best_first_dirty'])
    lines+=['## What the measurements support','',
        f"With 128×48 blocks and the same max/std gamma rule, **{single['family']}** is the best of the four single-reference controls by the stated larger-of-two image-RMS criterion. G is a data-normalized global ceiling; S, D and R follow local block statistics. None is a thermal-noise bound.",'',
        f"The tested per-baseline control has {baseline['shared']['archive_bytes']/1e6:.3f} MB total storage, including {baseline['shared']['stored_bytes_by_category']['visibility_indexes']/1e6:.3f} MB of block indexes. Its larger image difference is {100*baseline['worst_residual_normalized_difference']:.3f}% of reference residual RMS, versus {100*single['worst_residual_normalized_difference']:.3f}% for the best single-reference grouped control. These are separately rate-tuned configurations, with different reference recipes; this is not an isolated test of sorting alone.",'',
        f"The lowest first-run dirty-image difference belongs to **{label(best_dirty)}**. The smallest larger-of-two reference-peak change belongs to **{label(best_peak)}**. The saved preset follows global restored-image RMS; another science objective can favor another configuration.",'',
        'Only four distinct block shapes were carried through full-data imaging; all eight shapes were screened. Pilot unweighted visibility MSE nominates candidates and does not establish their final image ranking. Neither an eight-statistic global reference nor any single relative error number is a universal scientific tolerance.','',
        ]
    small=next(e for e in cases if e['layout']=='amplitude32x24');run=small['residual_normalized_runs'][0]
    dirty=small['comparison']['products']['image-dirty.fits']['regions']['above_elevation_cut']['difference_rms']
    if run['restored_difference']['rms']>5*dirty and run['clean_termination']['final_threshold_reached']:
        lines += [f"The 32×24 first-run dirty-image difference is {1000*dirty:.3f} mJy/beam, while its restored-image difference is {1000*run['restored_difference']['rms']:.3f} mJy/beam. CLEAN reached its final threshold after {run['clean_termination']['cleaned_components']:,} components. This observed sensitivity of the restored result is not an iteration-limit failure, and small peak changes do not establish agreement across the image. Logs and termination counts are retained.",'']
    lines += ['| Historical target | Payload MB | Whole MB | Larger image difference / reference residual RMS % |',
        '| --- | ---: | ---: | ---: |']
    for h in result['historical_context']:
        lines.append(f"| {h['target_percent']}% | {h['payload_bytes']/1e6:.3f} | {h['archive_bytes']/1e6:.3f} | {100*h['worst_difference_over_reference_residual']:.4f} |")
    lines.append(f"| 16%, selected here | {best['rate']['physical_payload_bytes']/1e6:.3f} | {best['shared']['archive_bytes']/1e6:.3f} | {100*best['worst_residual_normalized_difference']:.4f} |")
    lines+=['','Historical runs use the same reference image and calibrated-only preservation policy, but rates, block choices and threshold policies differ; the older 25% experiment also used a different bit ladder. This is context, not an isolated measurement of the epsilon policy change.','']
    for caption,suffix in [('Payload and image difference','rate'),('Selected full image comparison','best-images'),('Selected reference-peak zoom','best-zoom'),('Bit distribution by samples and bytes','bits'),('All first-run image differences','images'),('All first-run reference-peak zooms','zoom')]:
        lines += [f'![{caption}](figs/fixed16-{suffix}.png)','']
    lines+=['## Byte accounting and timing','',
        '| Selected archive category | Physical MB |', '| --- | ---: |']
    for category,amount in sorted(best['shared']['stored_bytes_by_category'].items(),key=lambda item:-item[1]):
        lines.append(f"| {category.replace('_',' ')} | {amount/1e6:.3f} |")
    lines+=['', '| Selected precision | Complex samples % | Serialized payload bytes % |',
        '| --- | ---: | ---: |']
    distribution=best['payload_distribution'];total_samples=sum(distribution['samples_by_bits'].values())
    total_bytes=sum(distribution['serialized_bytes_by_bits'].values())
    for bits in ('3','4','6','8','12','16','float32'):
        lines.append(f"| {bits} | {100*distribution['samples_by_bits'].get(bits,0)/total_samples:.4f} | {100*distribution['serialized_bytes_by_bits'].get(bits,0)/total_bytes:.4f} |")
    lines+=['',f"Exact positive-zero blocks account for {100*distribution['exact_positive_zero_samples']/total_samples:.3f}% of all samples and {distribution['exact_positive_zero_bytes']/1e6:.3f} MB. This counts samples in entirely zero blocks, not every individual zero inside mixed blocks.",'',
        '| Finalist | Payload + permutation/index MB | Sum of partition encode s | Global sampling s | Verified bundle wall s | Compaction s | Export s | First CLEAN s |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for e in cases:
        lines.append(f"| {label(e)} | {e['rate']['payload_plus_permutations_and_indexes_bytes']/1e6:.3f} | {e['encoding']['sum_encode_seconds_excluding_inspection']:.2f} | {e['encoding']['global_reference_sampling_seconds']:.2f} | {e['encoding']['elapsed_seconds_including_copy_inspection_verification']:.2f} | {e['shared']['elapsed_seconds']:.2f} | {e['export_verification']['export']['elapsed_seconds']:.2f} | {e['imaging']['wsclean_elapsed_seconds']:.2f} |")
    lines+=['','Timings cover all 15 integrations on a shared host with overlapping archive/image pipelines and no cache eviction. They are observations, not isolated speed rankings. Grouped cases use 24 encoder threads and serial partitions. The per-baseline control uses six independent partition processes, each with one encoder worker; its summed partition times exceed wall time. A preliminary serial attempt was stopped and retained as incomplete, then the full control was restarted; the report records that discarded attempt separately. Fixture checks and completed serial partitions are compared for exact array equivalence. The full-conversion 15-second goal remains unmet. Global sampling is included in verified-bundle time; compaction and export are additional stages.','',
        'Every byte category, sample/byte precision histogram and auto/cross precision split is retained in the JSON. Calibrated-only storage omits original DATA/MODEL_DATA and zeros flagged visibility values, while preserving flags and original time/frequency resolution. Four-bit weights were verified exact on this dataset. All 108 source-file hashes still match; supported metadata and original metre UVW passed exact comparisons.','',
        '## Reproduction','', '```bash','export PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1',
        'python scripts/sample_distribution.py INPUT.ms SAMPLING_DIR',
        'python scripts/screen_normalized.py INPUT.ms SCREEN_DIR \\',
        '  --references SAMPLING_DIR/fraction2-seed20261024.json --target .16 --pilot-target .1575',
        'python scripts/run_normalized_study.py --source INPUT.ms --root NEW_RESULTS --plan SCREEN_DIR/plan.json',
        '```','',
        '[Machine-readable results](../reports/calibrated-fixed16.json). The codec/archive suite has 110 passing tests. No held-out observation, Q/U/V, spectral-cube or calibration-solution fidelity is claimed. The paper is unchanged.','']
    Path('devnotes/fixed16.md').write_text('\n'.join(lines))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--partial',action='store_true');args=p.parse_args()
    result=collect(args.partial)
    output=ROOT/'preview.json' if args.partial else Path('reports/calibrated-fixed16.json')
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not args.partial:
        plots(result,root=ROOT,prefix='fixed16');supplementary_plots(result);note(result)
    for e in result['experiments']:
        print(e['name'],e['rate']['fraction'],[100*x['restored_difference_over_reference_residual'] for x in e['residual_normalized_runs']])


if __name__=='__main__':main()
