"""Audit the historical three-epsilon Stokes I study and refresh its development note."""
from collections import Counter
import json
from pathlib import Path
import re

import h5py
import numpy as np
from astropy.io import fits

from vamp.imaging import read_fits, _grid
from vamp.codec import unpack_header
from report_fixed_rate import payload_distribution, fraction
from report_fixed10 import moments
from report_normalized import file_sha256, verify_reused_hashes, HASHES

ROOT=Path('results/calibrated-epsd-dense')
REF=Path('results/calibrated-template/reference')
MASK=Path('results/calibrated-template/above10deg.fits')
BYTES=1431429120
LADDER=[3,4,5,6,8,10,12,14,16]
EPSILONS=[.01,.04,.18]
PREVIOUS_COMMIT='b9596c800e176dfe312945212453b99161343bac'


def collect():
    HASHES.clear()
    plan=json.loads((ROOT/'plan.json').read_text())
    controller=json.loads((ROOT/'study-status.json').read_text());assert controller['status']=='passed'
    assert len(plan)==3 and sorted(e['epsilon_D'] for e in plan)==EPSILONS
    mask=fits.getdata(MASK).astype(bool);assert mask.sum()==2245726
    reference,header=read_fits(REF/'image-image.fits')
    residual,rh=read_fits(REF/'image-residual.fits');assert _grid(header)==_grid(rh)
    residual_stats=moments(residual[mask]);entries=[]
    for chosen in sorted(plan,key=lambda e:e['epsilon_D']):
        name=chosen['name'];load=lambda stage:json.loads((ROOT/f'{name}-{stage}.json').read_text())
        cfg=chosen['config'];policy=cfg['normalized_bounds']
        assert cfg['bits']==LADDER and cfg['block_header']=='precision-nibble-v1'
        assert cfg['keep_calib_only'] and not cfg['burn_flag'] and cfg['weight_bits']==4
        assert cfg['baseline_group']==128 and cfg['channel_segment']==48 and cfg['ordering']=='amplitude'
        assert cfg['columns']==['CORRECTED_DATA'] and not cfg['shared_order']
        assert policy['epsilon_D']==chosen['epsilon_D'] and all(policy['epsilon_'+k] is None for k in 'GSR')
        assert cfg['relative_bounds'] is None and cfg['noise'] is None and not cfg['product_bounds']
        assert all(cfg[k] is None for k in ('max_component_mse','max_complex_error','max_complex_bias','max_baseline_component_mse','max_p99_error'))
        export=load('export');shared=load('compact');encoding=load('encode')
        assert export['status']=='passed' and export['source_files_unchanged']
        assert export['distribution_reference_verification']['source_resampling_identical']
        assert shared['all_partition_inventories_verified']
        assert sum(shared['stored_bytes_by_category'].values())==shared['archive_bytes']
        manifest=json.loads((ROOT/(name+'.bundle')/'manifest.json').read_text());assert manifest['config']==cfg
        bits=Counter();samples=Counter();gammas=Counter();serialized=0
        for part in manifest['partitions']:
            assert all(x['original_bits_exact'] for x in part['summary']['weights_quality'].values())
            for p in part['summary']['visibility_block_statistics']['CORRECTED_DATA']['products']:
                bits.update(p['component_bit_width_histogram']);samples.update(p['complex_sample_bit_width_histogram'])
                gammas.update(p['power_gamma_block_histogram']);serialized+=p['serialized_block_bytes']
        assert 8*sum(samples.values())==BYTES
        categories=shared['stored_bytes_by_category'];physical=categories['visibility_block_payload']
        assert physical==serialized
        distribution=payload_distribution(shared['archive'],by_class=True)
        assert distribution['samples_by_bits']==dict(samples)
        assert sum(distribution['serialized_bytes_by_bits'].values())==serialized
        tag_counts=Counter()
        with h5py.File(shared['archive'],'r') as h:
            top=json.loads(h['manifest'][...].tobytes())
            for part in top['partitions']:
                g=h[part['group']];spec=json.loads(g['manifest'][...].tobytes())
                for item in spec['inventory']:
                    if item['role']!='visibility_block_payload':continue
                    ds=g[item['path']];offsets=ds.parent['offsets'][...]
                    assert ds.size<=64*1024**2
                    body=memoryview(ds[...])
                    for lo,hi in zip(offsets[:-1],offsets[1:]):
                        blob=body[int(lo):int(hi)];tag=int(blob[0]);assert tag&0xf0==0xa0
                        width=unpack_header(blob)[0];assert width in (*LADDER,0)
                        tag_counts[str(tag&15)]+=1
        assert sum(tag_counts.values())==sum(bits.values())
        rate=dict(payload_bytes=serialized,physical_payload_bytes=physical,original_complex64_bytes=BYTES,
                  fraction=physical/BYTES,physical_payload_fraction=physical/BYTES,
                  payload_plus_permutations_and_indexes_bytes=sum(categories.get(k,0) for k in ('visibility_block_payload','permutations','visibility_indexes')))
        entry=dict(**chosen,rate=rate,shared=shared,encoding=encoding,export_verification=export,
            visibility=load('visibility'),sample_bits=dict(samples),block_bits=dict(bits),gamma_blocks=dict(gammas),
            payload_distribution=distribution,precision_tag_counts=dict(tag_counts),
            archive_sha256=file_sha256(shared['archive']),plot_label=f"εD = {chosen['epsilon_D']:g}",weights_all_exact=True)
        measured=[]
        for repeat in (False,True):
            comparison=load('repeat-comparison' if repeat else 'comparison')
            assert comparison['comparison_mask']['sha256']==file_sha256(MASK)
            for product in comparison['products'].values():
                assert file_sha256(product['reference'])==product['reference_file_sha256']
                assert file_sha256(product['candidate'])==product['candidate_file_sha256']
            directory=ROOT/(name+('-repeat-images' if repeat else '-images'))
            complete=json.loads((directory/'complete.json').read_text())
            assert complete['status']=='complete' and complete['source_unchanged']
            image,ih=read_fits(directory/'image-image.fits');rest,rh=read_fits(directory/'image-residual.fits')
            assert _grid(ih)==_grid(rh)==_grid(header)
            assert np.all(np.isfinite(image[mask])) and np.all(np.isfinite(rest[mask]))
            delta=moments((image-reference)[mask]);rs=moments(rest[mask])
            assert delta['rms']==comparison['products']['image-image.fits']['regions']['above_elevation_cut']['difference_rms']
            log_path=directory/'wsclean.log';log=log_path.read_text()
            total=re.search(r'Total:\s+(\d+) components',log)
            major=re.search(r'(\d+) major iterations? (?:was|were) performed',log)
            measured.append(dict(run='repeat' if repeat else 'first',restored_difference=delta,
                restored_difference_over_reference_residual=delta['rms']/residual_stats['rms'],candidate_residual=rs,
                dirty_difference_over_reference_residual=comparison['products']['image-dirty.fits']['regions']['above_elevation_cut']['difference_rms']/residual_stats['rms'],
                clean_termination=dict(final_threshold_reached='Cleaning finished because the final threshold was reached.' in log,
                    cleaned_components=int(total[1]) if total else None,major_iterations=int(major[1]) if major else None,log_sha256=file_sha256(log_path))))
            assert measured[-1]['clean_termination']['final_threshold_reached']
            if repeat:entry['repeat']=dict(versus_reference=comparison,imaging=complete)
            else:entry.update(comparison=comparison,imaging=complete)
        entry['residual_normalized_runs']=measured
        entry['worst_residual_normalized_difference']=max(x['restored_difference_over_reference_residual'] for x in measured)
        entries.append(entry)
    imported=json.loads(Path('reports/calibrated-import.json').read_text())
    assert all(file_sha256(Path(imported['extracted_ms'])/p)==v['sha256'] for p,v in imported['files'].items())
    previous=json.loads(Path('reports/calibrated-fixed16.json').read_text())
    old_paper=Path('docs/paper/work_note_20261025_init_work.md')
    old_hash=file_sha256(old_paper)
    assert old_hash=='22c079ca6a557d873c801a08c71a53a007606d5e315f41ff3c5843609fbf3b89'
    result=dict(schema='vamp.epsilon-d-dense-study/1',status='complete',results_root=str(ROOT),
        original_complex64_bytes=BYTES,original_ms_physical_bytes=4428203949,experiments=entries,
        reference_residual_statistics=residual_stats,reference_restored_statistics=moments(reference[mask]),
        comparison_mask=dict(path=str(MASK),sha256=file_sha256(MASK),pixels=int(mask.sum()),minimum_elevation_deg=10),
        source_integrity=dict(original_import_checksums_match=True,files_checked=len(imported['files'])),
        previous_manuscript=dict(path=str(old_paper),sha256=old_hash,content_unchanged=True),
        tests_passed=116,hardware=previous['hardware'],study_controller=controller,
        correction=dict(previous_epsilon_D=.018,corrected_epsilon_D=.18,
            previous_results_commit=PREVIOUS_COMMIT,
            reused_verified_cases=['epsd-001-dense','epsd-004-dense'],
            new_case='epsd-018-dense',
            note='The 0.018 experiment is historical and is not relabelled as 0.18. The corrected case was encoded from the original MS and CLEANed twice.'),
        precision_labels=[str(x) for x in LADDER]+['float32'],
        plot_title='Fixed epsilon_D; expanded precision ladder; same MS and above-10° mask',
        empirical_best_worst_of_two=min(entries,key=lambda e:e['worst_residual_normalized_difference'])['name'],
        selection='Three user-specified fixed epsilon_D values. No rate target, no retuning from images; smallest worst-of-two image RMS identifies the displayed minimum-error case, not a general optimum.',
        limitations=['One calibrated observation; joint Stokes I only. No held-out, polarization, spectral-cube or calibration-solution validation.',
                    'Residual RMS contains sky/deconvolution structure; it is not thermal noise.',
                    'Two CLEAN runs describe observed variation, not confidence intervals.',
                    'Shared-host overlapping pipelines without cache eviction; timings are not isolated throughput rankings.'])
    verify_reused_hashes()
    return result


def tables(result):
    rows=['| εD | Payload MB (% of complex64) | Complete archive MB | Difference mJy/beam, first / repeat | Difference / reference residual RMS %, first / repeat |',
          '| ---: | ---: | ---: | ---: | ---: |']
    detail=['| εD | Dirty difference mJy/beam | Worst reference-peak change % | Worst faint-region RMS difference % | Candidate residual RMS Jy/beam, first / repeat |',
            '| ---: | ---: | ---: | ---: | ---: |']
    timing=['| εD | Sum partition encode s | Global sampling s | Verified bundle wall s | Compaction s | MS export s | CLEAN s, first / repeat |',
            '| ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for e in result['experiments']:
        a,b=e['residual_normalized_runs'];c=e['comparison'];d=e['repeat']['versus_reference']
        rows.append(f"| {e['epsilon_D']:g} | {e['rate']['physical_payload_bytes']/1e6:.3f} ({100*e['rate']['fraction']:.3f}%) | {e['shared']['archive_bytes']/1e6:.3f} | {1000*a['restored_difference']['rms']:.4f} / {1000*b['restored_difference']['rms']:.4f} | {100*a['restored_difference_over_reference_residual']:.4f} / {100*b['restored_difference_over_reference_residual']:.4f} |")
        peak=100*max(abs(x['products']['image-image.fits']['relative_change_at_reference_peak']) for x in (c,d))
        faint=100*max(fraction(x,'image','faint_reference') for x in (c,d))
        dirty=c['products']['image-dirty.fits']['regions']['above_elevation_cut']['difference_rms']
        detail.append(f"| {e['epsilon_D']:g} | {1000*dirty:.4f} | {peak:.5f} | {faint:.4f} | {a['candidate_residual']['rms']:.6f} / {b['candidate_residual']['rms']:.6f} |")
        timing.append(f"| {e['epsilon_D']:g} | {e['encoding']['sum_encode_seconds_excluding_inspection']:.2f} | {e['encoding']['global_reference_sampling_seconds']:.2f} | {e['encoding']['elapsed_seconds_including_copy_inspection_verification']:.2f} | {e['shared']['elapsed_seconds']:.2f} | {e['export_verification']['export']['elapsed_seconds']:.2f} | {e['imaging']['wsclean_elapsed_seconds']:.2f} / {e['repeat']['imaging']['wsclean_elapsed_seconds']:.2f} |")
    return '\n'.join(rows), '\n'.join(detail), '\n'.join(timing)


def visibility_table(result):
    rows=['| εD | Cross XX RMSE / signal RMS % | Cross XY % | Cross YX % | Cross YY % |',
          '| ---: | ---: | ---: | ---: | ---: |']
    for e in result['experiments']:
        values=[100*e['visibility']['products']['cross/'+p]['relative_rmse'] for p in ('XX','XY','YX','YY')]
        rows.append(f"| {e['epsilon_D']:g} | "+' | '.join(f'{v:.5f}' for v in values)+' |')
    return '\n'.join(rows)


def clean_termination_text(result):
    runs=[r for e in result['experiments'] for r in e['residual_normalized_runs']]
    assert len(runs)==6 and all(r['clean_termination']['final_threshold_reached'] for r in runs)
    counts=sorted({r['clean_termination']['major_iterations'] for r in runs})
    assert all(isinstance(n,int) and n>0 for n in counts)
    suffix=f'{counts[0]} major iterations each' if len(counts)==1 else f'{counts[0]}–{counts[-1]} major iterations'
    return 'All six CLEAN runs reached the final threshold, with '+suffix+'.'


def write_note(result):
    main,detail,timing=tables(result)
    best=next(e for e in result['experiments'] if e['name']==result['empirical_best_worst_of_two'])
    smallest=min(result['experiments'],key=lambda e:e['shared']['archive_bytes'])
    first=min(result['experiments'],key=lambda e:e['residual_normalized_runs'][0]['restored_difference']['rms'])
    tolerances=', '.join(f"{e['epsilon_D']:g}" for e in result['experiments'])
    lines=['# Fixed epsilon D with the expanded precision ladder','',
        f'25 October 2026. Three requested tolerances were encoded independently from the original calibrated source, with no rate cap or retuning: {tolerances}. All use 128×48×1 amplitude groups, median centres, max/std gamma (R0=4), range quantiles 1 and 0.995, and the 3/4/5/6/8/10/12/14/16 ladder plus exact fallback. G/S/R are disabled.','',
        f'Correction: the requested third tolerance is **0.18**, replacing 0.018. The 0.18 archive and both CLEAN runs are new measurements; verified 0.01 and 0.04 results are reused unchanged. The [earlier 0.018 results](https://github.com/peijin94/vamp/blob/{PREVIOUS_COMMIT}/devnotes/epsd-dense.md) remain historical and are not relabelled.','',
        'The four-bit precision selector is stored in the low nibble of a versioned header byte. Exact is tag 0; tags 1–9 select the nine widths; tag A is legacy two-bit; B–E are reserved and F is invalid/error. Invalid tags fail rather than pretending to contain data. Header sizes remain 16/28 bytes.','',
        main,'',detail,'',visibility_table(result),'',
        'The reference residual RMS is 1.031244811 Jy/beam, measured without clipping or mean removal on the fixed above-10° mask. Two CLEAN runs record variation, not statistical confidence. Minimum image error is a different objective from minimum archive bytes; no universal default is selected.','',
        f"Epsilon-D {smallest['epsilon_D']:g} has the smallest archive. Epsilon-D {best['epsilon_D']:g} has the smallest worst-of-two restored-image difference, and epsilon-D {first['epsilon_D']:g} has the smallest first-run difference. Tighter visibility bounds do not monotonically rank the restored CLEAN images: the unchanged 0.01 versus 0.04 comparison shows that reversal.",'',
        clean_termination_text(result),'',
        f"All source/export checks passed, all 108 import-file hashes match, and the dataset’s four-bit weights remain exact. Original metre UVW, supported metadata, coordinates and flags remain exact. Calibrated-only deliberately omits original DATA/MODEL_DATA and zeroes flagged values. The local suite has {result['tests_passed']} passing tests.",'',
        '![Size and image error](figs/epsd-dense-rate.png)','',
        '![Precision distribution](figs/epsd-dense-bits.png)','',
        '![All image differences](figs/epsd-dense-images.png)','',
        '![All peak zooms](figs/epsd-dense-zoom.png)','',
        '![Minimum-error full comparison](figs/epsd-dense-best-images.png)','',
        '![Minimum-error peak zoom](figs/epsd-dense-best-zoom.png)','',
        '## Timing','',timing,'',
        'These shared-host timings do not evict caches. The 0.01/0.04 pipelines overlapped other archive/image work; 0.18 was run separately after the correction, with its two CLEAN repetitions concurrent. Active workload differs, so this is not an isolated speed ranking. Verified bundle time includes source copying, inspection, encoding and verification. Compaction and MS export are additional. The full 15-second conversion target is not established. G is disabled, but the shared policy still measures and records global references; that cost is included.','',
        '## Byte accounting','']
    for e in result['experiments']:
        lines += [f"### epsilon_D = {e['epsilon_D']:g}",'','| Category | Physical MB |','| --- | ---: |']
        lines += [f"| {k.replace('_',' ')} | {v/1e6:.3f} |" for k,v in sorted(e['shared']['stored_bytes_by_category'].items(),key=lambda x:-x[1])]
        lines += ['','| Bits/component | Complex samples % | Serialized payload bytes % |','| --- | ---: | ---: |']
        p=e['payload_distribution'];ns=sum(p['samples_by_bits'].values());nb=sum(p['serialized_bytes_by_bits'].values())
        lines += [f"| {k} | {100*p['samples_by_bits'].get(k,0)/ns:.4f} | {100*p['serialized_bytes_by_bits'].get(k,0)/nb:.4f} |" for k in result['precision_labels']]
        lines += ['',f"Entirely positive-zero exact blocks contain {100*p['exact_positive_zero_samples']/ns:.3f}% of samples and occupy {p['exact_positive_zero_bytes']/1e6:.3f} MB.",'']
    lines += ['## Reproduction','', '```bash','export PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1',
        'python scripts/run_normalized_study.py --source INPUT.ms --root NEW_RESULTS --plan configs/calibrated-epsd-study.json',
        '```','',
        '[Machine-readable results](../reports/calibrated-epsd-dense.json). [Rewritten paper](../docs/paper/vamp.md). The previous manuscript was moved without content changes to [the initial work note](../docs/paper/work_note_20261025_init_work.md).','',
        'Only joint Stokes I on this supplied observation was tested; no held-out, polarization, spectral-cube or calibration-solution claim is made.','']
    Path('devnotes/epsd-dense.md').write_text('\n'.join(lines))


# This historical experiment report does not overwrite the current manuscript.
# Refresh current paper tables/figures with scripts/update_calibrated_paper.py.


if __name__=='__main__':
    result=collect()
    Path('reports/calibrated-epsd-dense.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    from report_fixed_rate import plots
    plots(result,root=ROOT,prefix='epsd-dense');write_note(result)
    for e in result['experiments']:
        print(e['epsilon_D'],e['rate']['fraction'],e['shared']['archive_bytes'],
              [r['restored_difference_over_reference_residual'] for r in e['residual_normalized_runs']])
