"""Report the ten-percent payload experiment with the six-width ladder."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from varquanms.imaging import _grid, file_sha256, read_fits
from report_fixed_rate import fraction, payload_distribution, plots

ROOT=Path('results/calibrated-fixed10')
REFERENCE=Path('results/calibrated-template/reference')
MASK=Path('results/calibrated-template/above10deg.fits')
ORIGINAL_BYTES=1431429120
BUDGET_BYTES=143142912


def moments(values):
    return dict(pixels=int(values.size),rms=float(np.sqrt(np.mean(values*values))),
                mean=float(values.mean()),standard_deviation=float(values.std()))


def profile_label(entry):
    return ('Gamma search' if entry['config']['power_gamma_policy']=='search' else 'Max/std')+' + '+entry['family']


def collect(partial=False):
    plan=json.loads((ROOT/'screen/plan.json').read_text())
    interim=json.loads((ROOT/'interim-plan.json').read_text())
    assert plan[:3]==interim
    reference,header=read_fits(REFERENCE/'image-image.fits')
    residual,residual_header=read_fits(REFERENCE/'image-residual.fits')
    assert _grid(header)==_grid(residual_header)
    mask=fits.getdata(MASK).astype(bool)
    assert mask.shape==reference.shape and int(mask.sum())==2245726
    assert np.all(np.isfinite(reference[mask])) and np.all(np.isfinite(residual[mask]))
    reference_statistics=moments(residual[mask]);entries=[]
    for selected in plan:
        name=selected['name'];read=lambda label:json.loads((ROOT/(name+'-'+label+'.json')).read_text())
        if partial and not (ROOT/(name+'-comparison.json')).exists():continue
        rate=read('rate');shared=read('compact');export=read('export')
        manifest=json.loads((ROOT/(name+'.bundle')/'manifest.json').read_text())
        assert manifest['config']==selected['config']
        assert selected['config']['bits']==[3,4,6,8,12,16]
        assert export['status']=='passed' and export['source_files_unchanged']
        assert shared['all_partition_inventories_verified']
        assert rate['original_complex64_bytes']==ORIGINAL_BYTES and rate['within_target']
        assert rate['physical_payload_within_target']
        assert shared['stored_bytes_by_category']['visibility_block_payload']==rate['physical_payload_bytes']<=BUDGET_BYTES
        assert sum(shared['stored_bytes_by_category'].values())==shared['archive_bytes']
        bits=Counter();samples=Counter();gammas=Counter();serialized=0
        for part in manifest['partitions']:
            assert all(q['original_bits_exact'] for q in part['summary']['weights_quality'].values())
            for product in part['summary']['visibility_block_statistics']['CORRECTED_DATA']['products']:
                bits.update(product['component_bit_width_histogram'])
                samples.update(product['complex_sample_bit_width_histogram'])
                gammas.update(product['power_gamma_block_histogram'])
                serialized+=product['serialized_block_bytes']
        assert serialized==rate['payload_bytes'] and 8*sum(samples.values())==ORIGINAL_BYTES
        entry=dict(**selected,rate=rate,shared=shared,encoding=read('encode'),export_verification=export,
            visibility=read('visibility'),sample_bits=dict(samples),block_bits=dict(bits),gamma_blocks=dict(gammas),
            archive_sha256=file_sha256(Path(shared['archive'])),weights_all_exact=True)
        entry['payload_distribution']=payload_distribution(shared['archive'],by_class=True)
        assert entry['payload_distribution']['samples_by_bits']==dict(samples)
        assert sum(entry['payload_distribution']['serialized_bytes_by_bits'].values())==serialized
        measured=[]
        for repeat in (False,True):
            suffix='repeat-comparison' if repeat else 'comparison'
            if partial and not (ROOT/(name+'-'+suffix+'.json')).exists():continue
            comparison=read(suffix)
            assert comparison['comparison_mask']['sha256']==file_sha256(MASK)
            assert comparison['comparison_mask']['pixels']==int(mask.sum())
            for product in comparison['products'].values():
                assert file_sha256(product['reference'])==product['reference_file_sha256']
                assert file_sha256(product['candidate'])==product['candidate_file_sha256']
            directory=ROOT/(name+('-repeat-images' if repeat else '-images'))
            complete=json.loads((directory/'complete.json').read_text())
            assert complete['status']=='complete' and complete['source_unchanged']
            image,ih=read_fits(directory/'image-image.fits')
            image_residual,rh=read_fits(directory/'image-residual.fits')
            assert _grid(ih)==_grid(header)==_grid(rh)
            assert np.all(np.isfinite(image[mask])) and np.all(np.isfinite(image_residual[mask]))
            delta=moments((image-reference)[mask]);candidate_residual=moments(image_residual[mask])
            old=comparison['products']['image-image.fits']['regions']['above_elevation_cut']
            assert delta['rms']==old['difference_rms']
            stats=dict(run='repeat' if repeat else 'first',restored_difference=delta,
                restored_difference_over_reference_residual=delta['rms']/reference_statistics['rms'],
                candidate_residual=candidate_residual,
                candidate_residual_rms_relative_change=candidate_residual['rms']/reference_statistics['rms']-1,
                dirty_difference_over_reference_residual=comparison['products']['image-dirty.fits']['regions']['above_elevation_cut']['difference_rms']/reference_statistics['rms'])
            measured.append(stats)
            if repeat:entry['repeat']=dict(versus_reference=comparison,imaging=complete)
            else:entry.update(comparison=comparison,imaging=complete)
        entry['residual_normalized_runs']=measured
        if len(measured)==2:
            entry['worst_residual_normalized_difference']=max(x['restored_difference_over_reference_residual'] for x in measured)
        entries.append(entry)
    previous=json.loads(Path('reports/calibrated-fixed25.json').read_text())
    previous=next(e for e in previous['experiments'] if e['name']==previous['empirical_best_worst_of_two'])
    result=dict(schema='varquanms.fixed10-study/1',status='partial' if partial else 'complete',
        original_complex64_bytes=ORIGINAL_BYTES,payload_budget_bytes=BUDGET_BYTES,
        screen=json.loads((ROOT/'screen/screen.json').read_text()),experiments=entries,
        reference_residual_statistics=reference_statistics,
        reference_residual=dict(path=str(REFERENCE/'image-residual.fits'),sha256=file_sha256(REFERENCE/'image-residual.fits')),
        comparison_mask=dict(path=str(MASK),sha256=file_sha256(MASK),pixels=int(mask.sum()),minimum_elevation_deg=10),
        previous_25_percent_control=dict(name=previous['name'],bits=previous['config']['bits'],
            payload_bytes=previous['rate']['payload_bytes'],archive_bytes=previous['shared']['archive_bytes'],
            restored_difference_rms=[previous['comparison']['products']['image-image.fits']['regions']['above_elevation_cut']['difference_rms'],
                previous['repeat']['versus_reference']['products']['image-image.fits']['regions']['above_elevation_cut']['difference_rms']]),
        tests_passed=103,hardware=json.loads(Path('reports/calibrated-error-limits.json').read_text())['hardware'],
        selection_objective='Minimize the larger restored-difference RMS / reference residual RMS across two CLEAN runs among actual <=10% payload finalists.',
        limitations=['All 15 integrations, 48 channels, four correlations archived; images are joint Stokes I above 10 degrees only.',
            'Calibrated-only: original DATA/MODEL_DATA omitted, flagged visibilities zeroed and flags retained. Weights are checked exact here.',
            'Residual RMS includes sky and deconvolution structure, not a validated thermal-noise estimate. CLEAN runs can differ on identical inputs.',
            'The earlier 25% control used 4/8/12/16 widths; this is not a one-variable experiment isolating rate from the changed ladder.',
            'Pilot thresholds were chosen using bytes, before imaging. Layout variants sample different original blocks and are only a pilot.',
            'Shared host with overlapping stages and no cache eviction. Timings are observations, not isolated speed rankings.',
            'No held-out observation, Q/U/V, calibration-solution or spectral-cube validation.'])
    if not partial:
        assert len(entries)==len(plan) and all('repeat' in e for e in entries)
        result['empirical_best_worst_of_two']=min(entries,key=lambda e:e['worst_residual_normalized_difference'])['name']
        result['best_first_dirty']=min(entries,key=lambda e:e['residual_normalized_runs'][0]['dirty_difference_over_reference_residual'])['name']
        imported=json.loads(Path('reports/calibrated-import.json').read_text())
        assert all(file_sha256(Path(imported['extracted_ms'])/n)==info['sha256'] for n,info in imported['files'].items())
        result['source_integrity']=dict(original_import_checksums_match=True,files_checked=len(imported['files']))
    return result


def note(result):
    best=next(e for e in result['experiments'] if e['name']==result['empirical_best_worst_of_two'])
    old=result['previous_25_percent_control']
    old_error=max(old['restored_difference_rms'])/result['reference_residual_statistics']['rms']
    Path('configs/calibrated-fixed10-best.json').write_text(json.dumps(best['config'],indent=2)+'\n')
    lines=['# Ten-percent visibility payload with 3/4/6/8/12/16-bit selection','',
        f"The best worst-of-two restored-image comparison among the tested finalists is **{best['name']}**. Its physical visibility payload is **{best['rate']['physical_payload_bytes']/1e6:.3f} MB ({100*best['rate']['fraction']:.3f}%)** and its complete archive is **{best['shared']['archive_bytes']/1e6:.3f} MB**. [Exact opt-in configuration](../configs/calibrated-fixed10-best.json).",'',
        f"The smaller file has a measurable fidelity cost: its larger restored-image difference is **{100*best['worst_residual_normalized_difference']:.3f}% of reference residual RMS**, versus {100*old_error:.3f}% for the earlier 25%-budget selection. That is {best['worst_residual_normalized_difference']/old_error:.2f} times the difference RMS. The 10% target is achievable, but these measurements do not support calling its image agreement equivalent to the 25% result.",'',
        'The ceiling is **143,142,912 bytes**, 10% of the original 1,431,429,120 complex64 visibility bytes. It includes serialized factors, block headers, exceptions and exact fallback after secondary compression. Permutations, indexes, coordinates, flags, weights, metadata and quality arrays are counted separately.','',
        'All candidates use the [new 3/4/6/8/12/16 ladder](precision-ladder-3bit.md), calibrated-only CORRECTED_DATA, component medians, independent amplitude order per product, and unchanged time/channel resolution. Original DATA/MODEL_DATA are omitted and flagged values zeroed while flags remain. Four-bit weights are exact on this dataset.','',
        'The pilot screened three layouts/codec choices × three error rules with 14 byte-count refinements per rule, targeting 9.9%. Absolute, block-relative and hybrid 128×48 max/std profiles were retained, plus the alternative with the lowest sampled cross-parallel-hand MSE. Grouping pilots sample different rank blocks, so their visibility errors are nomination criteria, not imaging results. All full files must fit 10% before imaging. Actual rates differ slightly; this compares candidates under a common ceiling, not at identical bytes.','',
        '| Profile | Group | Gamma | Absolute component MSE cap | Relative component RMS fraction |',
        '| --- | --- | --- | ---: | ---: |']
    for e in result['experiments']:
        c=e['config'];relative=(c['relative_bounds'] or {}).get('max_component_rmse')
        absolute=c['max_component_mse']
        lines.append(f"| {profile_label(e)} | {c['baseline_group']}×{c['channel_segment']}×1 | {c['power_gamma_policy']} | {absolute if absolute is not None else '—'} | {relative if relative is not None else '—'} |")
    lines+=['','The [error-rule definitions](absolute_relative_error.md) are unchanged: component MSE ≤B, or ≤r²P, or their intersection; maximum complex error ≤10√budget, complex bias ≤√budget, and baseline component MSE ≤2×budget. P is block signal power, not thermal noise. Max/std uses R0=4; search tries gamma 1, 1/2 and 1/4.','',
        '## Images relative to the reference residual','',
        f"Reference residual RMS is **{result['reference_residual_statistics']['rms']:.9f} Jy/beam**, over the same 2,245,726 above-10° pixels. The numerator is RMS(quantized restored image − reference restored image), without mean removal or clipping. This denominator includes residual sky structure and is not a thermal-noise estimate.",'',
        '| Profile | Payload MB (%) | Whole MB | Difference RMS mJy/beam: run 1 / 2 | Difference / residual RMS %: run 1 / 2 | Dirty difference / residual RMS % |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for e in result['experiments']:
        a,b=e['residual_normalized_runs']
        lines.append(f"| {profile_label(e)} | {e['rate']['payload_bytes']/1e6:.3f} ({100*e['rate']['fraction']:.3f}%) | {e['shared']['archive_bytes']/1e6:.3f} | {1000*a['restored_difference']['rms']:.4f} / {1000*b['restored_difference']['rms']:.4f} | {100*a['restored_difference_over_reference_residual']:.4f} / {100*b['restored_difference_over_reference_residual']:.4f} | {100*a['dirty_difference_over_reference_residual']:.4f} |")
    lines+=['','| Profile | Candidate residual RMS Jy/beam: run 1 / 2 | Worst restored difference / restored-reference RMS % | Worst faint-region difference % | Largest absolute pixel difference Jy/beam | Largest reference-peak change % |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for e in result['experiments']:
        comparisons=(e['comparison'],e['repeat']['versus_reference'])
        products=[c['products']['image-image.fits'] for c in comparisons]
        a,b=e['residual_normalized_runs']
        lines.append(f"| {profile_label(e)} | {a['candidate_residual']['rms']:.6f} / {b['candidate_residual']['rms']:.6f} | {100*max(fraction(c) for c in comparisons):.4f} | {100*max(fraction(c,'image','faint_reference') for c in comparisons):.4f} | {max(p['regions']['above_elevation_cut']['difference_max_abs'] for p in products):.6f} | {100*max(abs(p['relative_change_at_reference_peak']) for p in products):.5f} |")
    lines+=['','Faint pixels are defined from the fixed restored reference below 1% of its masked absolute peak; they are not a blank-sky noise estimate. The peak comparison is at the brightest absolute reference pixel, not a solar identification or integrated-flux measurement. Similar candidate residual RMS does not establish image agreement: a different CLEAN model can redistribute sky structure while leaving this scalar nearly unchanged.','']
    lines+=['The hybrid and relative profiles are close in global restored RMS (2.510% versus 2.526% of the restored-reference RMS). Their small separation is not evidence of a universal ranking. Full gamma search gives the lowest dirty-image difference and the smallest peak change (0.092%, versus 0.517% for hybrid), despite a larger global restored difference. The saved hybrid preset follows the stated global-RMS objective; a peak-focused objective would favor a different tested profile.','',
        'At this 10% ceiling, relative and hybrid limits improve global restored agreement over the direct max/std absolute rule. The earlier 25% selection favored absolute limits. Neither absolute nor relative visibility thresholds alone determine final CLEAN quality; rate, gamma choice, weighting and deconvolution also affect the measured result.','']
    lines+=['',f"For context, the previous 25%-budget candidate used {old['payload_bytes']/1e6:.3f} MB and gave 5.238 / 4.967 mJy/beam differences (0.508 / 0.482% of this residual RMS). It used the earlier 4/8/12/16 ladder; the comparison changes both rate and available widths. Source-only CLEAN repeats previously differed by about 0.229–0.236% of restored-image RMS. These are limited repeat observations, not confidence intervals or corrections to subtract.",'',
        'All images use the requested 2048², 4 arcmin, eight-thread multiscale Stokes I template and reference restoring beam. Best means lowest larger error in these two runs on this observation, not a universal optimum or approval of a scientific tolerance.','']
    for label,file in [('Payload and residual-normalized image error','rate'),('Selected full image comparison','best-images'),
                       ('Selected reference-peak zoom','best-zoom'),('Bit distribution by samples and bytes','bits'),
                       ('All first-run differences','images'),('All reference-peak zooms','zoom')]:
        lines += [f'![{label}](figs/fixed10-{file}.png)','']
    lines+=['## Accounting and timing','',
        '| Profile | Payload + permutation/index MB | Encode only s | Verified bundle s | Compaction s | Export s | First CLEAN s |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for e in result['experiments']:
        lines.append(f"| {profile_label(e)} | {e['rate']['payload_plus_permutations_and_indexes_bytes']/1e6:.3f} | {e['encoding']['sum_encode_seconds_excluding_inspection']:.2f} | {e['encoding']['elapsed_seconds_including_copy_inspection_verification']:.2f} | {e['shared']['elapsed_seconds']:.2f} | {e['export_verification']['export']['elapsed_seconds']:.2f} | {e['imaging']['wsclean_elapsed_seconds']:.2f} |")
    categories=best['shared']['stored_bytes_by_category'];distribution=best['payload_distribution']
    all_three=sum(e['payload_distribution']['samples_by_bits'].get('3',0) for e in result['experiments'])
    cross_three=sum(e['payload_distribution']['samples_by_baseline_class_and_bits']['cross'].get('3',0) for e in result['experiments'])
    if all_three and not cross_three:
        lines+=['','Every selected three-bit sample is an autocorrelation; none is a cross-correlation. The cross-correlation streams use the other widths. This study therefore does not establish three-bit cross-correlation image fidelity. The JSON records the sample distribution separately for ordered auto and cross baseline classes.','']
    lines+=['',f"For the selected file, error-accounting arrays occupy {(categories['quantization_statistics']+categories['visibility_quality'])/1e6:.3f} MB, coordinates {categories['coordinates']/1e6:.3f} MB, and permutations/indexes {(categories['permutations']+categories['visibility_indexes'])/1e6:.3f} MB. These costs are outside the visibility-payload ceiling and inside the complete-file size.",'',
        f"Exact positive-zero blocks contain {100*distribution['exact_positive_zero_samples']/(ORIGINAL_BYTES/8):.3f}% of all samples but use only {distribution['exact_positive_zero_bytes']/1e6:.3f} MB. The large float32 sample fraction in the plot therefore does not imply a large float32 byte fraction. This includes the calibrated-only zeroing policy; the ratio is not a lossless preservation claim for flagged values.",'',
        'These durations cover all 15 integrations, with overlapping work on a shared host and no cache eviction. They are not isolated speed rankings. Full conversion is still above the 15-second goal. The JSON includes every byte category, sample/byte precision distributions, flags/weights checks, timings, input hashes and the full pilot.','',
        '## Reproduction','', '```bash','export PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1',
        'python scripts/screen_fixed_rate.py INPUT.ms NEW_SCREEN_DIR \\',
        '  --target-fraction .10 --pilot-fraction .099 --bits 3 4 6 8 12 16 \\',
        '  --layouts amplitude128x48 amplitude32x24 search128x48 --max-finalists 4 --prefix tenth',
        'python scripts/run_error_limits.py --source INPUT.ms --plan NEW_SCREEN_DIR/plan.json --root NEW_RESULTS',
        'python scripts/run_error_limits.py --source INPUT.ms --plan NEW_SCREEN_DIR/plan.json --root NEW_RESULTS --repeat-only',
        '```','',
        '[Machine-readable results](../reports/calibrated-fixed10.json). All 108 original source-file hashes still match. Export verification preserves original UVW and declared metadata exactly. The codec suite has 103 passing tests. No Q/U/V, calibration-solution, held-out-observation or spectral-cube fidelity is claimed. The paper is unchanged.','']
    Path('devnotes/fixed10.md').write_text('\n'.join(lines))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--partial',action='store_true');args=p.parse_args()
    result=collect(args.partial)
    output=ROOT/'preview.json' if args.partial else Path('reports/calibrated-fixed10.json')
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    if not args.partial:plots(result,root=ROOT,prefix='fixed10');note(result)
    for e in result['experiments']:
        print(e['name'],e['rate']['fraction'],[100*r['restored_difference_over_reference_residual'] for r in e['residual_normalized_runs']])


if __name__=='__main__':main()
