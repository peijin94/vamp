"""Collect the requested CLEAN-template benchmark without expanding the paper."""
from collections import Counter
import json
import os
from pathlib import Path
import platform

import numpy as np

from varquanms.imaging import compare_run, read_fits, file_sha256
from varquanms.partitioned import bundle_summary
from varquanms.container import Reader

ROOT = Path('results/calibrated-template')
OLD = Path('results/calibrated-20261023')
CASES = ['lossless', 'uniform', 'power', 'amplitude256', 'buckets256', 'tight001', 'stokesI-priority', 'sixbit2048']


def collect():
    result = dict(schema='varquanms.calibrated-template-tests/1',
                  issue='https://github.com/peijin94/varquanMS/issues/1',
                  mask=json.loads((ROOT/'above10deg.json').read_text()),
                  reference=json.loads((ROOT/'reference/complete.json').read_text()),
                  experiments=[])
    result['hardware'] = dict(host=platform.node(), cpu_model=next((line.split(':',1)[1].strip()
        for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')), 'unknown'),
        logical_cpus=os.cpu_count(), physical_ram_bytes=os.sysconf('SC_PHYS_PAGES')*os.sysconf('SC_PAGE_SIZE'),
        cache_methodology='No cache eviction; shared host with concurrent encoding/imaging and unrelated workloads')
    for label in ['reference-repeat', *CASES]:
        old = label in ('lossless', 'uniform', 'power')
        images = ROOT/(label if old or label == 'reference-repeat' else label+'-images')
        comparison = compare_run(ROOT/'reference', images, ROOT/'above10deg.fits')
        (ROOT/(label+'-comparison.json')).write_text(json.dumps(comparison, indent=2)+'\n')
        invocation = json.loads((images/'complete.json').read_text())
        entry = dict(name=label, comparison=comparison, imaging=invocation)
        if label != 'reference-repeat':
            bundle = OLD/(label+'-calib.bundle') if old else ROOT/(label+'.bundle')
            manifest = json.loads((bundle/'manifest.json').read_text())
            with Reader(bundle/manifest['partitions'][0]['archive'], verify=False) as reader:
                environment = reader.manifest['environment']
            summary = bundle_summary(bundle, manifest)
            counts, samples, accounting = Counter(), Counter(), Counter()
            moments = {}
            for part in manifest['partitions']:
                s = part['summary']
                for spec in s['visibility_block_statistics'].values():
                    for p in spec['products']:
                        counts.update(p['component_bit_width_histogram'])
                        samples.update(p['complex_sample_bit_width_histogram'])
                        accounting.update(p['logical_bytes_before_secondary_compression'])
                for spec in s['visibility_quality'].values():
                    for name,m in spec['by_product'].items():
                        q = moments.setdefault(name, dict(valid_samples=0, squared_error=0., bias_real_sum=0., bias_imag_sum=0., max_complex_error=0.))
                        n = m['valid_samples']; q['valid_samples'] += n
                        q['squared_error'] += m['component_mse']*n
                        q['bias_real_sum'] += m['bias_real']*n; q['bias_imag_sum'] += m['bias_imag']*n
                        q['max_complex_error'] = max(q['max_complex_error'],m['max_complex_error'])
            products = {name:dict(valid_samples=q['valid_samples'], component_mse=q['squared_error']/q['valid_samples'],
                                  bias_real=q['bias_real_sum']/q['valid_samples'],bias_imag=q['bias_imag_sum']/q['valid_samples'],
                                  max_complex_error=q['max_complex_error']) for name,q in moments.items()}
            verification = json.loads((Path('reports')/f'calibrated-{label}-export.json' if old else ROOT/(label+'-export.json')).read_text())
            if verification['status'] != 'passed' or not verification['source_files_unchanged']:
                raise ValueError('Failed export verification')
            entry.update(archive=summary, config=manifest['config'], environment=environment, block_bits=dict(counts), sample_bits=dict(samples),
                         visibility_accounting=dict(accounting), visibility_by_product=products,
                         weights_all_exact=all(q['original_bits_exact'] for p in manifest['partitions'] for q in p['summary']['weights_quality'].values()),
                         export_verification=dict(status='passed', source_files_unchanged=True, checked_columns=len(verification['column_checks']),
                                                  export=verification['export']),
                         peak_encode_rss_bytes=max(p['summary']['peak_process_rss_bytes'] for p in manifest['partitions']))
        result['experiments'].append(entry)
    result['screen'] = json.loads((ROOT/'screen/screen.json').read_text())
    result['first_integration_controls'] = json.loads((ROOT/'group-controls.json').read_text())
    result['first_integration_controls']['experiments'].append(json.loads((ROOT/'amplitude128-mse0.25-whole-first-integration.json').read_text()))
    result['first_integration_controls']['experiments'].append(json.loads((ROOT/'original32-mse0.25-whole-first-integration.json').read_text()))
    result['storage_recompression_probe'] = json.loads((ROOT/'repack-probe-final/report.json').read_text())
    result['stokesI_parallel_hands_control'] = json.loads((ROOT/'stokesI-parallel-hands.json').read_text())
    result['candidate_repeats'] = {}
    for label in ('amplitude256', 'buckets256'):
        result['candidate_repeats'][label] = dict(
            versus_reference=compare_run(ROOT/'reference', ROOT/(label+'-repeat-images'), ROOT/'above10deg.fits'),
            versus_first_candidate_run=compare_run(ROOT/(label+'-images'), ROOT/(label+'-repeat-images'), ROOT/'above10deg.fits'),
            imaging=json.loads((ROOT/(label+'-repeat-images')/'complete.json').read_text()))
    imported = json.loads(Path('reports/calibrated-import.json').read_text())
    source = Path(imported['extracted_ms'])
    if any(file_sha256(source/name) != metadata['sha256'] for name,metadata in imported['files'].items()):
        raise ValueError('Source differs from original import checksums')
    result['source_integrity'] = dict(matches_original_import=True, files_checked=len(imported['files']))
    result['limitations'] = [
        'Only joint Stokes I imaging of this 15-integration dataset; no Q/U/V, calibration solutions, light curves or held-out observations tested.',
        'Eight-thread lossless and repeated-source runs measure imaging variability; two controls do not establish its statistical distribution.',
        'All image metrics, peaks and apertures use the same all-integration >10 degree elevation mask.',
        'Archive ratios relative to original input include intentional raw/model/flagged-value removal; use the lossless calibrated-only comparator for incremental quantization.',
        'Shared host, warm/unspecified caches and concurrent workloads. Encode-only excludes inspection; verified-bundle includes partition copy, inspection and verification.',
        'The Stokes-I-priority profile relaxes XY/YX bounds; no polarization fidelity claim is supported.',
        'Selection uses this dataset and its images; no general archive defaults or scientific acceptance are established.'
    ]
    return result


def plot(result):
    os.environ.setdefault('MPLCONFIGDIR','/tmp/varquanms-mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from astropy.io import fits
    out = Path('devnotes/figs'); out.mkdir(exist_ok=True)
    fig, ax = plt.subplots(figsize=(8,5),constrained_layout=True)
    controls = []
    for e in result['experiments']:
        metric = e['comparison']['products']['image-image.fits']['regions']['above_elevation_cut']
        error = 100*metric['relative_rms_difference']
        if e['name'] in ('reference-repeat','lossless'): controls.append(error)
        if 'archive' in e:
            mb = e['archive']['archive_bytes']/1e6
            ax.scatter(mb,error)
            offset = {'stokesI-priority':(-35,7), 'lossless':(-38,7)}.get(e['name'], (4,5))
            ax.annotate(e['name'],(mb,error),xytext=offset,textcoords='offset points',fontsize=8)
            if e['name'] in result['candidate_repeats']:
                repeated = result['candidate_repeats'][e['name']]['versus_reference']['products']['image-image.fits']
                second = 100*repeated['regions']['above_elevation_cut']['relative_rms_difference']
                ax.plot([mb, mb], [error, second], color='black', linewidth=.8)
                ax.scatter(mb, second, marker='x', color='black', s=25)
    ax.axhspan(min(controls),max(controls),alpha=.15,color='gray',label='Two observed control differences')
    ax.set(xlabel='Whole archive (decimal MB)',ylabel='Restored image relative RMS difference (%)',
           title='2048², 4 arcmin, requested CLEAN; elevation >10° at all integrations')
    ax.set_yscale('log')
    ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.savefig(out/'calibrated-test-rate.png',dpi=170); plt.close(fig)
    mask = fits.getdata(ROOT/'above10deg.fits').astype(bool)
    reference, _ = read_fits(ROOT/'reference/image-image.fits')
    labels = CASES
    differences=[]
    for label in labels:
        path=ROOT/(label if label in ('lossless','uniform','power') else label+'-images')/'image-image.fits'
        candidate,_=read_fits(path); differences.append(candidate-reference)
    limit=max(float(np.quantile(np.abs(d[mask]),.999)) for d in differences)
    fig,axes=plt.subplots(3,3,figsize=(12,12),constrained_layout=True)
    peak=float(np.max(np.abs(reference[mask])))
    im=axes.flat[0].imshow(np.where(mask,reference,np.nan),origin='lower',cmap='RdBu_r',norm=SymLogNorm(linthresh=peak*.001,vmin=-peak,vmax=peak))
    axes.flat[0].set_title('Reference (symmetric log)');fig.colorbar(im,ax=axes.flat[0],shrink=.75)
    for ax,label,d in zip(list(axes.flat)[1:],labels,differences):
        im=ax.imshow(np.where(mask,d,np.nan),origin='lower',cmap='RdBu_r',
                     norm=SymLogNorm(linthresh=limit*.01,vmin=-limit,vmax=limit))
        e=next(e for e in result['experiments'] if e['name']==label)
        metric=e['comparison']['products']['image-image.fits']['regions']['above_elevation_cut']
        ax.set_title(f"{label}: ΔRMS {100*metric['relative_rms_difference']:.3f}%")
        fig.colorbar(im,ax=ax,shrink=.75)
    for ax in axes.flat: ax.set(xlabel='Pixel x',ylabel='Pixel y')
    fig.suptitle('Jy/beam; >10° mask; common symmetric log difference scale; extremes clipped at max(candidate p99.9)')
    fig.savefig(out/'calibrated-test-differences.png',dpi=170);plt.close(fig)
    plot_issue_panels(result, out)


def plot_issue_panels(result, out):
    """Show the source neighbourhood and both block/sample precision fractions."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    from astropy.io import fits
    mask = fits.getdata(ROOT/'above10deg.fits').astype(bool)
    reference, _ = read_fits(ROOT/'reference/image-image.fits')
    candidate, _ = read_fits(ROOT/'amplitude256-images/image-image.fits')
    selected = next(e for e in result['experiments'] if e['name'] == 'amplitude256')
    cy, cx = selected['comparison']['products']['image-image.fits']['reference_peak_yx']
    y0, y1 = max(0, cy-64), min(reference.shape[0], cy+64)
    x0, x1 = max(0, cx-64), min(reference.shape[1], cx+64)
    crop = np.s_[y0:y1, x0:x1]
    valid = mask[crop]
    a, b = reference[crop], candidate[crop]
    peak = max(float(np.max(np.abs(a[valid]))), float(np.max(np.abs(b[valid]))))
    delta = b-a
    limit = float(np.max(np.abs(delta[valid]))) or 1.
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.6), constrained_layout=True)
    for ax, values, label in zip(axes, [a, b, delta], ['Reference', 'Amplitude 128×48', 'Decoded − reference']):
        difference = ax is axes[2]
        norm = None if difference else SymLogNorm(linthresh=peak*.001, vmin=-peak, vmax=peak)
        options = dict(vmin=-limit, vmax=limit) if difference else dict(norm=norm)
        im = ax.imshow(np.where(valid, values, np.nan), origin='lower', cmap='RdBu_r',
                       extent=(x0-.5, x1-.5, y0-.5, y1-.5), **options)
        ax.set(title=label, xlabel='Pixel x', ylabel='Pixel y')
        fig.colorbar(im, ax=ax, shrink=.8, label='Jy/beam')
    fig.suptitle('128-pixel zoom at masked reference peak; 4 arcmin/pixel; elevation >10°\n'
                 'Reference/decoded share a symmetric-log scale; difference uses its full linear range')
    fig.savefig(out/'calibrated-test-zoom.png', dpi=180); plt.close(fig)

    entries = [e for e in result['experiments'] if e.get('block_bits')]
    widths = ['4', '6', '8', '12', '16', 'float32']
    colors = ['#9467bd', '#17becf', '#1f77b4', '#2ca02c', '#ff7f0e', '#7f7f7f']
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), constrained_layout=True, sharey=True)
    y = np.arange(len(entries))
    for ax, field, title in zip(axes, ['block_bits', 'sample_bits'], ['Fraction of blocks', 'Fraction of complex samples']):
        totals = np.array([sum(e[field].values()) for e in entries])
        left = np.zeros(len(entries))
        for width, color in zip(widths, colors):
            fractions = np.array([e[field].get(width, 0) for e in entries])/totals*100
            ax.barh(y, fractions, left=left, color=color, label=width+' bit' if width != 'float32' else 'float32 lossless')
            left += fractions
        ax.set(title=title, xlabel='Percent', xlim=(0, 100), yticks=y,
               yticklabels=[e['name'] for e in entries])
        ax.grid(axis='x', alpha=.2)
    axes[0].invert_yaxis()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=6, fontsize=9)
    fig.suptitle('Selected precision per real/imaginary component, all 15 integrations\n'
                 'Fractions are not byte shares; lossless includes cheaply compressed all-flagged zero blocks')
    fig.savefig(out/'calibrated-test-bit-distribution.png', dpi=180); plt.close(fig)


def write_note(result):
    def product(e, name='image'):
        return e['comparison']['products'][f'image-{name}.fits']
    def percent(p, region='above_elevation_cut'):
        return 100*p['regions'][region]['relative_rms_difference']
    entries = {e['name']:e for e in result['experiments']}
    lines = ['# Calibrated data compression tests', '',
        'Completed 2026-10-23/24. Tracking: [Calibrated data compression #1](https://github.com/peijin94/varquanMS/issues/1). The paper is unchanged.', '',
        'Input: `lwa-15slot-flat.ms`, 15 integrations, 931,920 rows, 48 channels (63.966–68.464 MHz), XX/XY/YX/YY. '
        'The source is 4,428,203,949 bytes. All 108 files still match the original import checksums. '
        'Every archive uses the same calibrated-only policy: omit original DATA/MODEL_DATA, replace flagged values with zero, '
        'preserve flags and original UVW/metadata, and encode weights with four bits. Weights are **exact** on this dataset. '
        'Default exports recreate calibrated DATA from decoded CORRECTED_DATA.', '',
        '## CLEAN template and comparison region', '', '```text',
        'wsclean -size 2048 2048 -scale 4amin -data-column CORRECTED_DATA',
        '  -j 8 -pol I -multiscale -multiscale-scale-bias 0.8',
        '  -niter 150000 -mgain 0.95 -auto-threshold 0.5 -auto-mask 3',
        '  -local-rms -horizon-mask 10deg -weight briggs 0 -mem 40',
        '  -reorder -no-update-model-required', '```', '',
        'All runs use the same pinned WSClean 3.7 container, private MS copies, and the reference restoring beam. '
        '`-mem 40` is a percentage of system memory. Exact commands, container IDs, source/FITS checksums, '
        'beam, iterations and software versions are recorded in the [machine-readable report](../reports/calibrated-template-tests.json).', '',
        'All image metrics, peaks, apertures and plots use **2,245,726 pixels above 10° elevation at every one of the 15 original TIME epochs**. '
        'The mask uses celestial FITS WCS, mean antenna ITRF position, UTC times, Astropy and zero refraction. '
        'Original UVW is never regenerated. The independently computed mask differs slightly from WSClean’s own horizon mask '
        '(2,250,505 pixels); the same independent mask is used for every comparison. IERS final values cover all epochs. '
        'The retained table from `https://maia.usno.navy.mil/ser7/finals2000A.all` has SHA256 '
        '`95a61292c957fa74ed624bd4f7a84b1b984077e412b43ed39b99e092f38de9da`.', '',
        r'Reported relative RMS is $100\sqrt{\langle(\widehat I-I)^2\rangle}/\sqrt{\langle I^2\rangle}$ over that mask. '
        r'It is not a thermal-noise measurement. Faint pixels have $|I|<0.01\max|I|$ in the masked reference. '
        'Apertures contain pixels within radius four of the brightest absolute reference pixel; they are not identified solar photometry.', '',
        '## Whole-dataset image results', '',
        '| Profile | Whole archive (MB) | Dirty difference (%) | Restored difference (%) | Faint-region difference (%) | Peak change (%) |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for e in result['experiments']:
        p = product(e)
        size = f"{e['archive']['archive_bytes']/1e6:.3f}" if 'archive' in e else '—'
        lines.append(f"| {e['name']} | {size} | {percent(product(e,'dirty')):.6f} | {percent(p):.4f} | {percent(p,'faint_reference'):.4f} | {100*p['relative_change_at_reference_peak']:+.6f} |")
    lines += ['', 'MB is decimal. `lossless` means lossless visibilities under the same destructive calibrated-only policy. '
        'Ratios against the original 4.428 GB include intentionally omitted columns and flagged values; '
        'use the 1.080 GB calibrated-only lossless archive to isolate the incremental benefit of quantization.', '',
        '| Repeated candidate | First versus reference (%) | Repeat versus reference (%) | Repeat versus first (%) |',
        '| --- | ---: | ---: | ---: |']
    for label, repeat in result['candidate_repeats'].items():
        a = repeat['versus_reference']['products']['image-image.fits']
        b = repeat['versus_first_candidate_run']['products']['image-image.fits']
        lines.append(f"| {label} | {percent(product(entries[label])):.4f} | {percent(a):.4f} | {percent(b):.4f} |")
    lines += ['', '**Working choice for this dataset and Stokes I:** `amplitude256`, 128 baselines × 48 channels × one product. '
        'It gives 509.7 MB and 0.55–0.69% restored-image differences in two runs, versus 1,079.6 MB for the same-policy lossless archive. '
        'Its first-run maximum image difference is 0.642 Jy/beam and faint-region relative RMS is 1.045%; small global RMS does not mean every weak feature is unchanged.', '',
        '`sixbit2048` is smaller at 447.1 MB but gives 1.270% restored and 1.933% faint-region differences. '
        'This explicitly adds a six-bit experimental rung; the default precision ladder remains 4/8/12/16. '
        '**Reject the tested `buckets256` configuration for this image-fidelity objective:** it saves only about 1 MB relative to amplitude256, '
        'but repeatedly produces about 9.2% restored and 14.1–14.2% faint-region differences. '
        'Visibility component MSEs are nearly equal (172.508 amplitude versus 174.335 buckets), and dirty-image differences are close. '
        'The large difference appears after nonlinear CLEAN; its detailed cause is not established. These configurations also differ in block size and channel segmentation, so this is not an isolated test of ordering.', '',
        '`stokesI-priority` saves about 62 MB relative to the earlier uniform archive by relaxing XY/YX bounds. '
        'Across all 15 partitions, 150 checksum comparisons confirm identical XX/YY payloads and reconstruction mappings to `uniform`. '
        'Their different CLEAN results are imaging variability, not an improvement in Stokes I quantization. '
        'Neither this profile nor these Stokes I tests establish Q/U/V fidelity. Two reference controls and two candidate repeats do not establish a full statistical distribution of CLEAN variability.', '',
        '![Size versus restored-image difference](figs/calibrated-test-rate.png)', '',
        'Vertical axis is logarithmic. Black crosses/segments show the two observed candidate runs, not confidence intervals. '
        'The shaded band spans the two lossless/reference controls.', '',
        '![Reference, amplitude256 reconstruction and difference above 10 degrees](figs/calibrated-test-amplitude.png)', '',
        '![Zoom around the brightest masked reference pixel](figs/calibrated-test-zoom.png)', '',
        'The 128-pixel zoom spans 8.53 degrees in image coordinates. Reference and decoded panels share their colour scale; '
        'the difference has a separate linear scale. All panels exclude pixels at or below the elevation cut.', '',
        '![Selected component bit widths by block and complex sample](figs/calibrated-test-bit-distribution.png)', '',
        'Precision fractions count blocks or complex samples, not stored bytes. Lossless includes all-flagged zero blocks; '
        'their compressed cost can be small. No four-bit blocks were selected in amplitude256.', '',
        '[All tested profiles, including the failed bucket configuration](figs/calibrated-test-differences.png). '
        'These diagnostic difference maps share a symmetric-log scale and clip extremes at the largest candidate 99.9th-percentile absolute difference; all numerical metrics are unclipped.', '',
        '## Grouping and precision controls', '',
        'A deterministic 28-case screen sampled up to 128 groups per product in the first integration. It covered per-baseline scaling, '
        'original-order groups, amplitude sorting, factor-of-two amplitude buckets, length/orientation grouping and median/power quantization, '
        'with component-MSE ceilings 0.01/1/64/256. Different groupings sampled different values, and these payload rates exclude metadata; '
        'the screen is a shortlist, not a whole-MS rate-distortion proof. All sampled results are in the JSON report.', '',
        'The following controls encode and verify **every sample of the first integration**, counting whole archives:', '',
        '| First-integration control | Whole archive (MB) | Achieved component MSE | Encode call including inspection (s) |',
        '| --- | ---: | ---: | ---: |']
    for e in result['first_integration_controls']['experiments']:
        s=e['summary']; q=s['visibility_quality']['CORRECTED_DATA']['by_product']; n=sum(v['valid_samples'] for v in q.values())
        mse=sum(v['component_mse']*v['valid_samples'] for v in q.values())/n
        lines.append(f"| {e['name']} | {s['archive_bytes']/1e6:.3f} | {mse:.6f} | {s['elapsed_seconds']:.2f} |")
    lines += ['', '`baseline` uses 1×24 blocks without sorting; the `32` controls use 32×24; `amplitude128` uses 128×48. '
        '`original32` keeps original order within auto/cross and no-valid-data classes, without amplitude sorting. '
        'Grouping can improve both size and error in this prototype: original32 at MSE ceiling 0.25 uses 50.34 MB with achieved MSE 0.01046, '
        'versus 56.15 MB and MSE 0.40026 for per-baseline scaling at ceiling 1. The latter also spends much longer encoding tiny blocks. '
        'For matched 32×24 shapes and ceiling 1, amplitude sorting is 1.8% smaller than buckets but has slightly higher achieved MSE; '
        'that does not establish a sorting advantage at equal error.', '',
        'All quantized full-dataset candidates use independent per-product orderings, 4/8/12/16-bit component rails (plus 6 in the explicit experiment), '
        'max/.995 range candidates, float32 scales, exceptions and exact float32 fallback. The uniform candidates have no offset; '
        'the power control uses component medians and gamma choices 1/.5. Bounds use stored CORRECTED_DATA units, with no thermal-noise claim. '
        'For MSE ceiling B, the tested baseline-MSE, maximum complex error and complex bias ceilings are 2B, 10√B and √B. '
        'Names indicate B: uniform/power=1, amplitude256/buckets256=256, tight001=.01, sixbit2048=2048. '
        'Stokes-I priority uses B=1 for XX/YY and B=256 for XY/YX. Exact configurations are in `configs/calibrated-test-*.json`.', '',
        '## Time and byte accounting', '',
        '| Profile | Encoding only, summed over 15 partitions (s) | Complete verified bundle (s) | MS export (s) | WSClean only (s) |',
        '| --- | ---: | ---: | ---: | ---: |']
    for label in CASES:
        e=entries[label]; s=e['archive']; export=e['export_verification']['export']['elapsed_seconds']
        et=f'{export:.2f}' if export is not None else 'not retained'
        lines.append(f"| {label} | {s['sum_encode_seconds_excluding_inspection']:.2f} | {s['elapsed_seconds_including_copy_inspection_verification']:.2f} | {et} | {e['imaging']['wsclean_elapsed_seconds']:.2f} |")
    hw=result['hardware']
    lines += ['', f"Host: {hw['host']}, {hw['cpu_model']}, {hw['logical_cpus']} logical CPUs. Quantized full-dataset encodes use 24 workers; first-integration controls use eight; WSClean uses eight threads. "
        'Runs shared the host and overlapped other work; caches were not evicted. These are workload observations, not isolated latency rankings. '
        'The first three archives reuse encoding measurements from the earlier session, but all images use the new template. '
        'Complete bundle time includes partition copying, inspection and verification; MS export and imaging are separate. '
        'No 15-second complete conversion of this 15-integration MS is demonstrated.', '',
        'For amplitude256, 314.57 MB is allocated visibility block payload, 83.01 MB repeated subtable arrays, '
        '41.51 MB error-accounting arrays, 21.08 MB amplitude statistics, 19.58 MB coordinates, and the rest maps/container metadata. '
        'All categories sum to the whole archive. Its 29,224 blocks select 20,825×8-bit, 1,401×12-bit, 120×16-bit and 6,878 lossless; '
        'many lossless blocks are cheaply compressed all-flagged zero blocks, not failed error budgets. No four-bit blocks were selected. '
        'Complete per-profile histograms, error statistics, exception costs and bytes are in the JSON.', '',
        'A lossless recompression probe on a 42.248 MB first-partition archive produced 41.194 MB (gzip 1) and 40.791 MB (gzip 9), '
        'with every logical-array checksum unchanged. This 2.5–3.4% gain is modest. Repeated subtables and quality arrays are significant remaining costs. '
        'Repeated subtables are candidates for exact sharing; the error arrays generally differ and cannot be assumed redundant.', '',
        '## Reproduction and limits', '', '```bash',
        'export PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1',
        'python scripts/run_wsclean.py INPUT.ms --output REFERENCE --data-column CORRECTED_DATA --clean-profile calibrated --threads 8',
        'python scripts/make_elevation_mask.py INPUT.ms REFERENCE/image-image.fits above10deg.fits --iers-file finals2000A.all',
        'python scripts/timeslot_archive.py encode INPUT.ms NEW.bundle --config configs/calibrated-test-amplitude256.json --report encode.json',
        'python scripts/validate_partition_export.py NEW.bundle INPUT.ms NEW.ms --report export.json',
        'python scripts/run_wsclean.py NEW.ms --output IMAGES --data-column CORRECTED_DATA --clean-profile calibrated --threads 8 --beam-from REFERENCE/image-image.fits',
        'python scripts/compare_wsclean.py REFERENCE IMAGES --mask above10deg.fits --report comparison.json', '```', '',
        'Use new output paths. `scripts/screen_calibrated.py` reproduces the block screen; `scripts/probe_archive_recompression.py` measures lossless repacking. '
        'Local full-resolution FITS, archives, exports, commands and logs are under ignored `results/calibrated-template/`. '
        '`scripts/report_calibrated_tests.py` collects them into this note, the figures and the JSON report.', '',
        'This is joint Stokes I imaging of one 15-integration dataset. It does not validate individual time frames, transient light curves, '
        'weak polarization, spectra, calibration solutions, absolute flux or a held-out dataset. No Dysco comparison was run here. '
        'The recommendation is specific to this experiment, not a universal scientific archive default.', '',
        'References: [WSClean local RMS](https://wsclean.readthedocs.io/en/v3.6/local_rms_thresholding.html), '
        '[WSClean masks](https://wsclean.readthedocs.io/en/v3.7/masking.html), '
        '[Astropy AltAz](https://docs.astropy.org/en/stable/api/astropy.coordinates.AltAz.html).', '']
    shared_report = Path('reports/calibrated-shared-storage.json')
    if shared_report.exists():
        shared = json.loads(shared_report.read_text())
        entries = {e['name']:e for e in shared['experiments']}
        a = entries['amplitude256']
        proof = shared['amplitude256_matches_previously_imaged_export']
        phases = shared['original_encoding_phase_seconds']
        lines += ['## Follow-up: exact sharing across integrations', '',
            'A new single-file HDF5 representation shares only bitwise-identical arrays and compresses manifest JSON. '
            'Original scientific manifests, payloads, scales, codebooks, flags, weights, UVW and reconstruction maps are unchanged. '
            'This is a storage change with no second quantization pass.', '',
            '| Profile | Directory bundle (MB) | Shared HDF5 (MB) | Additional repack including checks (s) |',
            '| --- | ---: | ---: | ---: |']
        for e in shared['experiments']:
            lines.append(f"| {e['name']} | {e['storage_compaction']['source_bundle_bytes']/1e6:.3f} | {e['archive_bytes']/1e6:.3f} | {e['elapsed_seconds']:.2f} |")
        lines += ['', '![Whole-archive sizes before and after exact sharing](figs/calibrated-shared-size.png)', '',
            f"The working amplitude profile is now **{a['archive_bytes']/1e6:.1f} MB**, versus "
            f"{entries['lossless']['archive_bytes']/1e6:.1f} MB for lossless visibility under the same storage and preservation policies. "
            f"Its independent export matches the previously imaged export in all {len(proof['columns'])} supported column/metadata checks, bit for bit. "
            f"All {proof['source_files_checked']} original source-file checksums still match. "
            'The existing 0.55–0.69% Stokes I image differences therefore describe unchanged imaging inputs; no new CLEAN run was needed to assess this lossless storage transformation. '
            'CLEAN itself remains subject to the previously measured run-to-run variability.', '',
            'All three shared archives passed full inventory checks and exact scientific-manifest comparison. '
            'The amplitude profile additionally passed full MS export/source verification and equality against the prior export. '
            'The six-bit profile retains its earlier 1.27% global and 1.93% faint-region image differences.', '',
            f"For amplitude256, exact duplicate arrays remove {sum(a['storage_compaction']['duplicate_stored_bytes_by_category'].values())/1e6:.2f} MB, "
            f"including {a['storage_compaction']['duplicate_stored_bytes_by_category']['subtables']/1e6:.2f} MB of repeated subtables. "
            'Manifest compression and container repacking account for the remaining saving. '
            'Visibility payload and error-accounting arrays still occupy 314.57 MB and 41.51 MB, respectively; the latter were not assumed redundant.', '',
            f"Repacking took {a['elapsed_seconds']:.2f} s: {a['storage_compaction']['source_validation_seconds']:.2f} s input validation, "
            f"{a['storage_compaction']['copy_and_share_seconds']:.2f} s copying/sharing, and {a['output_validation_seconds']:.2f} s output validation, plus final accounting/publication. "
            'This time is additional to creating the directory bundle. It is not a 44-second MS encoder. '
            f"The earlier {sum(phases.values()):.1f} s MS-to-bundle run spent {phases['partition_copy']:.1f} s copying partitions, "
            f"{phases['inspection']:.1f} s inspecting, {phases['encoding']:.1f} s encoding, and {phases['verification_integrity_and_other']:.1f} s in verification/integrity/other work. "
            'About 79.5% was outside encoding; eliminating repeated copies and inspections is the next speed target. '
            'These phase timings are shared-host observations, not isolated CPU profiles. The 15-second full-conversion target remains unmet.', '',
            '```bash',
            'python scripts/timeslot_archive.py compact NEW.bundle SHARED.h5 --report compact.json',
            'python scripts/validate_partition_export.py SHARED.h5 INPUT.ms SHARED.ms --report shared-export.json',
            '```', '',
            'See the [machine-readable sharing report](../reports/calibrated-shared-storage.json) and '
            '[format contract](../docs/format.md#shared-partition-storage-varquanmsshared-partitions1). '
            'The source bundle is retained; the shared file can be relocated and decoded independently. '
            'Local fixtures cover special-value bits, calibrated DATA reconstruction, direct access to a later partition, '
            'and rejection of corrupt payloads, mappings, missing data, unsafe links and incomplete files. The local suite has 95 passing tests.', '']
    if Path('reports/calibrated-gamma-ratio.json').exists():
        lines += ['Follow-up: [gamma selection from maximum/std](gamma_max_std.md) compares two ratio rules with full gamma search, '
                  'including matched-block rates, complete shared archives and the same elevation-masked CLEAN test.', '']
    if Path('reports/calibrated-error-limits.json').exists():
        lines += ['Follow-up: [absolute versus relative error limits](absolute_relative_error.md) compares absolute, '
                  'block-relative and hybrid thresholds at matched sampled payload rates, with full-data archives and CLEAN images.', '']
    if Path('reports/calibrated-fixed25.json').exists():
        lines += ['Follow-up: [fixed 25% visibility-payload budget](fixed25.md) compares five complete finalists, '
                  'physical byte counts and repeated calibrated images, with an opt-in measured configuration.', '']
    if Path('reports/calibrated-fixed10.json').exists():
        lines += ['Follow-up: [10% visibility-payload budget](fixed10.md) uses the new 3/4/6/8/12/16-bit ladder '
                  'and compares four full-data archives with two CLEAN runs each, including differences normalized '
                  'to the reference residual RMS.', '']
    if Path('reports/calibrated-epsd-dense.json').exists():
        lines += ['Follow-up: [fixed epsilon D with the expanded bit ladder](epsd-dense.md) records '
                  'three prescribed tolerances, repeated calibrated images and the rewritten methods paper.','']
    if Path('reports/calibrated-fixed16.json').exists():
        lines += ['Follow-up: [16% payload with dimensionless G/S/D/R budgets](fixed16.md) replaces fixed B, '
                  'screens 333 coefficient/layout combinations and compares ten full-data finalists '
                  'with repeated calibrated images.', '']
    Path('devnotes/test_calibrated_dataset.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    result=collect()
    Path('reports/calibrated-template-tests.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    plot(result)
    write_note(result)
    for e in result['experiments']:
        p=e['comparison']['products'];m=p['image-image.fits']
        print(e['name'],e.get('archive',{}).get('archive_bytes'),
              'dirty%',100*p['image-dirty.fits']['regions']['above_elevation_cut']['relative_rms_difference'],
              'restored%',100*m['regions']['above_elevation_cut']['relative_rms_difference'],
              'peak%',100*m['relative_change_at_reference_peak'], 'aperture%',100*m['aperture']['relative_sum_change'])
