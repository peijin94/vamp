"""Refresh paper tables/figures from frozen measurements; never rerun experiments."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re

os.environ.setdefault('MPLCONFIGDIR', '/tmp/varquanms-mpl')
ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / 'docs/paper/varquantMS.md'
FIGS = PAPER.parent / 'figs'


def load(name):
    return json.loads((ROOT / 'reports' / (name + '.json')).read_text())


def replace_section(text, name, body):
    result, count = re.subn(r'<!-- ' + name + r'_START -->.*?<!-- ' + name + r'_END -->',
                           lambda _: f'<!-- {name}_START -->\n{body.strip()}\n<!-- {name}_END -->',
                           text, flags=re.S)
    if count != 1:
        raise ValueError(f'Expected one {name} section; found {count}')
    return result


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])


def runtime_table(timing):
    rows = [[p['phase'], f"{p['seconds']:.3f}", f"{p['percent']:.2f}"]
            for p in timing['phase_wall_seconds']]
    rows.append(['Total process wall time', f"{timing['total_process_wall_seconds']:.3f}", '100.00'])
    parallel = [[name, f'{value:.3f}'] for name,value in timing['parallel_partition_stage_seconds'].items()]
    return '''The following breakdown uses the recorded median run (run 2), whose full
process wall time is 72.685 s. Rows are non-overlapping elapsed intervals, or
explicit differences of those measured intervals, and sum to the complete run.

''' + table(['Phase', 'Wall seconds', '% of total'], rows) + '''

The remaining bundle interval includes discovery, overlapping partition work
and finalization. The instrumentation does not separate these into individual
additive wall-time contributions. The two overhead rows are measured remainders.
Input-bundle validation, array copying/sharing and final-archive validation are
sequential subphases of the 19.325-second compaction stage.

For context, the same run also records the following work summed across the
15 partitions. These totals overlap in wall time and must not be added
to the preceding table:

''' + table(['Parallel work', 'Summed partition seconds'], parallel) + '''

Encoding/writing is the largest recorded partition-work total. The parent
waits while partition processes run, so attributing that interval to parent CPU
or adding thread/partition durations would misrepresent the 72.685-second
latency. A finer additive split inside parallel encoding would require new
instrumentation; none is inferred from unrelated single-context profiles.
[Median-run timing audit](../../reports/calibrated-runtime-breakdown.json).
'''


def v_validation_section(study):
    rows=[]
    for case in sorted(study['experiments'],key=lambda e:e['epsilon_D']):
        n=[r['normalization'] for r in case['runs']]
        rows.append([f"{case['epsilon_D']:g}",
            ' / '.join(f"{1000*x['difference_rms']:.4f}" for x in n),
            ' / '.join(f"{100*x['difference_over_reference_residual_rms']:.4f}" for x in n)])
    return '''A follow-up on 27 October images Stokes V at all three epsilon values,
using the original verified exports and the same template, reference beam and
above-10° mask. Each case is CLEANed twice. All six runs reach the final
threshold and preserve their input MS hashes. These new runs are separate from
the preceding performance-study V pair at 0.18; no result has been relabelled or
substituted. The matching V reference residual RMS is 0.384439444 Jy/beam.

''' + table(['εD', 'V difference mJy/beam, first / repeat',
             'V difference / reference residual RMS %, first / repeat'],rows) + '''

[Full three-epsilon V measurements](../../reports/calibrated-epsd-v.json).
Figures 7–8 use the same layout as the I comparisons in Figures 5–6. The V
reference selects a source-centred crop independently of the reconstructions.

![Recovered Stokes V and RMS-normalized percentage differences](figs/epsd-v-images.png)

*Figure 7. Stokes V counterparts to Figure 5: recovered images above signed
percentage differences, with epsilon-D 0.01, 0.04 and 0.18 from left to right.
Both rows use bwr; the intensity scale is linear with fixed
limits of −10 to +10 Jy/beam, and the RMS-normalized difference scale is linear
from −30% to +30%. Extended difference-colourbar
ends indicate values outside that display range. Only the fixed above-10° mask
is shown. Top-row x axes and the last two columns' y axes are omitted.*

![Stokes V comparisons around Taurus A](figs/epsd-v-zoom.png)

*Figure 8. The same V comparisons in a fixed 64×64-pixel crop around a compact
V feature in the Taurus A (Crab Nebula) region. The reference peak is
+6.091 Jy/beam, coincident with the bright I source and within one restoring
beam of the [catalogue position](https://simbad.u-strasbg.fr/simbad/sim-id?Ident=Crab+Nebula).
Circles mark the selected reference peak. The centre stays above 50° elevation
and every crop pixel remains above 45°
throughout the observation. Selection uses only the reference; the full-field
residual RMS denominator and colour scales are retained. This is measured V
at a radio-source position; intrinsic circular polarization and instrumental
leakage have not been separated.*
'''


def sections(study, perf, lossless, audit, vstudy, timing):
    run = perf['optimized']['runs'][0]
    total = run['archive_bytes']
    payload = perf['visibility_payload_bytes']
    original = study['original_complex64_bytes']
    cases = study['experiments']
    for name, sha in audit['input_report_sha256'].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != sha:
            raise ValueError(f'Input report changed: {name}; refresh and review the storage audit')
    assert sum(run['byte_categories'].values()) == total
    assert sum(lossless['stored_bytes_by_category'].values()) == lossless['archive_bytes']
    assert sum(perf['sample_bits'].values()) == perf['complex_samples']
    selected = next(e for e in cases if e['epsilon_D'] == .18)
    assert selected['sample_bits'] == perf['sample_bits']
    assert sum(selected['payload_distribution']['serialized_bytes_by_bits'].values()) == payload
    iv_rows = []
    control_rows = []
    for pol in ('I', 'V'):
        runs = perf['imaging'][pol]
        n = [r['normalization'] for r in runs]
        control = perf['reference_repeat_controls'][pol]['image-image.fits']
        iv_rows.append([pol, f"{n[0]['reference_residual_rms']:.6f}",
                        ' / '.join(f"{1000*x['difference_rms']:.4f}" for x in n),
                        ' / '.join(f"{100*x['difference_over_reference_residual_rms']:.4f}" for x in n)])
        control_rows.append([pol, f"{1000*control['difference_rms']:.6f}",
                             f"{100*control['difference_over_reference_residual_rms']:.6f}"])
    iv = '''The subsequent test uses the epsilon-D 0.18 archive produced by the measured
optimized workflow. Its decoded visibilities equal those of the earlier 0.18
case: 147 column comparisons and the retained metadata match. An independent
frozen decoder exports the new archive, and original-source checks pass. Across
archive implementations, 5,040 scientific dataset references have matching
shapes, dtypes and hashes, with unchanged measured visibility errors. The extra
I runs below are distinct CLEAN realizations of the same decoded values; they
do not replace the earlier I measurements. In that performance study, V was
tested at 0.18; the subsequent three-epsilon V sweep is reported below.

''' + table(['Stokes', 'Reference residual RMS Jy/beam', 'Difference mJy/beam, first / repeat',
             'Difference / reference residual RMS %, first / repeat'], iv_rows) + '''

A repeat of each uncompressed reference quantifies observed CLEAN variation:

''' + table(['Stokes', 'Reference-repeat difference mJy/beam', 'Difference / reference residual RMS %'], control_rows) + '''

All four quantized I/V CLEAN runs reached the final threshold. The V discrepancy
is about 10.06% of its reference residual RMS, compared with 4.81–4.88% in I.
Although the absolute V difference is smaller, its normalization is also smaller.
Both exceed their corresponding reference-repeat differences; the V discrepancy
cannot be explained by the observed repeat alone. These controls are not
subtracted from the reported errors. Visual agreement or a small I difference
does not establish weak-polarization or leakage fidelity. The follow-up below
measures the two stricter V cases directly.
[Measured I/V and identity report](../../reports/calibrated-performance.json).

'''
    iv += v_validation_section(vstudy)
    size_rows = [
        ['Original full MS (physical)', f"{study['original_ms_physical_bytes']/1e6:.3f}", '—', 'Includes raw and model columns'],
        ['Original CORRECTED_DATA (logical)', '—', f'{original/1e6:.3f}', 'Denominator for payload fractions'],
        ['Calibrated-only, exact retained values', f"{lossless['archive_bytes']/1e6:.3f}",
         f"{lossless['stored_bytes_by_category']['visibility_lossless']/1e6:.3f}", 'Shared single-file lossless control']]
    for e in cases:
        size_rows.append([f"Quantized, εD={e['epsilon_D']:g}", f"{e['shared']['archive_bytes']/1e6:.3f}",
                          f"{sum(e['payload_distribution']['serialized_bytes_by_bits'].values())/1e6:.3f}", 'Initial three-epsilon study'])
    size_rows.append(['Quantized, εD=0.18, optimized', f'{total/1e6:.3f}', f'{payload/1e6:.3f}', 'I/V and speed benchmark, run 1'])
    labels = {
        'visibility_block_payload': 'Visibility blocks (codes, headers, factors, exceptions, exact blocks)',
        'visibility_quality': 'Measured visibility-error records',
        'quantization_statistics': 'Quantization/reference statistics',
        'coordinates': 'Original coordinates and identities',
        'permutations': 'Baseline permutations',
        'subtables': 'Shared subtables',
        'container_manifests_verification_and_links': 'Container, manifests, verification and links',
        'reconstruction': 'Reconstruction maps',
        'visibility_indexes': 'Visibility-block indexes',
        'flags_weights': 'Flags and weights',
        'auxiliary_columns': 'Other auxiliary columns'}
    cat_rows = [[label, f"{run['byte_categories'][key]:,}",
                 f"{100*run['byte_categories'][key]/total:.3f}"] for key, label in labels.items()]
    cat_rows.append(['Total', f'{total:,}', '100.000'])
    bit_rows = []
    dist = selected['payload_distribution']
    for bits in ['3','4','5','6','8','10','12','14','16','float32']:
        bit_rows.append(['Exact float32' if bits == 'float32' else bits,
                         f"{perf['block_bits'].get(bits,0):,}",
                         f"{perf['sample_bits'].get(bits,0):,}",
                         f"{100*perf['sample_bits'].get(bits,0)/perf['complex_samples']:.3f}",
                         f"{dist['serialized_bytes_by_bits'].get(bits,0)/1e6:.6f}"])
    pre = audit['optimized_precompression_accounting']
    nblocks = sum(perf['block_bits'].values())
    quantblocks = nblocks - perf['block_bits']['float32']
    assert pre['header_bytes'] == 28*quantblocks + 16*perf['block_bits']['float32']
    storage = table(['Representation', 'Whole MB', 'Visibility MB', 'Scope'], size_rows) + f'''

The full-source-to-archive change includes a deliberate preservation-policy
change: raw DATA and MODEL_DATA each contain {original/1e6:.3f} logical MB and are
omitted; original flagged visibility values are zeroed. Their combined
{2*original/1e6:.3f} logical MB must not be interpreted as separately measured
physical file savings. The exact control uses the same calibrated-only policy,
flags and exact-on-this-dataset four-bit weights, with no visibility quantization.
Its verified shared-file archive occupies {lossless['archive_bytes']/1e6:.3f} MB.
Before lossless structural sharing, its directory bundle occupied
{lossless['storage_compaction']['source_bundle_bytes']/1e6:.3f} MB; sharing and
container compaction save {lossless['saved_bytes']/1e6:.3f} MB without requantization.

Relative to this retained-data lossless control, the optimized 0.18 archive is
{lossless['archive_bytes']/total:.3f}× smaller ({100*total/lossless['archive_bytes']:.3f}% remaining).
Its visibility stream is {100*payload/original:.3f}% of logical complex64 bytes,
or {100*payload/lossless['stored_bytes_by_category']['visibility_lossless']:.3f}% of
the already losslessly compressed control stream. The complete archive is
{100*total/study['original_ms_physical_bytes']:.3f}% of the physical source, but that
ratio combines omission, lossless compression, deduplication and quantization.
The {total-selected['shared']['archive_bytes']:,}-byte increase relative to the earlier
0.18 file is provenance/container overhead; visibility bytes are unchanged.
[Lossless sharing measurements](../../reports/calibrated-shared-storage.json),
[storage audit](../../reports/calibrated-paper-storage.json).

The optimized run-1 archive has the following physical allocation. Every byte is
included, including container overhead; percentages may round independently.

''' + table(['Contributor', 'Bytes', '% of archive'], cat_rows) + f'''

Visibility payload contributes {100*payload/total:.2f}% of the file. Error records
and quantization statistics together occupy
{(run['byte_categories']['visibility_quality']+run['byte_categories']['quantization_statistics'])/1e6:.3f} MB
({100*(run['byte_categories']['visibility_quality']+run['byte_categories']['quantization_statistics'])/total:.2f}%).
These diagnostic arrays, followed by coordinates, are the main remaining costs;
metadata has not disappeared into a visibility-only compression claim.

The 0.18 precision distribution is shown numerically below. Each width applies
to both real and imaginary components; a four-bit component block has eight
nominal code bits per complex sample. Figure 4 compares all three epsilon cases.

''' + table(['Bits/component', 'Blocks', 'Complex samples', '% samples', 'Stored block MB'], bit_rows) + f'''

The {nblocks:,} blocks include {quantblocks:,} quantized and
{perf['block_bits']['float32']:,} exact blocks. All exact blocks in this case contain
positive-zero placeholders: {dist['exact_positive_zero_samples']:,} samples but only
{dist['exact_positive_zero_bytes']/1e6:.6f} MB after lossless compression. This does
not recover discarded flagged values. The remaining {pre['flagged_exception_samples']:,}
flagged placeholders are exact exceptions inside mixed blocks. Thus a large
exact sample fraction does not imply a large exact-byte cost. Widths unused here
remain available to other data.

Before secondary compression, packed integer codes occupy
{pre['packed_code_bytes']/1e6:.6f} MB, exact-block component arrays
{pre['lossless_component_bytes']/1e6:.6f} MB, exception positions
{pre['exception_position_bytes']/1e6:.6f} MB and exception values
{pre['exception_value_bytes']/1e6:.6f} MB. Headers total
{pre['header_bytes']/1e6:.6f} MB, including {16*quantblocks/1e6:.6f} MB for the four
float32 factors (two centre components, range and gamma). Factor bytes are a
subset of headers, not an additional category. These are logical precompression
counts: codes and exceptions share compressed block bodies, so their final
physical contributions cannot be uniquely separated. The measured combined
visibility stream is {payload/1e6:.6f} MB.

![Rate–distortion and physical byte categories](figs/epsd-rate-distortion.png)

*Figure 9. Initial three-epsilon Stokes I study: complete archive size versus
RMS(image − reference) / RMS(reference residual), in percent for both CLEAN runs
(left), and physical byte contributors (right). The 0.18 visibility payload is
identical in the subsequent I/V study. Full optimized-file accounting is given
above; nominal component widths alone do not determine whole-archive size.*
'''
    speed_rows = []
    for key, label in [('baseline','Baseline (1 process × 24 threads)'),('optimized','Optimized (8 processes × 6 threads)')]:
        e = perf[key]
        speed_rows.append([label, ' / '.join(f"{r['process_wall_seconds']:.3f}" for r in e['runs']),
                           f"{e['median_seconds']:.3f}", f"{e['minimum_seconds']:.3f}–{e['maximum_seconds']:.3f}"])
    performance = '''Three unprofiled, fresh-process runs of each implementation measure elapsed
wall time from the complete 15-integration MS to a verified single-file archive.
The boundary includes discovery, reference sampling, inspection, encoding,
integrity and source checks, and final compaction/validation. It excludes MS
export and imaging. Scientific settings, flags/weights policy and decoded values
are held fixed at epsilon-D 0.18; no error bound is relaxed for this comparison.

The host has two Intel Xeon Silver 4214 CPUs, 48 logical CPUs and 376.5 GiB RAM.
Runs use a warm OS cache without eviction. Task-owned timing runs execute
sequentially, separately from profiling, unit tests and CLEAN; unrelated activity
on this shared host is not controlled. Six one-run process/thread allocation
pilots preceded selection of eight processes with six encoder threads each.
The table reports the subsequent three repeats, not the fastest pilot.

''' + table(['Execution', 'Three full-workflow runs, s', 'Median s', 'Min–max s'], speed_rows) + f'''

The median falls from {perf['baseline']['median_seconds']:.3f} to
{perf['optimized']['median_seconds']:.3f} s, a {perf['median_speedup']:.3f}× speedup.
In the optimized repeats, verified-bundle creation takes 53.105–56.132 s and
single-file compaction/validation another 16.592–19.325 s. The complete-MS
15-second target remains unmet. Parallel stage-time sums are not elapsed latency.
This benchmark compares complete execution configurations; it does not isolate
algorithmic changes from the increased worker allocation.

Reported parent peak RSS is 0.991–1.016 GB for the baseline and approximately
0.136 GB for the optimized workflow; the largest individual optimized worker
peak is 0.714–0.738 GB. A simultaneous aggregate memory peak was not measured,
so these values are not a total-memory comparison. Separately, the frozen
independent decoder exported the optimized archive in 58.35 s; export plus
source verification and resampling took 118.08 s. That check overlapped later
imaging and is not a repeated decoder-performance benchmark.

All repeated archives preserve the same {payload:,}-byte visibility payload,
precision distribution and reconstructed values. Their few-hundred-byte total
size variation comes from provenance. The performance validation reports
{perf['tests_passed']} passing tests, independent decoding and exact retained-column checks.
Detailed profiles, software records and reproduction commands are kept in the
[benchmark note](../../devnotes/calibrated-performance.md) and
[measured report](../../reports/calibrated-performance.json).

'''
    performance += runtime_table(timing)
    abstract = '''We present varquanMS, an experimental archive for radio-interferometric
Measurement Sets using local scaling and adaptive integer precision. Each
baseline–frequency block stores a complex median, range, power-spacing parameter
and a selected component width. A dimensionless epsilon-D criterion bounds
reconstruction error relative to local component median deviation, with maximum,
bias and baseline guards. Nine widths from 3 to 16 bits plus exact float32 fit
within a four-bit selector. Original UVW, sampling and retained metadata are
preserved. On one calibrated 15-integration OVRO-LWA observation, epsilon-D
0.01, 0.04 and 0.18 yield visibility payloads of 22.859, 17.433 and 9.932% of
logical complex64 bytes. At 0.18, a complete archive is 242.769 MB versus
985.236 MB for the same-policy retained-data lossless archive. The calibrated
profile explicitly discards raw/model columns and flagged values while retaining
flags; its four-bit weights are exact on this dataset. Additional I/V imaging
above 10° gives difference RMS of 4.810–4.876% and 10.056% of the matching
reference residual RMS, respectively. A follow-up extends V validation to all
three epsilon values. Three repeated full-MS conversion runs
reduce median time from 485.870 to 72.685 s (6.685×), including verification;
the 15-second target remains unmet. These measurements establish agreement with
the supplied calibrated reference, not sky truth or general polarization fidelity.'''
    conclusion = '''The three fixed epsilon-D values produce archives of 242.766–427.429 MB;
the initial Stokes I study has worst-of-two differences of 1.151–5.584% of
reference residual RMS. Epsilon-D 0.01 gives a smaller dirty-image difference
than 0.04 (0.146% versus 0.529%), yet a larger worst restored difference
(1.260% versus 1.151%). Tighter visibility limits therefore do not impose a
monotonic ordering on these CLEAN results. Two repeats do not establish
significance or a causal explanation for this reversal.

The additional 0.18 study measures I differences of 4.810–4.876% and V differences
of 10.056% of their respective reference residual RMS. V is more affected by this
measure than I, and exceeds its observed reference-repeat variation. The
242.769 MB complete archive is 4.058× smaller than the same-policy lossless
control; visibility payload is only 58.56% of the final file, with quality and
reference-statistics arrays contributing a further 26.03%. Full conversion takes
a median 72.685 s, including verification, versus 485.870 s for the baseline.
These are measured size, fidelity and speed tradeoffs for one observation;
none defines a universal optimum or a scientific acceptance threshold.'''
    return {'ABSTRACT': abstract, 'IV': iv, 'STORAGE': storage,
            'PERFORMANCE': performance, 'CONCLUSION': conclusion}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--skip-images',action='store_true',help='Refresh paper text/tables without local FITS files')
    args = parser.parse_args()
    os.chdir(ROOT)
    study = load('calibrated-epsd-dense')
    perf = load('calibrated-performance')
    lossless = next(e for e in load('calibrated-shared-storage')['experiments'] if e['name']=='lossless')
    audit = load('calibrated-paper-storage')
    vstudy = load('calibrated-epsd-v')
    timing = load('calibrated-runtime-breakdown')
    text = PAPER.read_text()
    for name, body in sections(study,perf,lossless,audit,vstudy,timing).items():
        text = replace_section(text,name,body)
    if not args.skip_images:
        from build_epsd_figures import measured, image_comparisons
        measured(study)
        image_comparisons(vstudy)
    PAPER.write_text(text)
    print('Updated paper from frozen reports; scientific experiments were not rerun.')


if __name__ == '__main__':
    main()
