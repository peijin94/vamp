"""Audit new three-epsilon V CLEANs and the existing median-run timing breakdown."""
import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from vamp.imaging import file_sha256, read_fits, _grid

ROOT = Path('results/calibrated-epsd-dense')
REF = Path('results/calibrated-template/reference-v')


def load(path):
    return json.loads(Path(path).read_text())


def timing_breakdown():
    perf_path = Path('reports/calibrated-performance.json')
    perf = load(perf_path)
    timing_root = Path('results/calibrated-performance/optimized')
    processes = load(timing_root/'processes.json')
    process = sorted((p for p in processes if not p['profiled']), key=lambda p:p['process_wall_seconds'])[1]
    raw_path = timing_root/(process['name']+'.json')
    raw = load(raw_path)
    total = process['process_wall_seconds']
    if process['returncode'] or not np.isclose(total, perf['optimized']['median_seconds'], rtol=1e-14):
        raise ValueError('Expected the successful recorded median run')
    compact = raw['compact']
    internal = compact['storage_compaction']
    sampling = raw['bundle']['global_reference_sampling_seconds']
    phases = [
        ('Global reference sampling', sampling, 'bundle.global_reference_sampling_seconds'),
        ('Remaining bundle work: discovery, parallel inspection/encoding/source checks, finalization',
         raw['bundle_wall_seconds']-sampling, 'bundle_wall_seconds minus reference sampling'),
        ('Validate input bundle for compaction', internal['source_validation_seconds'],
         'compact.storage_compaction.source_validation_seconds'),
        ('Copy/share arrays into final archive', internal['copy_and_share_seconds'],
         'compact.storage_compaction.copy_and_share_seconds'),
        ('Validate final archive', compact['output_validation_seconds'], 'compact.output_validation_seconds'),
        ('Other compaction/finalization overhead', compact['elapsed_seconds']-internal['source_validation_seconds']-
         internal['copy_and_share_seconds']-compact['output_validation_seconds'], 'compaction elapsed minus the three timed subphases'),
        ('Process startup/shutdown and remaining orchestration', total-raw['bundle_wall_seconds']-compact['elapsed_seconds'],
         'process wall minus bundle wall and compaction elapsed'),
    ]
    assert all(v>=0 for _,v,_ in phases) and np.isclose(sum(v for _,v,_ in phases),total,rtol=1e-14)
    manifest_path = timing_root/(process['name']+'.bundle/manifest.json')
    manifest = load(manifest_path)
    verification_paths = [manifest_path.parent/p['verification'] for p in manifest['partitions']]
    verifications = [load(p) for p in verification_paths]
    assert all(v['status']=='passed' for v in verifications)
    parallel = {
        'Inspection':raw['bundle']['sum_inspection_seconds'],
        'Encoding and writing, including integrity finalization':raw['bundle']['sum_encode_seconds_excluding_inspection'],
        'Source/integrity verification':sum(v['elapsed_seconds'] for v in verifications),
    }
    return dict(schema='vamp.median-runtime-breakdown/1',status='passed',
        measurement_date='2026-10-25',report_date='2026-10-27',selected_run=process['name'],
        total_process_wall_seconds=total,workflow_seconds=raw['total_workflow_seconds'],
        bundle_wall_seconds=raw['bundle_wall_seconds'],compaction_wall_seconds=compact['elapsed_seconds'],
        phase_wall_seconds=[dict(phase=name,seconds=value,percent=100*value/total,derivation=source) for name,value,source in phases],
        parallel_partition_stage_seconds=parallel,
        scope='Complete 15-integration MS to verified single-file archive; existing unprofiled median run. No new timing experiment.',
        limitations=['Remaining bundle work includes overlapping partition stages and discovery/finalization; a finer additive wall-time split was not instrumented.',
                    'Parallel partition-stage times are sums across 15 integrations, not additive contributions to the 72.685-second latency.',
                    'The two remainder rows are differences of measured elapsed durations, not independent function profiles.'],
        input_sha256={str(p):file_sha256(p) for p in [perf_path,raw_path,timing_root/'processes.json',manifest_path,*verification_paths]})


def v_results():
    original = load('reports/calibrated-epsd-dense.json')
    mask_info = original['comparison_mask']
    mask = fits.getdata(mask_info['path']).astype(bool)
    assert file_sha256(mask_info['path'])==mask_info['sha256']
    reference, header = read_fits(REF/'image-image.fits')
    residual, rh = read_fits(REF/'image-residual.fits')
    assert _grid(header)==_grid(rh) and mask.shape==reference.shape and int(mask.sum())==mask_info['pixels']
    assert np.all(np.isfinite(reference[mask])) and np.all(np.isfinite(residual[mask]))
    values = residual[mask]
    reference_stats = dict(pixels=int(values.size),rms=float(np.sqrt(np.mean(values**2))),
                           mean=float(values.mean()),standard_deviation=float(values.std()))
    ref_complete = load(REF/'complete.json')
    assert ref_complete['status']=='complete' and ref_complete['source_unchanged']
    expected = ref_complete['parameters']
    cases = []
    for source_case in original['experiments']:
        name = source_case['name']
        runs=[]
        export_path = ROOT/(name+'-export.json')
        export = load(export_path)
        assert export['status']=='passed' and export['source_files_unchanged']
        i_complete_path = ROOT/(name+'-images/complete.json')
        i_complete = load(i_complete_path)
        config = load(ROOT/(name+'-config.json'))
        assert config==source_case['config']
        for repeat in (False,True):
            stem=name+'-v'+('-repeat' if repeat else '')
            directory=ROOT/(stem+'-images')
            comparison_path=ROOT/(stem+'-comparison.json')
            comparison=load(comparison_path)
            complete=load(directory/'complete.json')
            assert comparison['status']=='measured' and complete['status']=='complete' and complete['source_unchanged']
            assert complete['source_sha256']==i_complete['source_sha256']
            assert complete['container_id']==ref_complete['container_id']
            for key,value in expected.items():
                if key!='restoring_beam': assert complete['parameters'][key]==value,(name,key)
            assert comparison['comparison_mask']['sha256']==mask_info['sha256']
            for product in comparison['products'].values():
                for label in ('reference','candidate'):
                    assert file_sha256(product[label])==product[label+'_file_sha256']
            image=comparison['products']['image-image.fits']
            normalization=comparison['residual_normalized']['image-image.fits']
            assert np.isclose(normalization['reference_residual_rms'],reference_stats['rms'],rtol=1e-12)
            log=(directory/'wsclean.log').read_text()
            final='Cleaning finished because the final threshold was reached.' in log
            if not final:raise ValueError('CLEAN did not reach final threshold: '+stem)
            runs.append(dict(run='repeat' if repeat else 'first',comparison=comparison,
                comparison_file=str(comparison_path),comparison_sha256=file_sha256(comparison_path),
                normalization=normalization,final_threshold_reached=final,
                wsclean_seconds=complete['wsclean_elapsed_seconds'],source_unchanged=complete['source_unchanged'],
                source_ms=complete['source_ms'],source_sha256=complete['source_sha256'],
                complete_file=str(directory/'complete.json'),complete_sha256=file_sha256(directory/'complete.json')))
        assert runs[0]['source_sha256']==runs[1]['source_sha256']
        cases.append(dict(name=name,epsilon_D=source_case['epsilon_D'],config=config,
            export_verification_file=str(export_path),export_verification_sha256=file_sha256(export_path),
            initial_stokes_i_complete_sha256=file_sha256(i_complete_path),same_input_ms_as_initial_stokes_i=True,
            comparison=runs[0]['comparison'],runs=runs))
    return dict(schema='vamp.epsd-v-followup/1',status='passed',measurement_date='2026-10-27',
        results_root=str(ROOT),reference_directory=str(REF),stokes='V',comparison_mask=mask_info,
        reference_residual_statistics=reference_stats,
        reference_file_sha256=file_sha256(REF/'image-image.fits'),
        reference_residual_file_sha256=file_sha256(REF/'image-residual.fits'),
        reference_complete_sha256=file_sha256(REF/'complete.json'),
        container_id=ref_complete['container_id'],container_digests=ref_complete['container_digests'],
        imaging_parameters=expected,experiments=cases,
        initial_stokes_i_report_sha256=file_sha256('reports/calibrated-epsd-dense.json'),
        limitations=['Follow-up V imaging of existing verified quantized exports; no new encoding or threshold tuning.',
                    'All three cases have two new CLEAN runs; these are distinct from the earlier epsilon-D=0.18 performance-study V runs.',
                    'Differences measure agreement with the supplied calibrated reference, not absolute polarization accuracy.'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--timings-only',action='store_true')
    args=parser.parse_args()
    timing=timing_breakdown()
    Path('reports/calibrated-runtime-breakdown.json').write_text(json.dumps(timing,indent=2,allow_nan=False)+'\n')
    if not args.timings_only:
        result=v_results()
        Path('reports/calibrated-epsd-v.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        for e in result['experiments']:
            print(e['epsilon_D'],[100*r['normalization']['difference_over_reference_residual_rms'] for r in e['runs']])
    print('Median-run additive wall time:',timing['total_process_wall_seconds'])


if __name__=='__main__':main()
