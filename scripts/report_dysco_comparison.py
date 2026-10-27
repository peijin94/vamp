"""Collect verified calibrated Dysco/vamp rate-distortion experiments.

The input MS, FITS files and per-run reports are local, ignored experiment
outputs. This script refuses to report an incomplete or unverified case.
"""
import argparse
import hashlib
import json
from pathlib import Path

from run_wsclean import snapshot
from run_dysco_benchmark import case_name
from vamp.inspection import _open


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def binned_diagnostics(report):
    keys = ('scope','amplitude_edges','length_edges_m','by_amplitude','by_baseline_length','by_channel')
    return {key:report[key] for key in keys}


def ms_storage_categories(directory):
    """Account for physical file lengths by root storage manager and subtables."""
    with _open(directory) as table:
        managers = list(table.getdminfo().values())
    files = [p for p in directory.rglob('*') if p.is_file()]
    root_files = [p for p in files if p.parent == directory]
    claimed = set()
    result = []
    for manager in managers:
        prefix = f"table.f{manager['SEQNR']}"
        matching = [p for p in root_files if p.name == prefix or p.name.startswith(prefix+'_')]
        claimed.update(matching)
        result.append({'type':manager['TYPE'],'columns':list(manager['COLUMNS']),
                       'bytes':sum(p.stat().st_size for p in matching),
                       'files':[p.name for p in matching]})
    subtables = sum(p.stat().st_size for p in files if p.parent != directory)
    other_root = sum(p.stat().st_size for p in root_files if p not in claimed)
    return {'root_storage_managers':result,'subtables_bytes':subtables,
            'other_root_bytes':other_root,
            'total_bytes':sum(entry['bytes'] for entry in result)+subtables+other_root}


def image_metrics(comparison, final_threshold=None):
    if comparison.get('status') != 'measured':
        raise ValueError('Image comparison is incomplete')
    result = comparison['residual_normalized']['image-image.fits']
    if result['denominator_status'] != 'positive':
        raise ValueError('Reference residual normalization is undefined')
    dirty = comparison['products']['image-dirty.fits']['regions']['above_elevation_cut']['difference_rms']
    return {'difference_rms_jy_per_beam': result['difference_rms'],
            'reference_residual_rms_jy_per_beam': result['reference_residual_rms'],
            'difference_over_reference_residual_percent': 100 * result['difference_over_reference_residual_rms'],
            'dirty_difference_over_reference_residual_percent': 100 * dirty / result['reference_residual_rms'],
            'clean_final_threshold_reached': comparison.get('final_threshold_reached',final_threshold)}


def collect(root, normalizations, bits, rf15_bits, flagged_input, af_bits, zero_bits):
    source_i = Path('reports/calibrated-epsd-dense.json')
    source_v = Path('reports/calibrated-epsd-v.json')
    performance_path = Path('reports/calibrated-performance.json')
    original_i, original_v, performance = read(source_i), read(source_v), read(performance_path)
    identity = read(root / 'input-identity.json')
    if (identity['status'] != 'passed' or original_i['status'] != 'complete'
            or original_v['status'] != 'passed' or performance['status'] != 'complete'):
        raise ValueError('Reference and input verification must pass')
    for pol, expected in (('I',original_i['reference_residual_statistics']['rms']),
                          ('V',original_v['reference_residual_statistics']['rms'])):
        if performance['reference_repeat_controls'][pol]['image-image.fits']['reference_residual_rms'] != expected:
            raise ValueError(f'{pol} reference repeat uses a different residual normalization')
    if snapshot(Path(identity['original_ms'])) != identity['original_sha256']:
        raise ValueError('Original source MS changed')
    if snapshot(Path(identity['prepared_ms'])) != identity['prepared_sha256']:
        raise ValueError('Prepared calibrated MS changed')
    if sha(identity['verification']) != identity['verification_sha256']:
        raise ValueError('Prepared input verification report changed')
    calibrated_input_checks = read(identity['verification'])['column_checks']
    independent_path = root/'input-identity-independent.json'
    independent = read(independent_path) if independent_path.exists() else None
    if independent is not None and (independent.get('status')!='passed'
            or independent['original_sha256'] != identity['original_sha256']
            or independent['prepared_sha256'] != identity['prepared_sha256']):
        raise ValueError('Independent calibrated-input recheck does not match the imaged sources')
    mask = original_v['comparison_mask']
    if sha(mask['path']) != mask['sha256']:
        raise ValueError('Imaging comparison mask changed')
    v_by_name = {case['name']: case for case in original_v['experiments']}
    cases = []
    for case in original_i['experiments']:
        vcase = v_by_name[case['name']]
        if case['export_verification']['status'] != 'passed' or not case['weights_all_exact']:
            raise ValueError('vamp export or weights verification did not pass')
        if not all(r['clean_termination']['final_threshold_reached'] for r in case['residual_normalized_runs']):
            raise ValueError('Stokes I CLEAN did not reach its threshold')
        if not all(r['final_threshold_reached'] and r['source_unchanged'] for r in vcase['runs']):
            raise ValueError('Stokes V CLEAN or source verification did not pass')
        i_runs = [{'difference_rms_jy_per_beam': r['restored_difference']['rms'],
                   'reference_residual_rms_jy_per_beam': original_i['reference_residual_statistics']['rms'],
                   'difference_over_reference_residual_percent': 100*r['restored_difference_over_reference_residual'],
                   'dirty_difference_over_reference_residual_percent': 100*r['dirty_difference_over_reference_residual'],
                   'clean_final_threshold_reached':r['clean_termination']['final_threshold_reached']}
                  for r in case['residual_normalized_runs']]
        v_runs = [image_metrics(r['comparison'],r['final_threshold_reached']) for r in vcase['runs']]
        cases.append({'method': 'vamp', 'name': case['name'], 'epsilon_D': case['epsilon_D'],
                      'visibility_payload_bytes': case['rate']['payload_bytes'],
                      'visibility_with_permutations_and_indexes_bytes': case['rate']['payload_plus_permutations_and_indexes_bytes'],
                      'complete_output_bytes': case['shared']['archive_bytes'],
                      'stored_bytes_by_category': case['shared']['stored_bytes_by_category'],
                      'precision_sample_counts': case['sample_bits'],
                      'precision_block_counts': case['block_bits'],
                      'payload_distribution':case['payload_distribution'],
                      'visibility_products': case['visibility']['products'],
                      'visibility_binned_diagnostics': binned_diagnostics(case['visibility']),
                      'stokes_I': i_runs, 'stokes_V': v_runs,
                      'config': case['config'], 'export_verified': True})
    settings = [(n, b, 2.5, flagged_input) for n in normalizations for b in bits]
    settings += [('RF', b, 1.5, flagged_input) for b in rf15_bits]
    settings += [('AF', b, 2.5, 'nan') for b in af_bits]
    settings += [('RF', b, 2.5, 'zero') for b in zero_bits]
    if len(set(settings)) != len(settings):
        raise ValueError('Duplicate experiment settings')
    for normalization, bit, truncation, flags in settings:
        name = case_name(normalization, bit, truncation, flags)
        enc_path, process_path = root/(name+'-encode.json'), root/(name+'-process.json')
        diagnostic_path = root/(name+'-visibility.json')
        audit_path = root/(name+'-verification-v2.json')
        enc, process, diagnostic, audit = map(read, (enc_path,process_path,diagnostic_path,audit_path))
        if (enc['status'] != 'complete' or process['returncode'] != 0 or audit['status'] != 'passed'
                or enc['bits_per_real_component'] != bit or enc['normalization'] != normalization
                or enc['distribution_truncation'] != truncation
                or enc['flagged_input_values'] != ('complex NaN' if flags == 'nan' else 'positive complex zero')):
            raise ValueError(f'{name}: encoding or verification failed')
        if diagnostic['source_column'] != 'CORRECTED_DATA' or diagnostic['decoded_column'] != 'DATA':
            raise ValueError(f'{name}: wrong visibility columns')
        if (Path(audit['source']).resolve() != Path(identity['prepared_ms']).resolve()
                or Path(audit['decoded']).resolve() != (root/(name+'.ms')).resolve()
                or Path(diagnostic['source']).resolve() != Path(audit['source']).resolve()
                or Path(diagnostic['decoded']).resolve() != Path(audit['decoded']).resolve()):
            raise ValueError(f'{name}: visibility or audit source does not match the verified MS')
        imaging = {}
        for pol in ('i','v'):
            runs = []
            for suffix in ('','-repeat'):
                path = root/(name+'-'+pol+suffix+'-comparison.json')
                if suffix and not path.exists():
                    continue
                comp = read(path)
                if (not isinstance(comp.get('final_threshold_reached'),bool) or comp['comparison_mask']['sha256'] != mask['sha256']
                        or comp['dysco_verification_sha256'] != sha(audit_path)
                        or comp['input_identity_sha256'] != sha(root/'input-identity.json')):
                    raise ValueError(f'{name}: imaging threshold or mask mismatch')
                runs.append(image_metrics(comp))
            imaging['stokes_'+pol.upper()] = runs
        v3_path = root/(name+'-verification-v3.json')
        bit_audit = read(v3_path) if v3_path.exists() else audit
        if bit_audit['decoded_sha256'] != audit['decoded_sha256']:
            raise ValueError(f'{name}: verification snapshots disagree')
        storage_categories = ms_storage_categories(root/(name+'.ms'))
        if storage_categories['total_bytes'] != enc['storage']['total_bytes']:
            raise ValueError(f'{name}: MS file sizes changed after encoding')
        cases.append({'method':'Dysco', 'name':name, 'bits_per_component':bit,
                      'normalization':normalization, 'distribution':enc['distribution'],
                      'flagged_input':flags,
                      'distribution_truncation':enc['distribution_truncation'],
                      'visibility_manager_bytes':enc['storage']['dysco_data_manager_bytes'],
                      'complete_output_bytes':enc['storage']['total_bytes'],
                      'allocated_output_bytes':enc['storage']['allocated_bytes'],
                      'stored_bytes_by_category':storage_categories,
                      'visibility_products':diagnostic['products'],
                      'visibility_binned_diagnostics':binned_diagnostics(diagnostic), **imaging,
                      'audit':{'table_count':audit['table_count'],
                               'exact_column_count':audit['exact_column_count'],
                               'exact_columns':audit['exact_columns'],
                               'flagged_samples':audit['flagged_samples'],
                                   'flagged_input_policy':audit.get('flagged_input_policy',enc['flagged_input_values']),
                                   'flagged_values_policy_mismatch':audit.get('flagged_values_policy_mismatch',audit['flagged_values_not_zero']),
                               'flagged_values_not_zero':audit['flagged_values_not_zero'],
                               'flagged_values_not_nan':audit.get('flagged_values_not_nan'),
                               'flagged_noncanonical_zero_components':bit_audit.get('flagged_noncanonical_zero_components'),
                               'unflagged_nonfinite_values':audit['unflagged_nonfinite_values']},
                      'timing':{'process_wall_seconds':process['process_wall_seconds'],
                                'work_seconds':enc['work_seconds'],
                                'metadata_copy_seconds':enc['metadata_copy_seconds'],
                                'putcol_and_flush_seconds':enc['putcol_and_flush_seconds'],
                                'close_seconds':enc['close_seconds'],
                                'peak_rss_bytes':enc['peak_rss_bytes'],
                                'scope':process['timing_scope']},
                      'artifact_sha256':{str(p):sha(p) for p in (enc_path,process_path,diagnostic_path,audit_path)}})
    return {'schema':'vamp.dysco-comparison/1', 'status':'complete',
            'measurement_date':'2026-10-27', 'source':'15-integration calibrated OVRO-LWA MS',
            'source_complex64_bytes':original_i['original_complex64_bytes'],
            'source_physical_bytes_including_raw_and_model':original_i['original_ms_physical_bytes'],
            'mask':mask, 'container_id':original_v['container_id'],
            'policy':'Calibrated visibilities only; original raw DATA/MODEL_DATA omitted; flagged visibility values purged, FLAG and FLAG_ROW retained, WEIGHT and SIGMA exact on this input. Primary Dysco runs encode flagged values as NaN; vamp uses zero.',
            'visibility_rate_definition':'vamp payload plus permutations and visibility indexes versus Dysco DATA manager file including payload, scales and manager header; generic container/table overhead is excluded in both.',
            'complete_rate_definition':'Complete vamp HDF5 archive versus complete Dysco Measurement Set; different containers and metadata overhead are included.',
            'image_metric':'100 × RMS(candidate restored − reference restored) / RMS(reference residual), over the fixed mask above 10 degrees; not a thermal-noise guarantee.',
            'uncompressed_reference_repeat_percent':{
                pol:100*performance['reference_repeat_controls'][pol]['image-image.fits']['difference_over_reference_residual_rms']
                for pol in ('I','V')},
            'all_dysco_clean_final_thresholds_reached':all(
                run['clean_final_threshold_reached'] for case in cases if case['method']=='Dysco'
                for pol in ('stokes_I','stokes_V') for run in case[pol]),
            'source_identity_verified_at_report_time':True,
            'calibrated_input_column_checks':calibrated_input_checks,
            'independent_input_recheck':None if independent is None else {
                'exact_columns':independent['independent_exact_column_checks'],
                'report_sha256':sha(independent_path),
                'lossless_archive_sha256':independent['lossless_archive_sha256']},
            'provenance_sha256':{str(p):sha(p) for p in (source_i,source_v,performance_path,root/'input-identity.json')},
            'timing_limitations':'Dysco process timings exclude calibrated input preparation, subsequent verification and imaging; prior vamp verified archive timings have a different scope. NaN-flagged RF, RF 1.5 and AF runs overlap other jobs; AF containers are CPU-limited. No controlled cold-cache comparison.',
            'cases':cases}


def plots(report, root, destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.io import fits
    destination.mkdir(parents=True, exist_ok=True)
    styles = {'vamp':('o','#1d668c'), 'RF':('s','#c1512c'),
              'RF1.5':('D','#d59921'), 'RFzero':('x','#888888'), 'AF':('^','#8b5da7')}
    adaptive = next(case for case in report['cases'] if case.get('epsilon_D') == .18)
    sample_counts = adaptive['precision_sample_counts']
    if sample_counts.get('float32',0) != adaptive['payload_distribution']['exact_positive_zero_samples']:
        raise ValueError('Exact-zero and float32 block sample counts are not interchangeable')
    colours = {'3':'#d88946','4':'#3d8ab0','5':'#66ad99','6':'#e4b85a',
               '8':'#a386bd','10':'#8eaf4b','12':'#cf87a7','14':'#b88454',
               '16':'#6b73a5','float32':'#aaaaaa'}
    if set(sample_counts)-set(colours):
        raise ValueError('Unknown precision in the plotted experiment')
    total = sum(sample_counts.values())
    if total*8 != report['source_complex64_bytes']:
        raise ValueError('Precision counts do not cover the original visibility samples')
    fig,ax=plt.subplots(figsize=(10,2.9),constrained_layout=True)
    offset=0.
    for bit in ('3','4','5','6','8','10','12','14','16','float32'):
        n=sample_counts.get(bit,0)
        if not n:
            continue
        width=100*n/total
        label='exact-zero blocks' if bit=='float32' else f'{bit} bits/component'
        ax.barh(1,width,left=offset,color=colours[bit],label=label,height=.55)
        if width >= 3:
            ax.text(offset+width/2,1,f'{width:.1f}%',ha='center',va='center',fontsize=10)
        offset+=width
    ax.barh(0,100,color=colours['3'],height=.55)
    ax.text(50,0,'100% at 3 bits/component',ha='center',va='center',fontsize=10)
    ax.set(yticks=[1,0],yticklabels=['vamp εD=0.18','Dysco RF 3-bit'],
           xlim=(0,100),xlabel='Share of input complex samples (%)')
    ax.set_title('Sample counts by selected block precision; includes flagged sample slots')
    ax.legend(ncol=4,loc='lower center',bbox_to_anchor=(.5,-.5),fontsize=8)
    fig.text(.5,-.1,'Exact-zero blocks are stored compactly; these percentages are not byte shares.',ha='center',fontsize=9)
    fig.savefig(destination/'dysco-matched-bit-distribution.png',dpi=170,bbox_inches='tight');plt.close(fig)
    for scope in ('complete','visibility'):
        fig, axes = plt.subplots(1,3,figsize=(14.5,4.4),constrained_layout=True)
        for case in report['cases']:
            if case.get('flagged_input') == 'zero':
                continue  # Retained as a control in JSON; omit overlapping controls from the plot.
            key = case.get('normalization','vamp')
            if key == 'RF':
                if case.get('flagged_input') == 'zero':
                    key = 'RFzero'
                elif case['distribution_truncation'] == 1.5:
                    key = 'RF1.5'
            marker,color = styles[key]
            label = ('εD='+str(case['epsilon_D']) if key == 'vamp' else str(case['bits_per_component'])+' bit')
            field = ('complete_output_bytes' if scope=='complete' else
                     'visibility_with_permutations_and_indexes_bytes' if key=='vamp' else 'visibility_manager_bytes')
            x = case[field]/1e6
            for ax,field in zip(axes[:2],('stokes_I','stokes_V')):
                y = case[field][0]['difference_over_reference_residual_percent']
                ax.scatter(x,y,marker=marker,color=color,s=45,zorder=3)
                plotted_label = label + ('*' if case[field][0]['clean_final_threshold_reached'] is False else '')
                if key != 'RF1.5':
                    ax.annotate(plotted_label,(x,y),xytext=(4,4),textcoords='offset points',fontsize=8)
            y = 100*case['visibility_products']['cross/I-paired']['relative_rmse']
            axes[2].scatter(x,y,marker=marker,color=color,s=45,zorder=3)
            if key != 'RF1.5':
                axes[2].annotate(label,(x,y),xytext=(4,4),textcoords='offset points',fontsize=8)
        for ax,title in zip(axes,('Stokes I CLEAN','Stokes V CLEAN','Paired cross-baseline Stokes I visibility')):
            xlabel = 'Complete output size (MB)' if scope=='complete' else 'Visibility payload and decoding metadata (MB)'
            ax.set_title(title); ax.set_xlabel(xlabel); ax.set_yscale('log'); ax.grid(alpha=.25)
        axes[0].set_ylabel('Restored image difference / reference residual RMS (%)')
        axes[1].set_ylabel('Restored image difference / reference residual RMS (%)')
        axes[2].set_ylabel('Cross-baseline visibility relative RMS error (%)')
        present = {'RFzero' if x.get('normalization')=='RF' and x.get('flagged_input')=='zero'
                   else 'RF1.5' if x.get('normalization')=='RF' and x.get('distribution_truncation')==1.5
                   else x.get('normalization','vamp') for x in report['cases']
                   if x.get('flagged_input')!='zero'}
        labels = {'vamp':'vamp','RF':'Dysco RF, 2.5σ',
                  'RF1.5':'Dysco RF, 1.5σ','RFzero':'Dysco RF, zero control','AF':'Dysco AF, 2.5σ'}
        handles = [plt.Line2D([],[],linestyle='',marker=m,color=c,label=labels[k])
                   for k,(m,c) in styles.items() if k in present]
        fig.legend(handles=handles,loc='lower center',ncol=len(handles),bbox_to_anchor=(.5,-.08))
        title = 'Calibrated 15-integration MS: output size and measured error'
        if not report['all_dysco_clean_final_thresholds_reached']:
            title += '\n* CLEAN did not reach the configured final threshold'
        fig.suptitle(title,fontsize=12)
        suffix = '' if scope=='complete' else '-visibility'
        fig.savefig(destination/f'dysco-rate-distortion{suffix}.png',dpi=170,bbox_inches='tight'); plt.close(fig)

    mask = np.squeeze(fits.getdata(report['mask']['path'])).astype(bool)
    zoom = read('reports/calibrated-v-zoom.json')
    if (sha(zoom['reference_fits']) != zoom['reference_fits_sha256']
            or zoom['comparison_mask_sha256'] != report['mask']['sha256']
            or zoom['entire_crop_minimum_elevation_deg'] <= 45):
        raise ValueError('The reference V zoom selection does not match the comparison')
    bounds = zoom['bounds']
    view = (slice(bounds['y_start'],bounds['y_stop']),slice(bounds['x_start'],bounds['x_stop']))
    if not mask[view].all():
        raise ValueError('V zoom contains pixels outside the comparison mask')
    included = {case['name']:case for case in report['cases']}
    variants = []
    if 'dysco-rf-03-nan' in included:
        variants.append(('dysco-rf-03-nan','','RF',2.5))
    elif 'dysco-rf-03' in included:
        variants.append(('dysco-rf-03','','RF',2.5))
    if 'dysco-rf-03-t15-nan' in included:
        variants.append(('dysco-rf-03-t15-nan','-rf15','RF',1.5))
    if 'dysco-af-03-nan' in included:
        variants.append(('dysco-af-03-nan','-af','AF',2.5))
    ours_bytes = next(c['visibility_with_permutations_and_indexes_bytes'] for c in report['cases']
                      if c.get('epsilon_D') == 0.18)
    for stem,suffix,normalization,truncation in variants:
        pairs = [('I',Path('results/calibrated-template/reference'),
                  Path('results/calibrated-epsd-dense/epsd-018-dense-images'),
                  root/(stem+'-i-images')),
                 ('V',Path('results/calibrated-template/reference-v'),
                  Path('results/calibrated-epsd-dense/epsd-018-dense-v-images'),
                  root/(stem+'-v-images'))]
        fig, axes = plt.subplots(2,3,figsize=(13,8),constrained_layout=True)
        for row,(pol,ref,ours,dysco) in enumerate(pairs):
            original = np.squeeze(fits.getdata(ref/'image-image.fits')).astype(np.float64)
            a = np.squeeze(fits.getdata(ours/'image-image.fits')).astype(np.float64)
            b = np.squeeze(fits.getdata(dysco/'image-image.fits')).astype(np.float64)
            residual = np.squeeze(fits.getdata(ref/'image-residual.fits')).astype(np.float64)
            denominator = np.sqrt(np.mean(residual[mask]**2))
            limit = 2.5 if pol == 'V' else np.percentile(original[mask],99.7)
            refimg=axes[row,0].imshow(np.where(mask,original,np.nan),origin='lower',
                                     cmap='bwr' if pol=='V' else 'viridis',
                                     vmin=-limit if pol=='V' else 0,vmax=limit)
            axes[row,0].set_title(f'{pol} reference');fig.colorbar(refimg,ax=axes[row,0],shrink=.8,label='Jy/beam')
            for col,(image,title) in enumerate(((a,'vamp εD=0.18'),(b,f'Dysco {normalization} 3-bit, {truncation}σ')),1):
                difference=100*(image-original)/denominator
                rms=np.sqrt(np.mean(difference[mask]**2))
                panel=axes[row,col].imshow(np.where(mask,difference,np.nan),origin='lower',cmap='bwr',vmin=-30,vmax=30)
                axes[row,col].set_title(f'{title} − reference\nRMS = {rms:.2f}%')
                fig.colorbar(panel,ax=axes[row,col],shrink=.8,label='% of reference residual RMS',extend='both')
            for ax in axes[row]: ax.set_xticks([]); ax.set_yticks([])
        dysco_bytes = included[stem]['visibility_manager_bytes']
        fig.suptitle(f'Visibility storage: {ours_bytes/1e6:.2f} MB vs {dysco_bytes/1e6:.2f} MB; pixels above 10°',fontsize=12)
        fig.savefig(destination/f'dysco-matched-images{suffix}.png',dpi=150,bbox_inches='tight'); plt.close(fig)

        images = (original[view],a[view],b[view])  # last row above is Stokes V
        fig,axes = plt.subplots(2,3,figsize=(11,7),constrained_layout=True)
        titles = ('Reference','vamp εD=0.18',f'Dysco {normalization} 3-bit, {truncation}σ')
        for col,(v,title) in enumerate(zip(images,titles)):
            upper=axes[0,col].imshow(v,origin='lower',cmap='bwr',vmin=-8,vmax=8)
            axes[0,col].set_title(title)
            delta=100*(v-images[0])/denominator
            lower=axes[1,col].imshow(delta,origin='lower',cmap='bwr',vmin=-30,vmax=30)
            axes[1,col].set_title(f'Crop difference RMS = {np.sqrt(np.mean(delta**2)):.2f}%')
            for row in (0,1): axes[row,col].set_xticks([]); axes[row,col].set_yticks([])
        fig.colorbar(upper,ax=axes[0,:],shrink=.75,label='Stokes V (Jy/beam)')
        fig.colorbar(lower,ax=axes[1,:],shrink=.75,label='% of full-domain reference residual RMS',extend='both')
        fig.suptitle('Taurus A region: Stokes V; all crop pixels above 45°',fontsize=12)
        fig.savefig(destination/f'dysco-matched-v-zoom{suffix}.png',dpi=170,bbox_inches='tight'); plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path('results/dysco-comparison'))
    parser.add_argument('--normalizations',nargs='+',choices=('RF','AF'),required=True)
    parser.add_argument('--bits',nargs='+',type=int,default=[3,4,6,8])
    parser.add_argument('--rf15-bits',nargs='*',type=int,default=[])
    parser.add_argument('--flagged-input',choices=('zero','nan'),default='nan')
    parser.add_argument('--af-bits',nargs='*',type=int,default=[])
    parser.add_argument('--zero-bits',nargs='*',type=int,default=[])
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--figures',type=Path,required=True)
    args=parser.parse_args()
    result=collect(args.root,args.normalizations,args.bits,args.rf15_bits,args.flagged_input,args.af_bits,args.zero_bits)
    plots(result,args.root,args.figures)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'cases':len(result['cases']),
                      'report':str(args.output),'figures':str(args.figures)}),flush=True)


if __name__ == '__main__':
    main()
