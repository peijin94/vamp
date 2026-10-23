"""Aggregate the 15-slot calibrated-only archives, exports and matched images."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from run_wsclean import snapshot
from varquanms.container import Reader
from varquanms.imaging import read_fits, sky_mask, plot_comparison

ROOT = Path(__file__).resolve().parents[1]


def summarize(kind, root):
    bundle = root / (kind+'-calib.bundle')
    manifest = json.loads((bundle/'manifest.json').read_text())
    summary = json.loads((ROOT/'reports'/('calibrated-'+kind+'-encode.json')).read_text())
    bits, samples, gammas, logical = Counter(), Counter(), Counter(), Counter()
    products = {}; weights = Counter(); weight_bytes = Counter()
    for part in manifest['partitions']:
        s = part['summary']
        for column, q in s.get('weights_quality', {}).items():
            if not q['original_bits_exact']: raise ValueError('This experiment expects exact three-level weights')
            weights.update({k:q[k] for k in ('four_bit_chunks','float32_chunks','lossy_chunks','packed_code_bytes')})
        for spec in s['visibility_block_statistics'].values():
            for product in spec['products']:
                bits.update(product['component_bit_width_histogram'])
                samples.update(product['complex_sample_bit_width_histogram'])
                gammas.update(product.get('power_gamma_block_histogram', {}))
                logical.update(product['logical_bytes_before_secondary_compression'])
        for q in s['visibility_quality'].values():
            for name, values in q['by_product'].items():
                node = products.setdefault(name, Counter())
                n = values['valid_samples']; node['valid_samples'] += n
                for field in ('component_mse','real_mse','imag_mse','bias_real','bias_imag'):
                    node[field] += values[field]*n
                for field in ('max_complex_error','max_baseline_component_mse'):
                    node[field] = max(node[field], values[field])
        with Reader(bundle/part['archive'], verify=False) as reader:
            if part is manifest['partitions'][0]: environment = reader.manifest['environment']
            main = next(t for t in reader.manifest['tables'] if t['path']=='.')
            for name in ('WEIGHT','WEIGHT_SPECTRUM'):
                if name not in main['columns']: continue
                for chunk in main['columns'][name]['chunks']:
                    for key in ('codebook','codes','values'):
                        if key in chunk:
                            entry = reader.entries[chunk[key]]
                            weight_bytes[key+'_stored_bytes'] += entry['stored_bytes']
                            weight_bytes[key+'_logical_bytes'] += entry['logical_bytes']
    for node in products.values():
        n = node['valid_samples']
        for field in ('component_mse','real_mse','imag_mse','bias_real','bias_imag'): node[field] /= n
        node['complex_bias_abs'] = float(np.hypot(node['bias_real'], node['bias_imag']))
    n = sum(q['valid_samples'] for q in products.values())
    mse = sum(q['component_mse']*q['valid_samples'] for q in products.values())/n if n else 0.
    return {'name':kind, 'config':manifest['config'], 'encoding':summary, 'environment':environment,
            'by_product':products, 'component_mse':mse, 'component_bit_width_histogram':dict(bits),
            'complex_sample_bit_width_histogram':dict(samples), 'gamma_histogram':dict(gammas),
            'logical_bytes_before_secondary_compression':dict(logical),
            'weights':dict(weights, **weight_bytes, original_bits_exact=True),
            'selected_visibility_logical_bytes':sum(p['summary']['selected_visibility_logical_bytes'] for p in manifest['partitions']),
            'visibility_dataset_stored_bytes':sum(p['summary']['visibility_datasets_stored_bytes'] for p in manifest['partitions']),
            'partition_times_seconds':[p['summary']['elapsed_seconds'] for p in manifest['partitions']]}


def plot(report):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['svg.hashsalt'] = 'varquanms-calibrated-v1'
    records = report['experiments']; fig, axes = plt.subplots(1,2,figsize=(10,4.4),constrained_layout=True)
    names = ['Calibrated-only\nlossless','Uniform','Median/power']
    sizes = [r['encoding']['archive_bytes']/1e6 for r in records]
    axes[0].bar(names,sizes,color=['#4b6480','#419b9a','#d66b2a'])
    for x,y in enumerate(sizes): axes[0].text(x,y+10,f'{y:.2f}',ha='center')
    axes[0].set(ylabel='Complete bundle (MB)',title='Identical column/flag/weight policy',ylim=(0,max(sizes)*1.15))
    for name,r,color in zip(names,records,['#4b6480','#419b9a','#d66b2a']):
        axes[1].scatter(r['component_mse'],r['encoding']['archive_bytes']/1e6,s=65,color=color,label=name.replace('\n',' '))
    axes[1].set(xlabel='Unflagged per-component MSE (stored units²)',ylabel='Complete bundle (MB)',title='Visibility rate–distortion')
    axes[1].legend(fontsize=8)
    for extension in ('png','svg'):
        path = ROOT/'docs/paper/figs'/('calibrated-rate.'+extension)
        fig.savefig(path,dpi=160)
        if extension == 'svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)


def image_details(root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2,2,figsize=(10,9),constrained_layout=True)
    for row, kind in enumerate(('uniform','power')):
        for col, product in enumerate(('dirty','image')):
            a, header = read_fits(root/'reference'/('image-'+product+'.fits'))
            b, _ = read_fits(root/(kind+'-images')/('image-'+product+'.fits'))
            valid = sky_mask(a.shape,header)&np.isfinite(a)&np.isfinite(b)
            delta = b-a; limit = float(np.quantile(abs(delta[valid]),.999))
            im = axes[row,col].imshow(np.where(valid,delta,np.nan),origin='lower',cmap='RdBu_r',vmin=-limit,vmax=limit)
            axes[row,col].set_title(kind.capitalize()+' · '+('dirty' if product=='dirty' else 'restored')+' difference')
            axes[row,col].set(xlabel='Pixel x',ylabel='Pixel y')
            fig.colorbar(im,ax=axes[row,col],shrink=.8,label='Jy/beam')
            axes[row,col].text(.02,.02,f'RMS {np.sqrt(np.mean(delta[valid]**2)):.4g} Jy/beam\nmax |Δ| {abs(delta[valid]).max():.4g}',
                               transform=axes[row,col].transAxes,fontsize=9,bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
    fig.suptitle('Calibrated image differences · each scale clips the largest 0.1% of |Δ|',fontsize=12)
    fig.savefig(ROOT/'docs/paper/figs/calibrated-differences.png',dpi=160)
    plt.close(fig)
    for kind in ('lossless','uniform','power'):
        plot_comparison(root/'reference/image-image.fits',root/(kind+'-images')/'image-image.fits',
                        ROOT/'docs/paper/figs'/('calibrated-'+kind+'.png'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT/'results/calibrated-20261023')
    p.add_argument('--archives-only',action='store_true')
    args = p.parse_args()
    records = [summarize(k,args.root) for k in ('lossless','uniform','power')]
    report = {'schema':'varquanms.calibrated-study/1','experiments':records,
              'timing_scope':'Concurrent codec runs on a shared host; 24 quantization workers each. Whole-bundle wall time includes per-partition copying, inspection, encoding, verification and member hashing. No cold-cache or 15-second full-MS claim.',
              'policy':'CORRECTED_DATA only; DATA/FLOAT_DATA/MODEL_DATA removed. Original flags preserved; flagged values discarded and reconstructed as complex positive zero. Four-bit weights are measured exact in this dataset.',
              'limitations':['Stokes I joint imaging only; no per-frame light-curve or Q/U/V validation',
                             'No independent calibration-solution or absolute-flux assessment',
                             'Visibility and weight budgets are engineering settings, not approved science thresholds',
                             'CLEAN uses a 10000-iteration cap; sky-image RMS is not thermal noise',
                             'Four-bit weights are exact here; these images do not test lossy weight approximation',
                             'Whole-source ratios include intentional column/value removal; compare the same-policy lossless control to isolate visibility quantization']}
    if not args.archives_only:
        imported = json.loads((ROOT/'reports/calibrated-import.json').read_text())
        expected = {k:v['sha256'] for k,v in imported['files'].items()}
        if snapshot(ROOT/imported['extracted_ms']) != expected: raise ValueError('Imported source changed')
        report['imported_source_files_unchanged'] = True
        for record in records:
            kind = record['name']
            verification = json.loads((ROOT/'reports'/('calibrated-'+kind+'-export.json')).read_text())
            if verification['status'] != 'passed': raise ValueError('Export verification incomplete')
            record['export_verification'] = verification
            record['imaging'] = json.loads((ROOT/'reports'/('wsclean-calibrated-'+kind+'.json')).read_text())
            run = json.loads((args.root/(kind+'-images')/'complete.json').read_text())
            record['imaging_run'] = {k:run[k] for k in ('status','container_id','parameters','elapsed_seconds','source_unchanged')}
        report['compatible_export'] = json.loads((ROOT/'reports/calibrated-power-compatible-export.json').read_text())
        report['CASA_compatibility'] = json.loads((ROOT/'reports/casa-calibrated-power.json').read_text())
        image_details(args.root)
    output = args.root/'archive-summary.json' if args.archives_only else ROOT/'reports/calibrated-evaluation.json'
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    plot(report)
    for r in records:
        print(r['name'],r['encoding']['archive_bytes'],r['component_mse'],r['component_bit_width_histogram'],r['weights'])


if __name__ == '__main__': main()
