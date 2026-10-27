"""Audit two additional epsilon-D cases and build the Dysco paper figures.

``collect`` needs the local MS, archive and WSClean results. ``render`` uses
only the committed JSON reports and the three existing issue figures.
"""
import argparse
from collections import Counter
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import shutil

from report_dysco_comparison import binned_diagnostics, image_metrics
from run_wsclean import snapshot
from stokes_validation import image_paths, reference_path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BASE = ROOT/'reports/dysco-comparison.json'
DEFAULT_EXTRA = ROOT/'reports/dysco-paper-extra.json'
EXTRA_CASES = (('epsd-009-dysco', 0.09), ('epsd-030-dysco', 0.3))
REQUIRED_EXACT = {'UVW', 'FLAG', 'FLAG_ROW', 'WEIGHT', 'SIGMA'}


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):
            digest.update(chunk)
    return digest.hexdigest()


def scientific_config(config):
    """Ignore only execution parallelism and the requested epsilon value."""
    config=deepcopy(config)
    config.pop('workers',None)
    config['normalized_bounds']['epsilon_D']=None
    return config


def collect_extra(root, base_path):
    base=read(base_path)
    if base.get('status')!='complete':
        raise ValueError('A completed Dysco comparison is required')
    reference=next(c for c in base['cases'] if c.get('epsilon_D')==0.18)
    identity_path=ROOT/'results/dysco-comparison/input-identity.json'
    identity=read(identity_path)
    if (identity.get('status')!='passed'
            or sha(identity_path)!=base['provenance_sha256']['results/dysco-comparison/input-identity.json']):
        raise ValueError('The established input identity has changed')
    if sha(ROOT/base['mask']['path'])!=base['mask']['sha256']:
        raise ValueError('The established comparison mask has changed')
    source=Path(identity['original_ms'])
    if snapshot(source)!=identity['original_sha256']:
        raise ValueError('The calibrated reference MS changed')
    results=[]
    for name,epsilon in EXTRA_CASES:
        config_path=root/(name+'-config.json')
        compact_path=root/(name+'-compact.json')
        export_path=root/(name+'-export.json')
        visibility_path=root/(name+'-visibility.json')
        bundle=root/(name+'.bundle')
        manifest_path=bundle/'manifest.json'
        config,compact,export,visibility=map(read,(config_path,compact_path,export_path,visibility_path))
        if (config['normalized_bounds']['epsilon_D']!=epsilon
                or scientific_config(config)!=scientific_config(reference['config'])):
            raise ValueError(f'{name}: comparison policy or epsilon differs')
        if (not compact.get('all_partition_inventories_verified')
                or export.get('status')!='passed' or not export.get('source_files_unchanged')
                or export.get('source_file_sha256')!=identity['original_sha256']
                or 'CORRECTED_DATA' not in export.get('quantized_columns_match_archive_decoder_bitwise',[])):
            raise ValueError(f'{name}: archive or source verification is incomplete')
        exact={c['column'] for c in export['column_checks'] if c['table']=='.' and c['export']=='bitwise_exact'}
        if not REQUIRED_EXACT.issubset(exact):
            raise ValueError(f'{name}: flags, UVW or weights were not checked exactly')
        archive=Path(compact['archive'])
        if not archive.is_file() or archive.stat().st_size!=compact['archive_bytes']:
            raise ValueError(f'{name}: compact archive file length differs')
        if sha(manifest_path)!=(bundle/'manifest.sha256').read_text().strip():
            raise ValueError(f'{name}: bundle manifest checksum differs')
        manifest=read(manifest_path)
        if manifest['config']!=config or manifest['status']!='complete':
            raise ValueError(f'{name}: bundle configuration or completion differs')
        counts=Counter()
        for part in manifest['partitions']:
            for product in part['summary']['visibility_block_statistics']['CORRECTED_DATA']['products']:
                counts.update(product['complex_sample_bit_width_histogram'])
        if sum(counts.values())*8!=base['source_complex64_bytes']:
            raise ValueError(f'{name}: sample precision histogram is incomplete')
        categories=compact['stored_bytes_by_category']
        rate=sum(categories[k] for k in ('visibility_block_payload','permutations','visibility_indexes'))
        restored=Path(export['restored']).resolve()
        if (Path(visibility['source']).resolve()!=source.resolve()
                or Path(visibility['decoded']).resolve()!=restored):
            raise ValueError(f'{name}: visibility diagnostics use different MS files')
        restored_hashes=snapshot(restored)
        images={}
        artifacts={str(path.relative_to(ROOT)):sha(path) for path in
                   (config_path,compact_path,export_path,visibility_path,manifest_path)}
        for pol in ('I','V'):
            runs=[]
            for repeat in (False,True):
                image_dir,comp_path=image_paths(root,name,pol,repeat)
                complete=read(image_dir/'complete.json')
                comp=read(comp_path)
                ref=reference_path(ROOT/'results/calibrated-template',pol)
                if (complete.get('status')!='complete' or complete.get('source_unchanged') is not True
                        or Path(complete['source_ms']).resolve()!=restored
                        or complete['source_sha256']!=restored_hashes
                        or Path(comp['reference']).resolve()!=ref.resolve()
                        or Path(comp['candidate']).resolve()!=image_dir.resolve()
                        or comp['comparison_mask']['sha256']!=base['mask']['sha256']
                        or comp['parameters']['polarizations']!=pol
                        or comp['container_id']!=base['container_id']):
                    raise ValueError(f'{name}: {pol} image or reference identity differs')
                final='Cleaning finished because the final threshold was reached.' in (image_dir/'wsclean.log').read_text()
                runs.append(image_metrics(comp,final))
                artifacts[str(comp_path.relative_to(ROOT))]=sha(comp_path)
                artifacts[str((image_dir/'complete.json').relative_to(ROOT))]=sha(image_dir/'complete.json')
            images['stokes_'+pol]=runs
        results.append({'method':'vamp','name':name,'epsilon_D':epsilon,
                        'config':config,'visibility_payload_bytes':categories['visibility_block_payload'],
                        'visibility_with_permutations_and_indexes_bytes':rate,
                        'complete_output_bytes':compact['archive_bytes'],
                        'stored_bytes_by_category':categories,'precision_sample_counts':dict(counts),
                        'visibility_products':visibility['products'],
                        'visibility_binned_diagnostics':binned_diagnostics(visibility),
                        'export_verified':True,'exact_source_columns':sorted(exact),
                        'artifact_sha256':artifacts,**images})
    if snapshot(source)!=identity['original_sha256']:
        raise ValueError('The calibrated source MS changed during the audit')
    return {'schema':'vamp.dysco-paper-extra/1','status':'complete',
            'base_report':str(base_path.relative_to(ROOT)),'base_report_sha256':sha(base_path),
            'source_complex64_bytes':base['source_complex64_bytes'],
            'comparison_mask':base['mask'],
            'source_unchanged':True,'cases':results}


def collect_series(base, extra):
    """Return sorted actual measured points, separated by method and settings."""
    groups={'vamp':[], 'Dysco RF, 2.5σ':[], 'Dysco RF, 1.5σ':[], 'Dysco AF, 2.5σ':[]}
    for case in base['cases']+extra['cases']:
        if case['method'] in ('vamp', 'varquanMS'):
            key='vamp'
            rate=case['visibility_with_permutations_and_indexes_bytes']
            label=f"εD={case['epsilon_D']:g}"
        elif case.get('flagged_input')=='nan':
            if case['normalization']=='RF' and case['distribution_truncation'] in (1.5,2.5):
                key=f"Dysco RF, {case['distribution_truncation']:g}σ"
            elif case['normalization']=='AF' and case['distribution_truncation']==2.5:
                key='Dysco AF, 2.5σ'
            else:
                continue
            rate=case['visibility_manager_bytes']
            label=f"{case['bits_per_component']} bit"
        else:
            continue
        values=(case['stokes_I'][0]['difference_over_reference_residual_percent'],
                case['stokes_V'][0]['difference_over_reference_residual_percent'],
                100*case['visibility_products']['cross/I-paired']['relative_rmse'])
        if not math.isfinite(rate) or rate<=0 or any(not math.isfinite(x) or x<=0 for x in values):
            raise ValueError('A log-axis comparison requires positive measured values')
        groups[key].append({'x':rate/1e6,'y':values,'label':label,
                            'epsilon_D':case.get('epsilon_D'),
                            'bits_per_component':case.get('bits_per_component'),
                            'clean_I':case['stokes_I'][0]['clean_final_threshold_reached'],
                            'clean_V':case['stokes_V'][0]['clean_final_threshold_reached']})
    for points in groups.values():
        points.sort(key=lambda p:p['x'])
    expected={0.01,0.04,0.09,0.18,0.3}
    if len(groups['vamp'])!=len(expected) or {p['epsilon_D'] for p in groups['vamp']}!=expected:
        raise ValueError('The five requested epsilon values are not all present')
    for key,expected_bits in (('Dysco RF, 2.5σ',{3,4,6,8}),('Dysco RF, 1.5σ',{3,4,6,8}),
                              ('Dysco AF, 2.5σ',{3,8})):
        if len(groups[key])!=len(expected_bits) or {p['bits_per_component'] for p in groups[key]}!=expected_bits:
            raise ValueError('A Dysco comparison series is incomplete')
    return groups


def draw_rate_figure(groups):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    styles={'vamp':('#1d668c','o','-'), 'Dysco RF, 2.5σ':('#c1512c','s','-'),
            'Dysco RF, 1.5σ':('#d59921','D','--'), 'Dysco AF, 2.5σ':('#8b5da7','^',':')}
    fig=plt.figure(figsize=(10,7.2),constrained_layout=True)
    grid=fig.add_gridspec(2,4)
    axes=[fig.add_subplot(grid[0,:2]),fig.add_subplot(grid[0,2:]),
          fig.add_subplot(grid[1,1:3])]
    for series,points in groups.items():
        color,marker,style=styles[series]
        for index,ax in enumerate(axes):
            legend_label=r'vamp ($\epsilon_D$)' if series=='vamp' else series
            ax.plot([p['x'] for p in points],[p['y'][index] for p in points],
                    color=color,marker=marker,linestyle=style,linewidth=1.35,markersize=5,label=legend_label)
            for point in points:
                # The RF series share bit counts and x-coordinates, so annotate
                # one of them; both remain visible as distinct lines/markers.
                if series=='Dysco RF, 1.5σ' and (index==2 or point['clean_I' if index==0 else 'clean_V']):
                    continue
                label=point['label']
                offset=(4,4)
                align='left'
                if series=='vamp':
                    label=f"{point['epsilon_D']:g}"
                    offsets={0.3:(5,7),0.18:(5,6),0.09:(5,6),0.04:(5,6),0.01:(-5,-13)}
                    offset=offsets[point['epsilon_D']]
                    align='right' if point['epsilon_D']==0.01 else 'left'
                    if point['epsilon_D']==0.3 and index==1:
                        offset=(5,-13)
                    if point['epsilon_D']==0.04 and index==1:
                        offset=(5,-13)
                elif point['bits_per_component']==8:
                    offset=(-5,6)
                    align='right'
                if index<2 and not point['clean_I' if index==0 else 'clean_V']:
                    label+='*'
                ax.annotate(label,(point['x'],point['y'][index]),xytext=offset,
                            textcoords='offset points',fontsize=9,color=color,ha=align,
                            bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.5))
    for ax,title in zip(axes,('Stokes I CLEAN','Stokes V CLEAN','Paired cross-baseline I visibility')):
        ax.set_title(title,fontsize=12)
        ax.set_xlabel('Visibility storage (MB)',fontsize=10.5)
        ax.margins(x=.08,y=.12)
        ax.set_yscale('log');ax.invert_yaxis();ax.grid(alpha=.25)
        ax.tick_params(labelsize=9)
    axes[0].set_ylabel('Image difference / residual RMS (%)',fontsize=10.5)
    axes[1].set_ylabel('Image difference / residual RMS (%)',fontsize=10.5)
    axes[2].set_ylabel('Relative visibility RMS error (%)',fontsize=10.5)
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=4,bbox_to_anchor=(.5,-.045),fontsize=9)
    return fig,axes


def render(base_path, extra_path, destination, issue_figure):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    base,extra=read(base_path),read(extra_path)
    if (base.get('status')!='complete' or extra.get('status')!='complete'
            or extra['base_report_sha256']!=sha(base_path)):
        raise ValueError('Comparison reports are incomplete or do not match')
    groups=collect_series(base,extra)
    destination.mkdir(parents=True,exist_ok=True)
    fig,axes=draw_rate_figure(groups)
    if not all(ax.yaxis_inverted() for ax in axes):
        raise ValueError('Error axes must be inverted')
    for suffix in ('png','svg','pdf'):
        path=destination/f'dysco-comparison-rate.{suffix}'
        fig.savefig(path,dpi=190,bbox_inches='tight')
        if suffix=='svg':
            # Matplotlib leaves end-of-line spaces in SVG path data.
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)
    issue_figure.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(destination/'dysco-comparison-rate.png',issue_figure)
    sources={'dysco-comparison-images.png':'dysco-matched-images-rf15.png',
             'dysco-comparison-v-zoom.png':'dysco-matched-v-zoom-rf15.png',
             'dysco-comparison-bits.png':'dysco-matched-bit-distribution.png'}
    for output,source in sources.items():
        src=ROOT/'devnotes/figs'/source
        shutil.copyfile(src,destination/output)
        if sha(src)!=sha(destination/output):
            raise ValueError('A copied issue figure differs')
    return groups


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=('collect','render'))
    parser.add_argument('--root',type=Path,default=ROOT/'results/dysco-paper-extra')
    parser.add_argument('--base-report',type=Path,default=DEFAULT_BASE)
    parser.add_argument('--extra-report',type=Path,default=DEFAULT_EXTRA)
    parser.add_argument('--paper-figures',type=Path,default=ROOT/'docs/paper/figs')
    parser.add_argument('--issue-figure',type=Path,default=ROOT/'devnotes/figs/dysco-paper-rate-distortion.png')
    args=parser.parse_args()
    if args.operation=='collect':
        if args.extra_report.exists():
            raise ValueError('Extra-case report must be new')
        result=collect_extra(args.root.resolve(),args.base_report.resolve())
        args.extra_report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(json.dumps({'status':result['status'],'cases':len(result['cases'])}),flush=True)
    else:
        groups=render(args.base_report,args.extra_report,args.paper_figures,args.issue_figure)
        print(json.dumps({'status':'rendered','series':{k:len(v) for k,v in groups.items()}}),flush=True)


if __name__=='__main__':
    main()
