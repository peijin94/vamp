"""Assemble versioned MX evaluation evidence and publication figures."""

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
plt.rcParams['svg.hashsalt']='varquanms-mx-study-v1'


def save(fig, stem):
    target=ROOT/'docs/paper/figs'
    for suffix in ['png','svg','pdf']:
        path=target/f'{stem}.{suffix}'
        metadata={'Date':None} if suffix=='svg' else {'CreationDate':None,'ModDate':None} if suffix=='pdf' else {}
        fig.savefig(path,dpi=180,bbox_inches='tight',metadata=metadata)
        if suffix=='svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)


def diagram():
    fig,ax=plt.subplots(figsize=(11,6))
    ax.set(xlim=(0,11),ylim=(0,6));ax.axis('off')
    def box(x,y,w,h,text,color):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.08',fc='white',ec=color,lw=1.6))
        ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=10)
    box(.3,4.6,10.4,1,'ONE MEGABLOCK: 32 baselines × 24 channels × 1 correlation\n768 complex samples = 1,536 real components = 6,144 original bytes','#2563a6')
    box(.3,3.1,3.0,1,'One outer factor S\nfloat32 per megablock','#16877c')
    box(3.8,3.1,6.9,1,'48 MICROBLOCKS\n24 channels × {real, imaginary}; 32 baselines each','#16877c')
    box(.3,1.4,3.0,1.1,'One UE8M0 scale per microblock\npower of two; 1 byte','#8054a7')
    box(3.8,1.4,3.0,1.1,'Element codes\n4 / 6 / 8 / custom 16 bits\nreal float32 escape available','#8054a7')
    box(7.3,1.4,3.4,1.1,'Measured error checks\nchoose smallest valid record\nexact complex64 fallback','#cc7934')
    for x0,y0,x1,y1 in [(1.8,4.6,1.8,4.15),(7.2,4.6,7.2,4.15),(4.5,3.1,1.8,2.6),(6.5,3.1,5.3,2.6),(8.3,3.1,9,2.6)]:
        ax.annotate('',xy=(x1,y1),xytext=(x0,y0),arrowprops={'arrowstyle':'->','color':'#667785','lw':1.5})
    ax.text(5.5,.65,r'$\hat{x}=\mathrm{float32}\,[S\,2^e\,c]$',ha='center',fontsize=17)
    ax.text(5.5,.15,'Format can be shared across the megablock or chosen per microblock. No channel/time averaging.',ha='center',fontsize=10)
    fig.suptitle('Experimental two-level microscaling: baseline-direction microblocks',fontsize=14)
    save(fig,'mx-grouping')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT/'results/mx-20261023')
    parser.add_argument('--imaging',type=Path,default=ROOT/'reports/wsclean-mx.json')
    args=parser.parse_args()
    samples=[];full=[]
    for path in sorted(args.root.glob('*/report.json')):
        report=json.loads(path.read_text())
        report['experiment']=path.parent.name
        report['configuration_role']='Grouping/error bounds and uniform control; MX formats are defined by variant, not configuration.bits/codebook'
        if path.parent.name=='full-uniform-tight':
            report['timing_interference']='This supplemental fidelity-control run overlapped the WSClean MX imaging job; do not compare its latency with other runs.'
        if report['sample_groups_per_product']:
            samples.append(report)
        else:
            with h5py.File(path.parent/'visibility.mx.h5','r') as h5:
                cross=h5['quality/cross_product_error_covariance'][:]
                adjacent=h5['quality/adjacent_channel_error_covariance'][:]
                channels=h5['quality/channel_component_mse'][:]
                report['covariance_summary']={'cross_product_real':cross.real.tolist(),'cross_product_imag':cross.imag.tolist(),
                    'maximum_absolute_adjacent_channel_covariance_by_product':np.max(np.abs(adjacent),axis=0).tolist(),
                    'maximum_channel_component_mse_by_product':np.max(channels,axis=0).tolist()}
            full.append(report)
    before=json.loads((args.root/'source-before.json').read_text())
    after=json.loads((args.root/'source-after.json').read_text())
    result={'schema':'vamp.mx-study/1','date_utc':'2026-10-23',
            'ocp_specification':'https://www.opencompute.org/wp-content/uploads/2026/10/OCP_MicroXcaling_Spec_v1.1.0-1.pdf',
            'scope':'Experimental DATA component records, not whole standalone MS archives. Standard MX element formats plus custom outer scale, variable formats, binary16 and exact escapes.',
            'sample_runs':samples,'full_runs':full,'source_files_checked':len(before),'source_files_unchanged':before==after,
            'source_file_sha256':before,
            'reference_decoder_sha256':hashlib.sha256((ROOT/'src/vamp/experimental_mx.py').read_bytes()).hexdigest(),
            'sample_selection_caveat':'Stratified deterministic rank sample includes extremes; not a random estimate.',
            'timing_caveat':'Single diagnostic run per configuration, one worker, warm/unspecified shared-host cache. Includes per-block source checks; excludes separate file decoding/export/imaging. Not directly comparable to the production 24-worker 12.19 s encode benchmark.'}
    if args.imaging.exists():result['imaging']=json.loads(args.imaging.read_text())
    image_run=args.root/'wsclean/complete.json'
    if image_run.exists():
        invocation=json.loads(image_run.read_text())
        result['imaging_run']={key:invocation[key] for key in ['status','container_id','wsclean_version','parameters','elapsed_seconds','source_unchanged']}
    (ROOT/'reports/microscaling-evaluation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    diagram()
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
    colors=['#226dac','#b47b28','#23837d','#9b4899']
    widths=['4','6','8','12','16','32']
    for k,r in enumerate(full):
        label={'full-max-micro':'MX + outer factor','full-unit-micro':'MX (S = 1)',
               'full-uniform':'Uniform','full-uniform-tight':'Uniform (tighter bounds)'}[r['experiment']]
        axes[0].scatter(r['serialized_block_bytes']/1e6,r['measured_errors']['component_mse'],s=70,
                        marker=['o','s','^','D'][k%4],color=colors[k%len(colors)],label=label)
        counts={w:0 for w in widths}
        for fmt,n in r['real_component_format_counts'].items():
            w={'exact32':'32','E2M1':'4','E2M3':'6','E3M2':'6','E4M3':'8','E5M2':'8','custom-scaled-FP16':'16','0':'32'}.get(fmt,fmt)
            counts[w]+=n
        width=.8/len(full)
        axes[1].bar(np.arange(len(widths))-.4+width*(k+.5),[100*counts[w]/sum(counts.values()) for w in widths],width,label=label,color=colors[k%len(colors)])
    axes[0].set(xlabel='Serialized visibility blocks (MB; includes headers/scales/tags)',ylabel='Mean component MSE',yscale='log',title='Measured points; same 32 × 24 megablocks')
    axes[0].grid(alpha=.2);axes[0].legend(fontsize=8)
    axes[1].set(xticks=range(len(widths)),xticklabels=widths,xlabel='Bits per real component (32 = exact)',ylabel='Real components (%)',title='Actual selected precision')
    axes[1].legend(fontsize=8)
    fig.suptitle('Full 47,714,304-complex-sample observation; exact flags and auxiliary export\nComponent payload experiment, not whole-archive compression ratios',fontsize=11)
    save(fig,'mx-comparison')
    print(json.dumps({'report':'reports/microscaling-evaluation.json','sample_runs':len(samples),'full_runs':len(full),'source_unchanged':before==after}))


if __name__=='__main__':main()
