"""Deterministic vector method figures and measured epsilon-D validation plots."""
from pathlib import Path
import json
import os
import shutil

os.environ.setdefault('MPLCONFIGDIR','/tmp/varquanms-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

DEST=Path('docs/paper/figs')


def save(fig,name):
    for suffix in ('png','svg','pdf'):
        path=DEST/(name+'.'+suffix)
        fig.savefig(path,dpi=180,bbox_inches='tight')
        if suffix=='svg':
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)


def methods():
    DEST.mkdir(exist_ok=True)
    fig,ax=plt.subplots(figsize=(12,6.7));ax.set(xlim=(0,12),ylim=(0,7));ax.axis('off')
    def box(x,y,w,h,text,color='#e5eff7'):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.10',facecolor=color,edgecolor='#30566c',linewidth=1.2))
        ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=10)
    def arrow(a,b):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=14,color='#30566c',linewidth=1.3))
    box(.2,5.55,3.1,1.1,'Inspected calibrated MS\nPartition time / configuration\nSeparate product and auto/cross')
    box(4.05,5.55,3.1,1.1,'Stable amplitude order\nKeep each baseline spectrum\nStore inverse reconstruction map')
    box(7.9,5.55,3.8,1.1,'Block: 128 baselines × 48 channels\nOne correlation product\nExact edge bounds')
    arrow((3.4,6.1),(3.95,6.1));arrow((7.25,6.1),(7.8,6.1))
    box(.2,3.0,3.1,1.25,'Original valid samples\nComponent median deviations\nD = pooled raw MAD', '#e6f2e6')
    box(4.05,3.0,3.1,1.25,'Local factors\nComplex median C; range t\nGamma from max/std', '#e6f2e6')
    box(7.9,2.8,3.8,1.65,'Candidate encode → decode\n3, 4, 5, 6, 8, 10, 12, 14, 16 bits\nReal + imaginary, same width\nExact float32 always competes','#fff0d9')
    arrow((9.8,5.45),(9.8,4.55));arrow((8.0,5.45),(5.6,4.35));arrow((8.0,5.45),(1.75,4.35))
    box(.2,.6,3.1,1.35,'Measured error guards\nT = epsilon_D × D\nMSE, maximum, bias, baseline', '#e6f2e6')
    box(4.05,.6,3.1,1.35,'Smallest accepted bytes\nCount factors and exceptions\nFour-bit precision selector','#fff0d9')
    box(7.9,.6,3.8,1.35,'Independently readable block\nPacked codes + exact exceptions\nChecksums, offsets, metadata')
    arrow((1.75,2.9),(1.75,2.05));arrow((7.25,3.6),(7.8,3.6))
    arrow((9.7,2.7),(1.75,2.05));arrow((3.4,1.3),(3.95,1.3));arrow((7.25,1.3),(7.8,1.3))
    fig.suptitle('Grouping, local factors and adaptive precision',fontsize=15)
    save(fig,'epsd-method')
    fig,axes=plt.subplots(1,2,figsize=(11,3.8),constrained_layout=True)
    codes=np.arange(16);u=2*codes/15-1
    for i,(gamma,color) in enumerate(((1.,'#2166ac'),(.5,'#d6604d'),(.25,'#238b45'))):
        levels=np.sign(u)*abs(u)**(1/gamma)
        axes[0].scatter(levels,np.full(16,2-i),color=color,s=30)
        for x,k in zip(levels,codes):
            if k in ((0,7,8,15) if gamma==1 else (0,15)):
                axes[0].annotate(str(k),(x,2-i),xytext=(0,8),textcoords='offset points',ha='center',fontsize=7)
        dense=np.linspace(-1,1,401)
        axes[1].plot(dense,np.sign(dense)*abs(dense)**(1/gamma),color=color,label=f'γ = {gamma:g}')
    axes[0].axvline(0,color='grey',linewidth=.6,linestyle=':')
    axes[0].set(yticks=[2,1,0],yticklabels=['γ = 1','γ = ½','γ = ¼'],xlabel='Reconstructed residual / stored range',ylim=(-.4,2.5),title='All 16 unsigned codes; b = 4')
    axes[1].set(xlabel='Uniform code coordinate u = 2k / (2ᵇ−1) − 1',ylabel='Reconstructed residual / stored range',title='Inverse power mapping')
    axes[1].legend();axes[1].grid(alpha=.2)
    save(fig,'epsd-levels')


def measured(result):
    from report_fixed_rate import plots
    root=Path(result['results_root'])
    plots(result,root=root,prefix='epsd-dense')
    for suffix in ('rate','bits','images','zoom','best-images','best-zoom'):
        shutil.copyfile(Path('devnotes/figs')/f'epsd-dense-{suffix}.png',DEST/f'epsd-{suffix}.png')
    cases=sorted(result['experiments'],key=lambda e:e['epsilon_D'])
    fig,axes=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    for e in cases:
        x=e['shared']['archive_bytes']/1e6
        yy=[100*r['restored_difference_over_reference_residual'] for r in e['residual_normalized_runs']]
        axes[0].plot([x,x],yy,'o-',color='#2166ac')
        axes[0].annotate(f"εD={e['epsilon_D']:g}",(x,max(yy)),xytext=(0,7),textcoords='offset points',ha='center',fontsize=9)
    axes[0].set(xlabel='Complete archive (MB)',
                ylabel='RMS(image − reference) /\nRMS(reference residual) [%]',
                title='Two CLEAN runs per fixed epsilon')
    axes[0].margins(x=.16,y=.18)
    main=('visibility_block_payload','visibility_quality','quantization_statistics','coordinates')
    left=np.zeros(len(cases))
    for k,color in zip((*main,'other'),('#2166ac','#67a9cf','#ef8a62','#fddbc7','#aaa')):
        values=np.array([sum(v for key,v in e['shared']['stored_bytes_by_category'].items() if (key==k if k!='other' else key not in main))/1e6 for e in cases])
        axes[1].barh(range(len(cases)),values,left=left,color=color,label=k.replace('_',' '));left+=values
    axes[1].set(yticks=range(len(cases)),yticklabels=[f"εD={e['epsilon_D']:g}" for e in cases],xlabel='Physical archive bytes (MB)',title='Whole-archive byte categories')
    axes[1].legend(loc='upper center',bbox_to_anchor=(.5,-.2),fontsize=7,ncol=2)
    for ax in axes:ax.grid(axis='x',alpha=.2)
    save(fig,'epsd-rate-distortion')


if __name__=='__main__':
    methods()
    path=Path('reports/calibrated-epsd-dense.json')
    if path.exists():measured(json.loads(path.read_text()))
