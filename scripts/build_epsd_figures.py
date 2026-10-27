"""Deterministic vector method figures and measured epsilon-D validation plots."""
from pathlib import Path
import json
import datetime
import os

os.environ.setdefault('MPLCONFIGDIR','/tmp/varquanms-mpl')
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt'] = 'varquanms-paper'
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
import numpy as np

DEST=Path('docs/paper/figs')


def save(fig,name):
    for suffix in ('png','svg','pdf'):
        path=DEST/(name+'.'+suffix)
        stamp = datetime.datetime(2026, 10, 26, tzinfo=datetime.timezone.utc)
        metadata = ({'Date': stamp.isoformat()} if suffix == 'svg' else
                    {'CreationDate': stamp, 'ModDate': stamp} if suffix == 'pdf' else None)
        fig.savefig(path,dpi=180,bbox_inches='tight',metadata=metadata)
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


def image_comparisons(result):
    """Recovered I or V over RMS-normalized percent error, three epsilon columns."""
    from astropy.io import fits
    from matplotlib.colors import Normalize, SymLogNorm
    from matplotlib.ticker import PercentFormatter
    from varquanms.imaging import read_fits, _grid, file_sha256
    from report_fixed_rate import REF

    pol = result.get('stokes', 'I')
    if pol not in ('I', 'V'):
        raise ValueError('Expected a single Stokes I or V study')
    cases = sorted(result['experiments'], key=lambda e: e['epsilon_D'])
    if len(cases) != 3:
        raise ValueError('The paper comparison requires three epsilon cases')
    mask_info = result['comparison_mask']
    mask_path = Path(mask_info['path'])
    if file_sha256(mask_path) != mask_info['sha256']:
        raise ValueError('Comparison mask differs from the measured study')
    mask = fits.getdata(mask_path).astype(bool)
    reference_directory = Path(result.get('reference_directory', REF))
    reference_path = reference_directory / 'image-image.fits'
    reference, header = read_fits(reference_path)
    residual, residual_header = read_fits(reference_directory / 'image-residual.fits')
    if (mask.shape != reference.shape or residual.shape != reference.shape
            or _grid(header) != _grid(residual_header)
            or int(mask.sum()) != mask_info['pixels']):
        raise ValueError('Reference, residual and mask must match the measured domain')
    if not np.all(np.isfinite(residual[mask])):
        raise ValueError('Reference residual has nonfinite comparison pixels')
    # One scalar, uncentered residual RMS over the FULL mask;
    # never re-estimate it in a zoom or from the quantized residual.
    residual_rms = float(np.sqrt(np.mean(residual[mask]**2)))
    expected = result['reference_residual_statistics']['rms']
    if not np.isfinite(residual_rms) or residual_rms <= 0 or not np.isclose(residual_rms, expected, rtol=1e-12):
        raise ValueError('Reference residual RMS no longer matches the report')
    reference_sha = file_sha256(reference_path)
    images, differences = [], []
    for case in cases:
        measured = case['comparison']['products']['image-image.fits']
        path = Path(measured.get('candidate', Path(result['results_root']) / (case['name'] + '-images') / 'image-image.fits'))
        if (reference_sha != measured['reference_file_sha256']
                or file_sha256(path) != measured['candidate_file_sha256']):
            raise ValueError('Image inputs differ from the recorded first CLEAN run')
        recovered, recovered_header = read_fits(path)
        if (recovered.shape != reference.shape or _grid(recovered_header) != _grid(header)
                or not np.all(np.isfinite(recovered[mask]))):
            raise ValueError('Recovered image grid/beam or comparison domain differs')
        error = recovered - reference
        rms = float(np.sqrt(np.mean(error[mask] ** 2)))
        if not np.isclose(rms, measured['regions']['above_elevation_cut']['difference_rms'], rtol=1e-12):
            raise ValueError('Image difference disagrees with measured RMS')
        images.append(recovered)
        differences.append(100.0 * error / residual_rms)

    peak = np.unravel_index(np.argmax(np.where(mask, np.abs(reference), -np.inf)), reference.shape)
    zoom_view = (slice(max(0, peak[0]-64), min(reference.shape[0], peak[0]+64)),
                 slice(max(0, peak[1]-64), min(reference.shape[1], peak[1]+64)))
    if pol == 'V':
        zoom_record = json.loads(Path('reports/calibrated-v-zoom.json').read_text())
        if (zoom_record['reference_fits_sha256'] != reference_sha
                or zoom_record['comparison_mask_sha256'] != mask_info['sha256']
                or zoom_record['entire_crop_minimum_elevation_deg'] <= 45.):
            raise ValueError('V zoom does not match the measured reference or high-elevation requirement')
        bounds = zoom_record['bounds']
        zoom_view = (slice(bounds['y_start'], bounds['y_stop']), slice(bounds['x_start'], bounds['x_stop']))
        if reference[zoom_view].shape != tuple(zoom_record['crop_shape']) or not mask[zoom_view].all():
            raise ValueError('V zoom must match its complete, valid recorded crop')
    scale = 30. if pol == 'V' else (max(float(np.max(np.abs(a[mask]))) for a in [reference, *images]) or 1.)
    limit = 30.0  # Fixed display range in percent; underlying errors are unchanged.
    intensity_norm = (Normalize(vmin=-scale, vmax=scale) if pol == 'V'
                      else SymLogNorm(.01 * scale, vmin=-scale, vmax=scale))
    DEST.mkdir(exist_ok=True)
    for zoom in (False, True):
        view = zoom_view if zoom else np.s_[:,:]
        extent = ((view[1].start-.5, view[1].stop-.5, view[0].start-.5, view[0].stop-.5)
                  if zoom else (-.5, reference.shape[1]-.5, -.5, reference.shape[0]-.5))
        fig, axes = plt.subplots(2, 3, figsize=(12.4, 8.5), layout='constrained')
        for col, (case, recovered, difference) in enumerate(zip(cases, images, differences)):
            top = axes[0, col].imshow(np.where(mask, recovered, np.nan)[view], origin='lower',
                                      extent=extent, cmap='viridis' if pol=='I' else 'bwr', norm=intensity_norm)
            if zoom and pol == 'V':
                source_y, source_x = zoom_record['center_yx']
                axes[0, col].add_patch(Circle((source_x, source_y), 5, fill=False,
                                              edgecolor='black', linewidth=.8, linestyle='--'))
            bottom = axes[1, col].imshow(np.where(mask, difference, np.nan)[view], origin='lower',
                                         extent=extent, cmap='bwr', vmin=-limit, vmax=limit)
            axes[0, col].set_title(f"εD = {case['epsilon_D']:g} · Stokes {pol}", fontsize=11)
            axes[1, col].set_title('Difference', fontsize=11)
        intensity_bar = fig.colorbar(top, ax=axes[0, :].tolist(), shrink=.85, pad=.02,
                                     label=rf'${pol}_{{\rm quant}}$ [Jy/beam]')
        ticks = np.linspace(-scale, scale, 5) if pol == 'V' else [-scale, -.1*scale, 0, .1*scale, scale]
        intensity_bar.set_ticks(ticks)
        intensity_bar.set_ticklabels([f'{value:g}' if pol == 'V' else f'{value:.2g}' for value in ticks])
        difference_bar = fig.colorbar(bottom, ax=axes[1, :].tolist(), shrink=.85, pad=.02,
                                      extend='both', ticks=np.arange(-30, 31, 10),
                                      label=rf'$({pol}_{{\rm quant}}-{pol}_{{\rm ref}})/\mathrm{{RMS}}_{{\rm residual}}$')
        difference_bar.ax.yaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
        for row in range(2):
            for col in range(3):
                ax = axes[row, col]
                if row == 0:
                    ax.xaxis.set_visible(False)
                else:
                    ax.set_xlabel('Pixel x')
                if col > 0:
                    ax.yaxis.set_visible(False)
                else:
                    ax.set_ylabel('Pixel y')
        view_label = 'reference-peak crop' if zoom else 'full field'
        elevation_label = 'above 10°'
        if zoom and pol == 'V':
            view_label, elevation_label = zoom_record['source']['label'], 'above 45°'
        fig.suptitle(f'Stokes {pol} · {view_label} · {elevation_label} · reference residual RMS = {residual_rms:.6f} Jy/beam',
                     fontsize=12)
        prefix = 'epsd' if pol=='I' else 'epsd-v'
        fig.savefig(DEST / (prefix+('-zoom.png' if zoom else '-images.png')), dpi=180)
        plt.close(fig)


def measured(result):
    from report_fixed_rate import plots
    root=Path(result['results_root'])
    plots(result,root=root,prefix='epsd',output=DEST,include_best=False,include_images=False)
    image_comparisons(result)
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
    path=Path('reports/calibrated-epsd-v.json')
    if path.exists():image_comparisons(json.loads(path.read_text()))
