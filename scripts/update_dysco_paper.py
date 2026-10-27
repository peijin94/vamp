"""Insert the measured four-figure Dysco comparison into section 3.7."""
import argparse
import json
from pathlib import Path

from build_dysco_paper_comparison import collect_series, sha


ROOT=Path(__file__).resolve().parents[1]
START='<!-- DYSCO_COMPARISON_START -->'
END='<!-- DYSCO_COMPARISON_END -->'
HEADING='### 3.7 Comparison'


def row(case):
    if case['method']=='varquanMS':
        label=f"varquanMS $\\epsilon_D={case['epsilon_D']:g}$"
        visibility=case['visibility_with_permutations_and_indexes_bytes']
    else:
        label=f"Dysco {case['normalization']} {case['distribution_truncation']:g}σ, {case['bits_per_component']}-bit"
        visibility=case['visibility_manager_bytes']
    paired=100*case['visibility_products']['cross/I-paired']['relative_rmse']
    image=[]
    for pol in ('I','V'):
        run=case['stokes_'+pol][0]
        mark='' if run['clean_final_threshold_reached'] else '†'
        image.append(f"{run['difference_over_reference_residual_percent']:.3f}{mark}")
    return (f"| {label} | {visibility/1e6:.3f} | {case['complete_output_bytes']/1e6:.3f} | "
            f"{paired:.3f} | {image[0]} | {image[1]} |")


def section(base,extra):
    collect_series(base,extra)  # Enforce the complete five-epsilon and Dysco sets.
    cases=base['cases']+extra['cases']
    varq=sorted((c for c in cases if c['method']=='varquanMS'),key=lambda c:c['epsilon_D'])
    dysco=[c for c in cases if c['method']=='Dysco' and c.get('flagged_input')=='nan']
    dysco.sort(key=lambda c:(0 if c['normalization']=='RF' and c['distribution_truncation']==1.5 else
                             1 if c['normalization']=='RF' else 2,c['bits_per_component']))
    values=varq+dysco
    af8=next(c for c in dysco if c['normalization']=='AF' and c['bits_per_component']==8)
    auto_xx=100*af8['visibility_products']['auto/XX']['relative_rmse']
    auto_yy=100*af8['visibility_products']['auto/YY']['relative_rmse']
    cross_i=100*af8['visibility_products']['cross/I-paired']['relative_rmse']
    a=next(c for c in varq if c['epsilon_D']==.18)
    b=next(c for c in dysco if c['normalization']=='RF' and c['distribution_truncation']==1.5 and c['bits_per_component']==3)
    any_unconverged=any(not c['stokes_'+p][0]['clean_final_threshold_reached'] for c in values for p in ('I','V'))
    lines=[HEADING,'',START,
      'We compare the calibrated archive against one-pass Dysco compression on the',
      'same 15-integration dataset. Both paths retain the original flags, weights,',
      'UVW and other supported metadata, omit raw and model visibilities, and discard',
      'flagged visibility values. The varquanMS export uses zeros for discarded',
      'values; Dysco receives complex NaN markers as described by Offringa (2016).',
      'Unflagged calibrated inputs are identical. The Dysco output is a directly',
      'readable MS, whereas the varquanMS HDF5 archive is exported to an MS for',
      'imaging. Both use the same WSClean I/V settings, reference beams and fixed',
      'above-10° mask. The original source was not modified.',
      '',
      'Dysco uses truncated-Gaussian levels with row-frequency (RF) normalization',
      'at 3, 4, 6 and 8 bits per component and truncations 1.5σ or 2.5σ; the',
      'antenna-frequency (AF) control uses 3 and 8 bits at 2.5σ. VarquanMS',
      'uses the same 128-baseline × 48-channel grouping and adaptive ladder as',
      'above, with additional $\\epsilon_D=0.09$ and $0.3$ measurements. Table 1',
      'reports actual bytes. Visibility storage includes varquanMS permutations',
      'and indexes, or Dysco DATA-manager scales and headers; complete sizes also',
      'include the different container and metadata representations.',
      '',
      'The paired-I column is the relative visibility RMS error for',
      '$(XX+YY)/2$ over unflagged cross-baselines. The I/V columns use',
      '$100\\,\\mathrm{RMS}(I_{\\mathrm{quant}}-I_{\\mathrm{ref}})/\\mathrm{RMS}(I_{\\mathrm{residual}})$',
      'and the corresponding V expression:',
      'the denominators are reference residual RMS, not image flux or thermal',
      'noise. Values are first CLEAN runs.',
      '',
      '| Method and setting | Vis. MB | Total MB | Paired-I % | I % | V % |',
      '| -------- | ---: | ---: | ---: | ---: | ---: |']
    lines.extend(row(c) for c in values)
    lines.extend(['',
      '*Table 1. Measured calibrated-data sizes and errors. MB uses $10^6$ bytes;',
      'visibility and whole-output sizes have different accounting boundaries.',
      ('A dagger denotes a CLEAN that did not reach its final threshold.' if any_unconverged
       else 'All listed CLEANs reached their configured final threshold.'),
      'See the machine-readable reports for repeats, absolute errors, per-product',
      'and baseline/frequency diagnostics, exact-column checks and provenance.*',
      '',
      '![Visibility-storage rate and distortion](figs/dysco-comparison-rate.png)',
      '',
      '*Figure 10. I/V restored-image differences and paired cross-baseline I',
      'visibility error versus visibility storage. Points at $\\epsilon_D=0.09$',
      'and $0.3$ are newly measured; labels on the varquanMS points are $\\epsilon_D$.',
      'Lines join measured settings of the same',
      'method; they do not establish interpolation or monotonic CLEAN behavior.',
      'Error axes are reversed, placing smaller errors higher. RF truncations use',
      ('separate series; AF has only the two measured bit widths.'
       + (' An asterisk marks a CLEAN that did not reach its final threshold.' if any_unconverged else '')+'*'),
      '',
      'At nearly identical visibility storage, the $\\epsilon_D=0.18$ archive uses',
      f"{a['visibility_with_permutations_and_indexes_bytes']/1e6:.3f} MB, compared with {b['visibility_manager_bytes']/1e6:.3f} MB",
      'for 3-bit RF at 1.5σ. Their first I/V differences are',
      f"{a['stokes_I'][0]['difference_over_reference_residual_percent']:.2f}% / {a['stokes_V'][0]['difference_over_reference_residual_percent']:.2f}% and",
      f"{b['stokes_I'][0]['difference_over_reference_residual_percent']:.2f}% / {b['stokes_V'][0]['difference_over_reference_residual_percent']:.2f}% of the matching reference residual RMS.",
      'The comparison is specific to this observation; output images are not',
      'independent evidence of the true sky or a general fidelity guarantee.',
      '',
      '![Matched I and V image differences](figs/dysco-comparison-images.png)',
      '',
      '*Figure 11. Full-field reference I and V images and signed differences',
      'for $\\epsilon_D=0.18$ and Dysco RF 3-bit, 1.5σ. Each difference is',
      'expressed as a percentage of the matching full-mask reference residual',
      'RMS. Reference I uses viridis, reference V and signed differences use bwr;',
      'the difference scale is limited to ±30%, with extensions for clipped values.*',
      '',
      '![Taurus A Stokes V comparison](figs/dysco-comparison-v-zoom.png)',
      '',
      '*Figure 12. Stokes V reference and reconstructions around the selected',
      'Taurus A feature, with signed differences below. The 64×64 crop is above',
      '45° throughout the observation; differences use the same full above-10°',
      'reference residual RMS. Intrinsic circular polarization and instrumental',
      'leakage have not been separated.*',
      '',
      '![Precision distribution at the matched rate](figs/dysco-comparison-bits.png)',
      '',
      '*Figure 13. Complex-sample counts by selected precision for varquanMS',
      '$\\epsilon_D=0.18$ and 3-bit Dysco RF. Counts include flagged sample slots;',
      'varquanMS exact-zero blocks are stored compactly. These are sample shares,',
      'not byte shares; Table 1 includes the corresponding storage costs.*',
      '',
      'CLEAN can be nonmonotonic: for Dysco RF at 1.5σ, the 4-bit I dirty-image',
      'difference is smaller than at 3 bits, while its restored-image difference',
      'is larger. Repeated images and stochastic-encoding variation have distinct',
      'causes; the measured points should not be interpreted as a monotonic',
      'guarantee. Autocorrelations are excluded from these images. In the full',
      f"8-bit AF run, XX and YY autocorrelation relative RMS errors are {auto_xx:.2f}% and {auto_yy:.2f}%,",
      f'while the cross-baseline paired-I error is {cross_i:.2f}%. Image agreement alone',
      'therefore does not establish autocorrelation fidelity.',
      '',
      '[Dysco comparison data](../../reports/dysco-comparison.json);',
      '[additional epsilon measurements](../../reports/dysco-paper-extra.json);',
      '[detailed comparison notes](../../devnotes/compare-dysco.md).',
      END])
    return '\n'.join(lines)


def update(text,body):
    old='new, and a matched Dysco comparison has not yet been performed.'
    new='new. A matched Dysco comparison is reported in Section 3.7.'
    if old in text:
        text=text.replace(old,new,1)
    elif new not in text:
        raise ValueError('Expected introduction comparison statement was not found')
    old_conclusion=('not a CASA 6 or precision-astrometry validation. A matched Dysco benchmark,\n'
                    'additional observations and dedicated weak-polarization checks remain necessary')
    new_conclusion=('not a CASA 6 or precision-astrometry validation. The matched Dysco comparison\n'
                    'covers the same observation; additional observations and dedicated weak-polarization checks remain necessary')
    if old_conclusion in text:
        text=text.replace(old_conclusion,new_conclusion,1)
    elif new_conclusion not in text:
        raise ValueError('Expected conclusion comparison statement was not found')
    if START in text or END in text:
        if text.count(START)!=1 or text.count(END)!=1 or text.count(HEADING)!=1:
            raise ValueError('Comparison markers are malformed')
        start=text.index(HEADING);end=text.index(END,start)+len(END)
        text=text[:start]+body+text[end:]
    else:
        marker='\n## 4. Conclusion and discussion'
        if text.count(marker)!=1 or HEADING in text:
            raise ValueError('Expected conclusion boundary was not found exactly once')
        text=text.replace(marker,'\n'+body+'\n'+marker,1)
    return text


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-report',type=Path,default=ROOT/'reports/dysco-comparison.json')
    parser.add_argument('--extra-report',type=Path,default=ROOT/'reports/dysco-paper-extra.json')
    parser.add_argument('--paper',type=Path,default=ROOT/'docs/paper/varquantMS.md')
    args=parser.parse_args()
    base=json.loads(args.base_report.read_text());extra=json.loads(args.extra_report.read_text())
    if base.get('status')!='complete' or extra.get('status')!='complete' or extra['base_report_sha256']!=sha(args.base_report):
        raise ValueError('Comparison reports are incomplete or mismatched')
    updated=update(args.paper.read_text(),section(base,extra))
    args.paper.write_text(updated)
    print(json.dumps({'paper':str(args.paper),'section':'3.7 Comparison','status':'updated'}),flush=True)


if __name__=='__main__':
    main()
