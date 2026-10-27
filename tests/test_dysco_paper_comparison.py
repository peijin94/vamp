"""The paper plot must use the requested cases and honest axis semantics."""
import copy
import importlib.util
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import build_dysco_paper_comparison as paper
import update_dysco_paper as manuscript


def case(epsilon=None,bits=None,normalization=None,truncation=None,flags='nan',rate=100):
    result={'method':'vamp' if epsilon is not None else 'Dysco',
            'stokes_I':[{'difference_over_reference_residual_percent':1.,'clean_final_threshold_reached':True}],
            'stokes_V':[{'difference_over_reference_residual_percent':2.,'clean_final_threshold_reached':True}],
            'visibility_products':{'cross/I-paired':{'relative_rmse':.03},
                                   'auto/XX':{'relative_rmse':.97},'auto/YY':{'relative_rmse':.96}},
            'complete_output_bytes':rate+10}
    if epsilon is not None:
        result.update(epsilon_D=epsilon,visibility_with_permutations_and_indexes_bytes=rate)
    else:
        result.update(bits_per_component=bits,normalization=normalization,distribution_truncation=truncation,
                      flagged_input=flags,visibility_manager_bytes=rate)
    return result


class DyscoPaperTests(unittest.TestCase):
    def sample(self):
        original=[case(epsilon=e,rate=r) for e,r in ((.01,400),(.04,300),(.18,200))]
        for trunc in (1.5,2.5):
            original.extend(case(bits=b,normalization='RF',truncation=trunc,rate=100+b)
                            for b in (3,4,6,8))
        original.extend(case(bits=b,normalization='AF',truncation=2.5,rate=90+b) for b in (3,8))
        original.append(case(bits=3,normalization='RF',truncation=2.5,flags='zero'))
        extra=[case(epsilon=.09,rate=250),case(epsilon=.3,rate=150)]
        return {'cases':original},{'cases':extra}

    def test_additions_preserve_scientific_configuration(self):
        original={'workers':24,'baseline_group':128,'normalized_bounds':{'epsilon_D':.18,'max_error_factor':10}}
        changed=copy.deepcopy(original);changed['workers']=6;changed['normalized_bounds']['epsilon_D']=.09
        self.assertEqual(paper.scientific_config(original),paper.scientific_config(changed))
        changed['baseline_group']=32
        self.assertNotEqual(paper.scientific_config(original),paper.scientific_config(changed))
        self.assertEqual(original['normalized_bounds']['epsilon_D'],.18)

    def test_series_use_both_new_epsilons_and_exclude_zero_control(self):
        base,extra=self.sample(); groups=paper.collect_series(base,extra)
        self.assertEqual([p['epsilon_D'] for p in groups['vamp']],[.3,.18,.09,.04,.01])
        self.assertEqual(len(groups['Dysco RF, 2.5σ']),4)
        self.assertEqual(len(groups['Dysco RF, 1.5σ']),4)
        self.assertEqual(len(groups['Dysco AF, 2.5σ']),2)
        extra['cases'][0]['epsilon_D']=.9
        with self.assertRaisesRegex(ValueError,'five requested'):
            paper.collect_series(base,extra)

    def test_legacy_measurement_labels_render_as_vamp(self):
        base,extra=self.sample()
        for case in base['cases']+extra['cases']:
            if case['method']=='vamp':
                case['method']='varquanMS'
        groups=paper.collect_series(base,extra)
        self.assertEqual(len(groups['vamp']),5)
        self.assertEqual(manuscript.section(base,extra).count('\n| vamp'),5)

    def test_paper_inserts_four_figures_and_both_corrected_epsilons_once(self):
        base,extra=self.sample(); body=manuscript.section(base,extra)
        self.assertIn(r'$\epsilon_D=0.09$',body)
        self.assertIn(r'$\epsilon_D=0.3$',body)
        self.assertEqual(body.count('![',),4)
        self.assertEqual(body.count('\n| vamp'),5)
        self.assertEqual(body.count('\n| Dysco'),10)
        intro='new, and a matched Dysco comparison has not yet been performed.'
        conclusion=('not a CASA 6 or precision-astrometry validation. A matched Dysco benchmark,\n'
                    'additional observations and dedicated weak-polarization checks remain necessary')
        text=intro+'\n\n## 4. Conclusion and discussion\n\n'+conclusion
        updated=manuscript.update(text,body)
        self.assertLess(updated.index('### 3.7 Comparison'),updated.index('## 4. Conclusion'))
        self.assertEqual(manuscript.update(updated,body),updated)

    def test_nonconverged_clean_is_explicitly_marked(self):
        base,extra=self.sample()
        extra['cases'][1]['stokes_I'][0]['clean_final_threshold_reached']=False
        body=manuscript.section(base,extra)
        row=next(line for line in body.splitlines() if r'vamp $\epsilon_D=0.3$' in line)
        self.assertIn('1.000†',row)
        self.assertIn('A dagger denotes a CLEAN that did not reach its final threshold.',body)

    @unittest.skipUnless(importlib.util.find_spec('matplotlib'),'Install plotting support')
    def test_each_error_axis_is_reversed_and_series_are_connected(self):
        import matplotlib.pyplot as plt
        base,extra=self.sample()
        extra['cases'][1]['stokes_I'][0]['clean_final_threshold_reached']=False
        rf=next(c for c in base['cases'] if c.get('normalization')=='RF'
                and c.get('distribution_truncation')==1.5 and c.get('bits_per_component')==3)
        rf['stokes_V'][0]['clean_final_threshold_reached']=False
        groups=paper.collect_series(base,extra)
        fig,axes=paper.draw_rate_figure(groups)
        try:
            self.assertEqual(len(axes),3)
            self.assertIn('0.3*',[t.get_text() for t in axes[0].texts])
            self.assertIn('3 bit*',[t.get_text() for t in axes[1].texts])
            for ax in axes:
                self.assertTrue(ax.yaxis_inverted())
                self.assertEqual(ax.get_yscale(),'log')
                self.assertEqual(len(ax.lines),4)
                for line in ax.lines:
                    self.assertNotEqual(line.get_linestyle(),'None')
                    self.assertEqual(list(line.get_xdata()),sorted(line.get_xdata()))
        finally:
            plt.close(fig)


if __name__=='__main__':
    unittest.main()
