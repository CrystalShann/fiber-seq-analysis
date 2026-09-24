"""Validate pooled focal-allele summaries and REF/ALT-only plots."""
import importlib
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
s = importlib.import_module('05_summarize_autocorrelations')
REGION = SimpleNamespace(region_id='test', focal_snp='rs1', ref='G', alt='A',
                         chr='chr1', analysis_start=100, analysis_end=399)


def fixture(groups, width=300):
    rows, profiles = [], []
    for sample, allele, n, value in groups:
        for _ in range(n):
            rid = f'r{len(rows)}'
            rows.append(dict(row_id=sample+'__'+rid, RID=rid, sample_name=sample,
                haplotype='HP1' if allele=='REF' else 'HP2', focal_genotype='0|1',
                allele_status='phased_focal_genotype', allele_label='rs1: '+('G' if allele=='REF' else 'A'),
                acf_valid=True, status='clustered', cluster='0'))
            p = np.zeros(width)
            p[0] = 1
            if width > 180:
                p[180] = value
            profiles.append(p)
    metadata = pd.DataFrame(dict(sample_name=sorted(set(r['sample_name'] for r in rows))))
    metadata['cell_line'] = metadata.sample_name.str.split('_').str[0]
    records = s.annotate_alleles(pd.DataFrame(rows), REGION, metadata)
    return np.array(profiles), records


class AlleleSummaries(unittest.TestCase):
    def test_hp_labels_can_reverse_without_reversing_ref_alt(self):
        _,r = fixture([('NA1_a','REF',1,.1),('NA1_a','ALT',1,.2)])
        flipped = r.copy()
        flipped['haplotype'] = ['HP2','HP1']
        flipped['focal_genotype'] = '1|0'
        table=pd.DataFrame(dict(sample_name=['NA1_a'],cell_line=['NA1']))
        self.assertEqual(s.annotate_alleles(flipped,REGION,table).allele_group.tolist(), ['REF','ALT'])
        flipped.loc[0,'allele_label']='rs1: A'
        with self.assertRaisesRegex(ValueError,'disagrees'):
            s.annotate_alleles(flipped,REGION,table)

    def test_pooled_molecule_weighting_and_no_sample_coverage_filter(self):
        p,r=fixture([('NA1_a','REF',1,.1),('NA1_a','ALT',1,.3),
                     ('NA2_a','REF',9,.1),('NA2_a','ALT',9,.7)])
        t=s.summarize_alleles(p,r,REGION)
        d=t['allele_acf_summary'].query("allele_group == 'ALT' and lag_bp == 180").iloc[0]
        self.assertAlmostEqual(d.mean_acf,.66)
        self.assertEqual(d.n_reads,10)
        renamed=r.assign(cell_line='one',sample_name='anonymous')
        pd.testing.assert_frame_equal(t['allele_acf_summary'],
            s.summarize_alleles(p,renamed,REGION)['allele_acf_summary'])
        for name,table in t.items():
            if name != 'allele_read_audit':
                self.assertFalse({'cell_line','sample_name','sample_names'} & set(table.columns))

    def test_zero_cluster_count_differs_from_missing_allele(self):
        p,r=fixture([('NA1_a','REF',2,.1),('NA1_a','ALT',2,.2)])
        r.loc[r.allele_group.eq('ALT'),'cluster']='1'
        c=s.summarize_alleles(p,r,REGION)['allele_cluster_proportions']
        self.assertEqual(c.query("allele_group == 'REF' and cluster == '1'").iloc[0].fraction,0)
        c=s.summarize_alleles(p[:2],r.iloc[:2],REGION)['allele_cluster_proportions']
        self.assertTrue(c.query("allele_group == 'ALT'").fraction.isna().all())

    def test_valid_unclustered_reads_contribute_to_acf_not_cluster_denominator(self):
        p,r=fixture([('NA1_a','REF',2,.1),('NA1_a','ALT',2,.2)])
        r.loc[0,'status']='insufficient_reads'
        t=s.summarize_alleles(p,r,REGION)
        c=t['allele_coverage'].set_index('allele_group')
        self.assertEqual(c.loc['REF','n_valid'],2)
        self.assertEqual(c.loc['REF','n_clustered'],1)
        self.assertEqual(c.loc['REF','n_valid_not_clustered'],1)

    def test_no_peak_is_zero_fraction_but_undefined_lag(self):
        p,r=fixture([('NA1_a','REF',2,0),('NA1_a','ALT',2,.2)])
        row=s.summarize_alleles(p,r,REGION)['allele_peak_features'].query("allele_group == 'REF'").iloc[0]
        self.assertEqual(row.peak_fraction,0)
        self.assertTrue(np.isnan(row.median_peak_lag))

    def test_truncated_band_is_not_reported_as_no_peak(self):
        p,r=fixture([('NA1_a','REF',2,.1),('NA1_a','ALT',2,.2)],width=251)
        t=s.summarize_alleles(p,r,REGION)
        self.assertTrue(t['allele_peak_features'].peak_fraction.isna().all())
        self.assertTrue(np.isnan(s.repeat_peak(p[0])[0]))

    def test_duplicate_original_and_merged_molecule_count_once(self):
        p,r=fixture([('NA1_a','REF',2,.1),('NA1_b','ALT',2,.2)])
        copy=r.iloc[[0]].copy()
        copy['sample_name']='NA1_b'
        copy['row_id']='NA1_b__'+copy.RID
        r=pd.concat([r,copy],ignore_index=True)
        p=np.vstack([p,p[0]])
        t=s.summarize_alleles(p,r,REGION)
        self.assertEqual(t['allele_coverage'].n_unique.sum(),4)
        self.assertEqual(t['allele_coverage'].n_duplicates_removed.sum(),1)
        r.loc[len(r)-1,'cluster']='1'
        with self.assertRaisesRegex(ValueError,'Conflicting duplicate'):
            s.summarize_alleles(p,r,REGION)

    def test_invalid_acf_is_excluded_from_peak_denominator(self):
        p,r=fixture([('NA1_a','REF',3,.1),('NA1_a','ALT',2,.2)])
        p[0]=np.nan
        r.loc[0,['acf_valid','status','cluster']]=[False,'zero_variance','']
        row=s.summarize_alleles(p,r,REGION)['allele_peak_features'].query("allele_group == 'REF'").iloc[0]
        self.assertEqual(row.n_valid,2)
        self.assertEqual(row.peak_fraction,1)

    def test_figures_have_no_sample_labels_and_acf_has_full_lags_and_two_zooms(self):
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib.figure import Figure
        from matplotlib.text import Text
        observed=[]
        def inspect(fig,path,**kwargs):
            observed.append(Path(path).name)
            labels=' '.join(text.get_text() for text in fig.findobj(Text))
            self.assertNotIn('NA18489',labels)
            self.assertNotIn('NA18508',labels)
            self.assertNotIn('cell line',labels.lower())
            if Path(path).name=='allele_average_acf.pdf':
                self.assertEqual(len(fig.axes),3)
                expected = [np.arange(width), np.arange(min(width, 501)), np.arange(500, min(width, 2001))]
                for ax, lags in zip(fig.axes, expected):
                    self.assertEqual(len(ax.lines),2)
                    self.assertEqual(len(ax.collections),0)
                    for line, allele in zip(ax.lines, ('REF','ALT')):
                        np.testing.assert_array_equal(line.get_xdata(),lags)
                        data=t['allele_acf_summary'].query('allele_group == @allele').set_index('lag_bp')
                        np.testing.assert_array_equal(line.get_ydata(),data.loc[lags,'mean_acf'].to_numpy())
                self.assertEqual(fig.axes[1].get_xlim(),(0,500))
                self.assertEqual(fig.axes[2].get_xlim(),(500,2000))
        for width in (300, 2000, 4000):
            with self.subTest(width=width):
                p,r=fixture([('NA18489_a','REF',1,.1),('NA18508_a','ALT',9,.2)],width=width)
                t=s.summarize_alleles(p,r,REGION)
                with tempfile.TemporaryDirectory() as directory, patch.object(Figure,'savefig',inspect):
                    s.plot_alleles(Path(directory),REGION,t)
        self.assertEqual(set(observed),{'allele_coverage.pdf','allele_cluster_proportions.pdf',
                                      'allele_average_acf.pdf','allele_peak_features.pdf'})


if __name__ == '__main__':
    unittest.main()
