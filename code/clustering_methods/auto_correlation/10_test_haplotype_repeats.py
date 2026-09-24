"""Tests of annotation-only haplotype ACF aggregation, peaks and PDF output."""
import importlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

s = importlib.import_module('05_summarize_autocorrelations')
h = importlib.import_module('09_haplotype_repeat_plots')
f = importlib.import_module('06_test_allele_summaries')


class HaplotypeRepeats(unittest.TestCase):
    def test_periodic_signal_has_primary_and_supported_multiples(self):
        lag = np.arange(1000)
        profile = np.cos(2*np.pi*lag/180) * np.exp(-lag/2500)
        primary, higher = h.repeat_candidates(profile)
        self.assertEqual(primary,180)
        np.testing.assert_array_equal(higher,[360,540,720,900])

    def test_local_positive_prominent_peak_required_and_no_forced_period(self):
        for width in (200,251,1000):
            self.assertTrue(np.isnan(h.repeat_candidates(np.zeros(width))[0]))
        p = np.zeros(1000);p[180] = .009
        self.assertTrue(np.isnan(h.repeat_candidates(p)[0]))
        p[180] = .03;p[125] = .04;p[250] = .02;p[400] = .06
        primary,higher = h.repeat_candidates(p)
        self.assertEqual(primary,125)
        np.testing.assert_array_equal(higher,[250])
        p = np.full(1000,-.2);p[180] = -.01
        self.assertTrue(np.isnan(h.repeat_candidates(p)[0]))

    def test_median_includes_unclustered_fibers_and_is_not_mean(self):
        p,r=f.fixture([('NA1_a','REF',3,.1),('NA1_a','ALT',2,.2)])
        p[:3,180] = [.1,.2,.9]
        r.loc[2,['status','cluster']] = ['insufficient_reads','']
        before = p.copy(); metadata = r.copy(deep=True)
        rows=h.summarize_repeats(p,r,f.REGION,s.summarize_alleles(p,r,f.REGION))
        self.assertAlmostEqual(rows[0]['median'][180],.2)
        self.assertEqual(rows[0]['n_valid'],3)
        self.assertEqual(rows[0]['primary'],180)
        np.testing.assert_array_equal(p,before)
        pd.testing.assert_frame_equal(r,metadata)

    def test_invalid_missing_alleles_duplicates_and_order(self):
        p,r=f.fixture([('NA1_a','REF',4,.1)])
        p[0]=np.nan;r.loc[0,['acf_valid','status','cluster']] = [False,'zero_variance','']
        p[1]=0;p[1,0]=1;p[1,220]=.2
        p[2]=0;p[2,0]=1;p[2,150]=.2
        p[3]=0;p[3,0]=1
        duplicate=r.iloc[[2]].copy();duplicate['row_id']='NA1_b__'+duplicate.RID;duplicate['sample_name']='NA1_b'
        r=pd.concat([r,duplicate],ignore_index=True);p=np.vstack([p,p[2]])
        rows=h.summarize_repeats(p,r,f.REGION,s.summarize_alleles(p,r,f.REGION))
        self.assertEqual(rows[0]['n_invalid'],1)
        self.assertEqual(rows[0]['n_valid'],3)
        self.assertEqual(rows[0]['n_duplicates'],1)
        np.testing.assert_array_equal(rows[0]['indices'],[2,1,3])
        self.assertEqual(rows[1]['n_valid'],0)
        self.assertTrue(np.isnan(rows[1]['median']).all())

    def test_all_new_pdfs_and_matrix_values_without_tables(self):
        import matplotlib
        matplotlib.use('Agg')
        from matplotlib.figure import Figure
        p,r=f.fixture([('NA1_a','REF',3,.1),('NA1_a','ALT',2,.2)],width=1000)
        rows=h.summarize_repeats(p,r,f.REGION,s.summarize_alleles(p,r,f.REGION))
        observed=[]
        original=Figure.savefig
        def inspect(fig,path,**kwargs):
            name=Path(path).name
            if name=='locus_haplotype_median_acf.pdf':
                image=fig.axes[0].images[0].get_array()
                np.testing.assert_allclose(image,np.vstack([row['median'][120:251] for row in rows]))
            if name=='haplotype_median_acf.pdf':
                for i,row in enumerate(rows):
                    np.testing.assert_allclose(fig.axes[i*2].lines[0].get_ydata(),row['median'])
            observed.append(name)
            return original(fig,path,**kwargs)
        with tempfile.TemporaryDirectory() as directory, patch.object(Figure,'savefig',inspect):
            h.plot_haplotype_repeats(Path(directory),f.REGION,p,rows)
            self.assertEqual(len(list(Path(directory).glob('*.pdf'))),4)
            self.assertTrue(all(path.suffix=='.pdf' for path in Path(directory).iterdir()))
        self.assertEqual(set(observed),{'haplotype_median_acf.pdf','haplotype_fiber_acf_heatmap.pdf',
            'locus_haplotype_median_acf.pdf','haplotype_repeat_length_distribution.pdf'})

    def test_empty_and_short_acfs_render_without_inventing_peaks(self):
        for width, valid in [(3,True),(251,True),(300,False)]:
            with self.subTest(width=width,valid=valid):
                p,r=f.fixture([('NA1_a','REF',1,.1)],width=width)
                if not valid:
                    p[:]=np.nan;r.loc[:,['acf_valid','status','cluster']]=[False,'zero_variance','']
                rows=h.summarize_repeats(p,r,f.REGION,s.summarize_alleles(p,r,f.REGION))
                self.assertTrue(np.isnan(rows[0]['primary']))
                with tempfile.TemporaryDirectory() as directory:
                    h.plot_haplotype_repeats(Path(directory),f.REGION,p,rows)
                    self.assertEqual(len(list(Path(directory).glob('*.pdf'))),4)


if __name__=='__main__':
    unittest.main()
