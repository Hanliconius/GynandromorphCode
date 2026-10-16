import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import yb_uniform_frequency as m


class UniformYbTests(unittest.TestCase):
    def test_missing_GT_and_quality_are_never_reference(self):
        for call in [dict(GT=(0,0),DP=20),dict(GT=(0,None),DP=20,GQ=50),
                     dict(GT=(0,1),DP=1,GQ=50),dict(GT=(0,2),DP=20,GQ=50)]:
            self.assertTrue(np.isnan(m.valid_dosage(call)))
        self.assertEqual(m.valid_dosage(dict(GT=(0,1),DP=2,GQ=20)),.5)

    def test_control_scoring_and_REF_orientation(self):
        y=np.array([0,0,1,1])
        result=m.score(np.array([1,1,0,.5]),y)
        self.assertEqual(result['balanced_accuracy'],1)
        self.assertFalse(result['pattern_is_alt'])
        self.assertIsNone(m.score(np.array([1,np.nan,0,.5]),y))

    def test_group_frequency_is_equal_individual_weight_and_missing_gate(self):
        self.assertEqual(m.group_mean([0,.5,1]),(.5,3))
        freq,n=m.group_mean([0,.5,1,np.nan,np.nan])
        self.assertTrue(np.isnan(freq));self.assertEqual(n,3)

    def test_user_phenotypes_override_metadata_and_PK_not_in_controls(self):
        cfg=dict(comparative=['c','v','p'],samples=[
            dict(sample='c',label='c',group='cythera'),
            dict(sample='v',label='v',group='vulcanus'),
            dict(sample='p',label='p',group='PK_dark')])
        with tempfile.TemporaryDirectory() as d:
            out=Path(d);meta=out/'meta.tsv'
            pd.DataFrame([dict(Accession='c',**{'Yb-dorsal':1,'Yb-ventral':1}),
                          dict(Accession='p',**{'Yb-dorsal':1,'Yb-ventral':1})]).to_csv(meta,sep='\t',index=False)
            t=m.phenotype_manifest(cfg,meta,out)
            np.testing.assert_array_equal(t[['Ybd','Ybv']].to_numpy(),[[0,0],[0,1],[1,1]])
            groups=m.group_rows(t,'Ybv',['c','v','p'])
            self.assertEqual(groups[0]['samples'],['v'])
            self.assertEqual(groups[1]['samples'],['c'])
            self.assertEqual(groups[2]['samples'],['p'])

    def test_posterior_calls_retain_opposing_states_and_low_depth_missing(self):
        modelpath=Path(__file__).with_name(m.MODEL)
        if not modelpath.is_file():
            modelpath=Path(__file__).resolve().parents[2]/'scripts/reciprocal_mixing'/m.MODEL
        model=m.load_module(modelpath,'test_original')
        # Nearly noiseless synthetic depths from the model's L=0, R=1 state.
        df=pd.DataFrame(dict(L_ref=[800,1],L_alt=[200,0],R_ref=[200,1],R_alt=[800,0]))
        result=m.posterior_calls(model,df,[.2,.2,.0005,0,5000])
        self.assertEqual(result.loc[0,'left_dosage'],0)
        self.assertEqual(result.loc[0,'right_dosage'],1)
        self.assertTrue(np.isnan(result.loc[1,'left_dosage']))
        self.assertEqual(result.loc[1,'right_GT'],'./.')


if __name__=='__main__':unittest.main()
