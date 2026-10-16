import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
import plot_optix_frequency_genoplots as m

class FrequencyTests(unittest.TestCase):
    def test_equal_individual_frequency_and_missing_call_threshold(self):
        f,n=m.group_frequency([0,.5,1])
        self.assertEqual(f,.5);self.assertEqual(n,3)
        f,n=m.group_frequency([0,.5,1,np.nan,np.nan])
        self.assertTrue(np.isnan(f));self.assertEqual(n,3)
        f,n=m.group_frequency([0,.5,1,.5,np.nan])
        self.assertEqual(f,.5);self.assertEqual(n,4)
    def test_race_PK_colour_and_phenotypes_stay_separate(self):
        meta=pd.DataFrame([dict(sample=n,group=g,B=1,D=d,R=0,is_PK=p) for n,g,d,p in
            [('a','cythera',0,False),('b','cythera',0,False),('c','PK_dark',0,True),
             ('d','PK_dark',1,True),('e','PK_pale',0,True)]])
        groups=m.frequency_groups(meta,['a','b','c','d','e'])
        self.assertEqual([g['n'] for g in groups],[2,1,1,1])
        self.assertEqual(len(groups),4)
    def test_user_annotations_and_Yb_are_explicit(self):
        meta=pd.DataFrame([dict(sample=n,group=g,B=np.nan,D=np.nan,R=np.nan,is_PK=False,
            Band_selection_eligible=False,Dennis_selection_eligible=False,Ray_selection_eligible=False)
            for n,g in [('a','cythera'),('b','vulcanus')]])
        with patch.object(m,'ORIGINAL_MANIFEST',return_value=(meta,[])):
            new,_=m.annotated_manifest({},None,None)
        np.testing.assert_array_equal(new[['B','D','R']].to_numpy(),[[1,0,0],[1,0,0]])
        self.assertEqual(new.iloc[1]['Yb_dorsal'],0)
        self.assertEqual(new.iloc[1]['Yb_ventral'],1)
        self.assertTrue(np.isnan(new.iloc[0]['Yb_dorsal']))
        self.assertTrue(new['Band_selection_eligible'].all())
    def test_REF_orientation_applies_to_group_and_corrected_focal(self):
        sites=pd.DataFrame([dict(chrom='chr18',pos=10,ref='A',alt='G',site_index=0,pattern_is_alt=False)])
        recall={('chr18',10,'A','G'):dict(mixture_left_GT='0/0',mixture_right_GT='0/1')}
        groups=[dict(indices=[0,1],n=2)]
        matrix,called=m.frequency_matrix(sites,np.array([[0,.5]]),recall,groups)
        np.testing.assert_allclose(matrix[:,0],[1,.5,.75])
        self.assertEqual(called[0,0],2)

if __name__=='__main__':unittest.main()
