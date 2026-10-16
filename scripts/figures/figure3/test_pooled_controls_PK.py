import unittest
import numpy as np
import pandas as pd
import plot_optix_pooled_controls_PK as m

class PoolTests(unittest.TestCase):
    def setUp(self):
        self.meta=pd.DataFrame([dict(sample=n,group=g,B=b,D=d,R=0,is_PK=pk,
            Band_selection_eligible=not pk,Dennis_selection_eligible=not pk,Ray_selection_eligible=not pk)
            for n,g,b,d,pk in [('a','rosina',1,0,False),('b','cythera',1,0,False),
            ('c','vulcanus',1,0,False),('d','florencia',0,1,False),
            ('p','PK_dark',1,0,True),('q','PK_dark',1,1,True),('r','PK_pale',1,0,True)]])
        self.names=self.meta['sample'].tolist()
    def test_one_control_row_each_and_no_double_counting(self):
        groups=m.frequency_groups(self.meta,self.names,'Band')
        self.assertEqual(groups[0]['samples'],['a','b','c'])
        self.assertEqual(groups[-1]['samples'],['d'])
        self.assertEqual([g['samples'] for g in groups[1:-1]],[['p'],['q'],['r']])
        self.assertEqual(sum(g['n'] for g in groups),len(self.names))
    def test_pools_individuals_not_means_of_races(self):
        groups=m.frequency_groups(self.meta,self.names,'Band')
        f,n=m.group_frequency(np.array([0,0,1,1,1,1,1])[groups[0]['indices']])
        self.assertAlmostEqual(f,1/3);self.assertEqual(n,3)
    def test_membership_changes_by_trait_but_PK_remains_separate(self):
        groups=m.frequency_groups(self.meta,self.names,'Dennis')
        self.assertEqual(groups[0]['samples'],['d'])
        self.assertEqual(groups[-1]['samples'],['a','b','c'])
        self.assertTrue(all(g['is_PK'] for g in groups[1:-1]))
    def test_duplicate_individual_fails(self):
        with self.assertRaisesRegex(RuntimeError,'Duplicated'):
            m.frequency_groups(pd.concat([self.meta,self.meta.iloc[:1]]),self.names,'Band')
    def test_saved_matrix_order_matches_control_focal_PK_labels(self):
        groups=m.frequency_groups(self.meta,self.names,'Band')
        t=pd.DataFrame([dict(chrom='c',pos=1,ref='A',alt='G',site_index=0,pattern_is_alt=True)])
        recall={('c',1,'A','G'):dict(mixture_left_GT='0/0',mixture_right_GT='0/1')}
        matrix,_=m.frequency_matrix(t,np.array([[1,1,1,0,0,.5,1]]),recall,groups)
        order=[2,0,1,*range(3,len(groups)+2)]
        np.testing.assert_allclose(matrix[order,0],[1,0,.5,0,.5,1,0])

if __name__=='__main__':unittest.main()
