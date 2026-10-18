import unittest
import numpy as np
import optix_arnaud_review as m

class ReviewTests(unittest.TestCase):
    def test_ref_alt_orientation(self):
        r={'mixture_left_GT':'1/1','mixture_right_GT':'0/1'}
        self.assertEqual(m.pair_doses(r,False),[0,.5])
        self.assertEqual(m.pair_doses(r,True),[1,.5])
    def test_unknown_focal_stays_missing(self):
        self.assertTrue(np.isnan(m.pair_doses({},True)).all())
    def test_confounded_classes_not_scored(self):
        results=m.conditional_scores(np.array([0,0,1,1]),np.array([0,0,1,1]),np.array([0,0,1,1]))
        self.assertTrue(all('balanced_accuracy' not in x for x in results))
    def test_conditional_signal_and_missingness(self):
        y=np.array([0,0,1,1,0,0,1,1]);z=np.array([0,0,0,0,1,1,1,1])
        r=m.conditional_scores(y.astype(float),y,z)
        self.assertEqual([x['balanced_accuracy'] for x in r],[1,1])
        d=y.astype(float);d[0]=np.nan
        self.assertNotIn('balanced_accuracy',m.conditional_scores(d,y,z)[0])
if __name__=='__main__':unittest.main()
