import tempfile
import unittest
from pathlib import Path
import numpy as np
import yb_previous_calls_new_samples as m

class PreviousDisplayTests(unittest.TestCase):
    def test_original_unmixed_PK_call_orientation_and_depth(self):
        path=Path(__file__).with_name('134_fit_reciprocal_L_R_mixture_genotypes.py')
        if not path.is_file():path=Path(__file__).resolve().parents[2]/'scripts/reciprocal_mixing'/path.name
        model=m.module(path,'test_model')
        p=[.2086930887,.2098653549,.0005,-.1195160302,56.99495196]
        value,confidence=m.pk_dosage(model,p,100,0,True)
        self.assertEqual(value,1);self.assertGreaterEqual(confidence,.95)
        value,confidence=m.pk_dosage(model,p,100,0,False)
        self.assertEqual(value,0)
        self.assertTrue(np.isnan(m.pk_dosage(model,p,1,0,False)[0]))
    def test_race_display_retains_original_available_individual_mean(self):
        np.testing.assert_equal(m.mean_available([[0,.5,np.nan],[1,np.nan,np.nan]]),[.5,.5,np.nan])
    def test_baseline_requires_actual_completed_bundle(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);log=d/'run.out';log.write_text('BUNDLE_PDF\t'+str(d/'missing.pdf')+'\n')
            with self.assertRaises(RuntimeError):m.parse_baseline(log)

if __name__=='__main__':unittest.main()
