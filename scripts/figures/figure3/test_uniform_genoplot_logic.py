import importlib.util
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).parent
spec=importlib.util.spec_from_file_location('plotter',HERE/'plot_optix_uniform_genoplots.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class LogicTests(unittest.TestCase):
    def test_namedtuple_labels_for_controls_PK_and_missing_annotations(self):
        rows=pd.DataFrame([
            dict(sample='ERR1',label='control',group='rosina',is_PK=False,PK_colour=np.nan,B=0,D=1,R=1),
            dict(sample='SRR1',label='P1001',group='PK_dark',is_PK=True,PK_colour='dark',B=1,D=0,R=0),
            dict(sample='ERR2',label='new',group='cythera',is_PK=False,PK_colour=np.nan,B=np.nan,D=np.nan,R=np.nan)])
        labels=[m.individual_label(r) for r in rows.itertuples()]
        self.assertEqual(labels[0],'control rosina | B0 D1 R1 | ERR1')
        self.assertEqual(labels[1],'P1001 PK dark | B1 D0 R0 | SRR1')
        self.assertEqual(labels[2],'new cythera | B? D? R? | ERR2')
    def test_orientation_and_dominant_heterozygotes(self):
        self.assertEqual(m.score([0,0,.5,1],[0,0,1,1])['balanced_accuracy'],1)
        self.assertTrue(m.score([0,0,.5,1],[0,0,1,1])['pattern_is_alt'])
        self.assertFalse(m.score([1,1,.5,0],[0,0,1,1])['pattern_is_alt'])
    def test_missing_calls_and_unidentifiable_orientation(self):
        self.assertIsNone(m.score([0,np.nan,.5,1],[0,0,1,1]))
        self.assertIsNone(m.score([.5,.5,.5,.5],[0,0,1,1]))
    def test_display_retains_conflicts_and_does_not_rank_by_focal_GT(self):
        scores=pd.DataFrame([dict(chrom='chr18',pos=1112000+i,ref='A',alt='G',module='Band',
            balanced_accuracy=.95,site_index=i,pattern_is_alt=True) for i in range(10)])
        recall={tuple(r[k] for k in m.KEY):dict(pair_posterior=.99,
            mixture_left_GT='1/1',mixture_right_GT='0/0') for r in scores.to_dict('records')}
        chosen,counts=m.display_sites(scores,recall,.9,True,'Band',5)
        self.assertEqual(counts['focal_callable'],10)
        self.assertEqual(chosen.pos.tolist(),[1112000,1112002,1112004,1112007,1112009])
        self.assertEqual(counts['displayed'],5)
    def test_PK_never_used_and_unknown_race_not_inferred(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            meta=root/'data/derived';meta.mkdir(parents=True)
            (meta/'region_panel_metadata.tsv').write_text('Accession\tTaxon Name\tB\tD\tR\nERR1\trosina\t0\t1\t1\n')
            old=root/'analysis/figures/key_regions_population_context_genoplots';old.mkdir(parents=True)
            (old/'exact_BA_comparative_sample_manifest.tsv').write_text('sample\tBA_pattern\tphenotype\nERR1\tBand\t0\nERR1\tDennis\t1\nERR1\tRay\t1\nPK1\tBand\t1\n')
            band=root/'analysis/original_502/checkpoints';band.mkdir(parents=True)
            pd.DataFrame([{'sample':'ERR1','plot_group':'Band absent (B=0)'}]).to_csv(band/'optix_band_original_region_reciprocal_genoplot.long.tsv.gz',sep='\t',index=False)
            pk=root/'pk.tsv';pk.write_text('sample\tlabel\tphenotype\tB\tD\tR\nPK1\tP1001\tdark\t1\t0\t0\n')
            cfg=dict(comparative=['ERR1','PK1','newRace'],samples=[dict(sample=n,label=n,group=g)
                for n,g in [('ERR1','rosina'),('PK1','PK_dark'),('newRace','cythera')]])
            table,_=m.phenotype_manifest(cfg,root,pk)
            a=table.set_index('sample')
            self.assertTrue(a.loc['ERR1','Band_selection_eligible'])
            self.assertFalse(a.loc['PK1','Band_selection_eligible'])
            self.assertFalse(a.loc['newRace','Dennis_selection_eligible'])
            self.assertTrue(pd.isna(a.loc['newRace','D']))

if __name__=='__main__':unittest.main()
