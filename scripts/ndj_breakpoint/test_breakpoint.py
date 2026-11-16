import argparse
import hashlib
import json
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch
import sequence_features as sf
import run_breakpoint as rb
from score_junctions import assess


class Tests(unittest.TestCase):
    def test_union_coordinates(self):
        t=sf.Track([(3,8),(5,10),(10,12),(20,22)])
        self.assertEqual(t.coverage(0,100),11)
        self.assertEqual(t.coverage(8,21),5)
        self.assertEqual(t.coverage(12,20),0)
        self.assertEqual(t.intervals(8,21),[(8,12),(20,21)])

    def test_planted_motifs(self):
        g4='GGGAGGGAGGGAGGG'
        t=sf.scan('NN'+g4+'NN'+sf.revcomp(g4)+'NN'+'TTAGG'*4+'NN'+'AT'*10)
        self.assertGreaterEqual(t['G4_canonical'].coverage(0,100),30)
        self.assertEqual(t['telomere_TTAGG_3'].cumulative[-1],20)
        self.assertEqual(t['RY_alternating12_proxy'].coverage(58,100),20)
        self.assertEqual(t['ambiguous'].coverage(0,2),2)

    def test_seed_orientation_gap(self):
        arm='AGGCTATGCA'
        for key,right in [('direct10',arm),('inverted10',sf.revcomp(arm)),('mirror10',arm[::-1])]:
            s=arm+'TTC'+right
            hits=sf.repeat_pairs(s)
            self.assertIn((0,len(s)),hits[key])
        self.assertFalse(sf.repeat_pairs(arm+'N'*101+sf.revcomp(arm)).get('inverted10'))

    def test_no_ambiguous_structure(self):
        t=sf.scan('N'*100)
        self.assertTrue(all(v.cumulative[-1]==0 for k,v in t.items() if k!='ambiguous'))
        self.assertIsNone(sf.score('NN',{},0,2)['GC'])

    def test_low_complexity_bounded(self):
        t=sf.scan('G'*3000)
        self.assertEqual(t['G4_canonical'].coverage(0,3000),3000)
        self.assertEqual(t['homopolymer8'].coverage(0,3000),3000)

    def test_junction_requires_observed_sequence(self):
        self.assertEqual(assess('ACGT','GTAA','','ACGTAA')['microhomology_bp'],2)
        self.assertEqual(assess('ACGT','CCAA','TT','ACGTTTCCAA')['insertion_bp'],2)
        with self.assertRaises(ValueError):assess('ACGT','GTAA','','ACGTGG')

    def test_matched_background_and_exclusion(self):
        rng=random.Random(4);s=''.join(rng.choices('ACGT',k=300000))
        cfg=dict(start=145000,end=146000,padding=[0],controls=20,min_controls=10,gc_caliper=.03,valid_caliper=.02)
        target,groups=rb.background(s,{},145000,146000,cfg,random.Random(1))
        self.assertEqual(len(groups['same_chr']),20)
        self.assertEqual(len(groups['GC_valid_matched']),20)
        for r in groups['GC_valid_matched']:
            self.assertLessEqual(abs(r['GC']-target['GC']),.03)
            self.assertTrue(r['end']<=45000 or r['start']>=246000)
        rows=rb.compare(target,groups,cfg,'pad_0');rb.bh(rows)
        gc=next(r for r in rows if r['feature']=='GC' and r['background']=='GC_valid_matched')
        self.assertIsNone(gc['p_two_sided'])

    def test_bh_and_ties(self):
        target={'GC':.5,'valid_fraction':1,'x':0}
        groups={'same_chr':[target]*20}
        rows=rb.compare(target,groups,{'min_controls':10},'pad_0');rb.bh(rows)
        self.assertTrue(all(r['p_two_sided']==1 for r in rows))
        self.assertTrue(all(r['percentile']==50 for r in rows))
        self.assertTrue(all(r['q_BH_exploratory']==1 for r in rows))

    def test_synthetic_bam_and_worker(self):
        import pysam
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'out';out.mkdir()
            rng=random.Random(3);seq=''.join(rng.choices('ACGT',k=240000))
            ref=root/'ref.fa';ref.write_text('>'+rb.CHROM+'\n'+seq+'\n');pysam.faidx(str(ref))
            bams={}
            for name in ['NDJ_left','NDJ_right']:
                path=root/(name+'.bam')
                header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':rb.CHROM,'LN':len(seq),
                    'M5':hashlib.md5(seq.encode()).hexdigest()}],'RG':[{'ID':name,'SM':name}]}
                with pysam.AlignmentFile(str(path),'wb',header=header) as bam:
                    for i in range(4):
                        r=pysam.AlignedSegment();r.query_name='read'+str(i);r.query_sequence='A'*100
                        r.flag=0 if i<3 else 1024;r.reference_id=0;r.reference_start=120000+i
                        r.mapping_quality=60 if i!=2 else 0;r.cigar=[(4,20),(0,80)]
                        r.query_qualities=pysam.qualitystring_to_array('I'*100)
                        r.set_tag('SA',f'{rb.CHROM},180001,+,80M20S,60,0;');bam.write(r)
                pysam.index(str(path));bams[name]=str(path)
            source=root/'source';source.mkdir();fig=root/'figure4';fig.mkdir()
            rb.js(source/'config.json',dict(left='NDJ_left',right='NDJ_right',reference=str(ref),
                bam_left=bams['NDJ_left'],bam_right=bams['NDJ_right'],
                chromosomes=[dict(contig=rb.CHROM,length=len(seq))]))
            rb.js(source/'COMPLETE.json',dict(status='exploratory_scan_complete_not_biological_validation'))
            rb.js(fig/'run.json',dict(source=str(source)))
            rb.js(fig/'chr15_descriptive_segments.json',dict(transition_bracket_bp=[120000,121000]))
            args=argparse.Namespace(source=source,figure4=fig,project=root,interval=None,control_bam=[],
                annotations=None,controls=10,min_controls=10,gc_caliper=.03,padding=[0],seed=5,mapq=30,skip_reads=False)
            resolved=rb.resolve(args)
            self.assertEqual(resolved['start'],120000)
            self.assertEqual(resolved['bams'],bams)
            rb.js(fig/'run.json',dict(source=str(root/'wrong')))
            with self.assertRaises(ValueError):rb.resolve(args)
            cfg=dict(inputs=[],code=[],reference=str(ref),bams=bams,start=120000,end=121000,
                length=len(seq),padding=[0],controls=10,min_controls=10,gc_caliper=.03,valid_caliper=.02,
                seed=5,annotations=None,mapq=30,clip_min=20,skip_reads=False,origin='synthetic')
            rb.js(out/'run.json',cfg)
            # Plot generation belongs on Pegasus; verify worker tables here only.
            with patch.object(rb,'plot') as mocked:rb.work(out)
            self.assertTrue(mocked.called)
            self.assertEqual(rb.read_json(out/'complete.json')['status'],'exploratory_context_complete')
            text=(out/'clip_clusters.tsv').read_text()
            self.assertIn('120000\tleft\t1',text)
            self.assertNotIn('120002\tleft',text)
            self.assertNotIn('120003\tleft',text)
            self.assertTrue((out/'homology_queries.fa').is_file())


if __name__=='__main__':unittest.main()
