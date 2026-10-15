import argparse
import csv
import gzip
import importlib.util
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam
import pysam.bcftools

HERE=Path(__file__).parent
spec=importlib.util.spec_from_file_location('post',HERE/'extend_optix_consistent_postcall.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class PostcallTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.parent=self.root/'parent';self.out=self.root/'output'
        (self.parent/'integration_Optix').mkdir(parents=True);self.out.mkdir()
        self.cfg=dict(parent=str(self.parent),output=str(self.out),comparative=['control'],
            focal=['L2','L3','R2','R3'],samples=[dict(sample=n,side=s,bam=n+'.bam')
                for n,s in [('L2','Left'),('L3','Left'),('R2','Right'),('R3','Right')]],
            fasta='reference.fa',bcftools='bcftools',limits=[1,2000])
    def tearDown(self):self.tmp.cleanup()
    def vcf(self,n=1,missing_gq=False):
        h=pysam.VariantHeader();h.contigs.add('NC_140406.1',length=5000)
        for key,typ,num in [('GT','String',1),('DP','Integer',1),('GQ','Integer',1),('AD','Integer','R')]:
            h.formats.add(key,num,typ,key)
        h.add_sample('control')
        p=self.parent/'integration_Optix/Optix_window.vcf.gz'
        with pysam.VariantFile(str(p),'wz',header=h) as f:
            for i in range(n):
                r=f.new_record(contig='NC_140406.1',start=i+9,alleles=('A','G'))
                gt=[(0,0),(0,1),(1,1)][i%3];ad=[(40,0),(20,20),(0,40)][i%3]
                r.samples['control']['GT']=gt;r.samples['control']['DP']=40
                r.samples['control']['AD']=ad
                if not missing_gq:r.samples['control']['GQ']=60
                f.write(r)
    def test_completed_parent_avoids_stale_slurm_dependencies(self):
        cfg=dict(self.cfg,no_hmel25=True,comparative=['c'+str(i) for i in range(42)])
        (self.parent/'run.json').write_text(json.dumps(cfg))
        (self.parent/'jobs.json').write_text(json.dumps({'extract':['old_job']}))
        (self.parent/'complete.json').write_text('{}')
        commands=[]
        def submit_command(cmd,**kwargs):
            commands.append(cmd);return str(100+len(commands))+'\n'
        args=argparse.Namespace(run=self.parent,model_script=HERE/m.MODEL)
        with patch.object(m.subprocess,'check_output',side_effect=submit_command):m.submit(args)
        self.assertEqual(len(commands),3)
        self.assertFalse(any('--dependency=' in x for cmd in commands[:2] for x in cmd))
        self.assertIn('--dependency=afterok:101:102',commands[2])
        self.assertIn('--array=0-3',commands[1])
    def test_missing_reference_gq_blocks(self):
        self.vcf(missing_gq=True)
        with self.assertRaisesRegex(RuntimeError,'audit failed'):m.qc(self.cfg,0)
        report=json.loads((self.out/'comparative_qc.json').read_text())
        self.assertFalse(report['passed'])
    def test_allele_identity_not_alt_index(self):
        self.vcf()
        commands=[]
        def fake(cmd):
            commands.append(cmd)
            if cmd[1]=='index':
                pysam.bcftools.index('-f','-c',cmd[-1])
                return
            h=pysam.VariantHeader();h.contigs.add('NC_140406.1',length=5000)
            h.formats.add('DP',1,'Integer','depth');h.formats.add('AD','R','Integer','depth')
            h.add_sample('L2')
            h.info.add('INDEL',0,'Flag','Indel record')
            with pysam.VariantFile(str(self.out/'L2.pileup.bcf'),'wb',header=h) as f:
                r=f.new_record(contig='NC_140406.1',start=9,alleles=('A','C','G'))
                r.samples['L2']['AD']=(11,3,7);r.samples['L2']['DP']=21;f.write(r)
                r=f.new_record(contig='NC_140406.1',start=9,alleles=('AT','A'))
                r.info['INDEL']=True;r.samples['L2']['AD']=(13,8);r.samples['L2']['DP']=21;f.write(r)
        with patch.object(m,'run',side_effect=fake):m.counts(self.cfg,0)
        with gzip.open(self.out/'L2.counts.tsv.gz','rt') as f:r=next(csv.DictReader(f,delimiter='\t'))
        self.assertEqual([r[k] for k in ['ref_depth','alt_depth','other_depth']],['11','7','3'])
        self.assertIn('-d',commands[0]);self.assertIn('250',commands[0])
        self.assertEqual(json.loads((self.out/'L2.complete.json').read_text())['skipped_overlapping_indel_records'],1)
    def test_real_archived_model_refits_and_pools(self):
        self.vcf(n=540);m.qc(self.cfg,0)
        shutil.copy2((HERE/m.MODEL) if (HERE/m.MODEL).is_file() else HERE.parents[1]/'reciprocal_mixing'/m.MODEL,self.out/m.MODEL)
        self.cfg['model_sha256']=m.sha(self.out/m.MODEL)
        for j,name in enumerate(self.cfg['focal']):
            with gzip.open(self.out/(name+'.counts.tsv.gz'),'wt',newline='') as f:
                w=csv.writer(f,delimiter='\t');w.writerow(['chrom','pos','ref','alt','ref_depth','alt_depth','DP','other_depth'])
                for i in range(540):
                    l=(i%9)//3/2;r=(i%9)%3/2
                    af=.8*l+.2*r if name.startswith('L') else .25*l+.75*r
                    na=round(40*af)
                    w.writerow(['NC_140406.1',i+10,'A','G',40-na,na,40,0])
            (self.out/(name+'.complete.json')).write_text('{}')
        m.fit(self.cfg,0)
        p=self.out/'reciprocal_mixture/optix_reciprocal_L_R_mixture.summary.tsv'
        with p.open() as f:stats={r['statistic']:float(r['value']) for r in csv.DictReader(f,delimiter='\t')}
        self.assertAlmostEqual(stats['estimated_c_left_from_right'],.2,delta=.025)
        self.assertAlmostEqual(stats['estimated_c_right_from_left'],.25,delta=.025)
        with gzip.open(self.out/'uniform_Optix_model_input.tsv.gz','rt') as f:
            r=next(r for r in csv.DictReader(f,delimiter='\t') if r['sample']=='TN_Left')
        self.assertEqual(int(r['DP']),80)
        self.assertNotIn('posterior95_changed_from_original',stats)
        self.assertEqual(json.loads((self.out/'corrected_calls.complete.json').read_text())['status'],
                         'AUDITED_CALLS_AND_MIXTURE_COMPLETE')

if __name__=='__main__':unittest.main()
