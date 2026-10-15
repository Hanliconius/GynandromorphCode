#!/usr/bin/env python3
"""Append Optix QC, bcftools focal depths and the archived reciprocal model.
Does not edit frozen workers, cancel jobs, or change comparative SNP calling.
Completion here covers audited calls and mixture inference, not trait ascertainment.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import shlex
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

MODEL = '134_fit_reciprocal_L_R_mixture_genotypes.py'
FIELDS = ['chrom','pos','ref','alt','sample','raw_GT','ref_depth','alt_depth',
          'AD_total','DP','GQ','QC_pass']


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True)+'\n')


def run(cmd):
    subprocess.run([str(x) for x in cmd], check=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def biallelic(rec):
    return (len(rec.ref)==1 and rec.ref.upper() in 'ACGT' and rec.alts and
            len(rec.alts)==1 and len(rec.alts[0])==1 and rec.alts[0].upper() in 'ACGT')


def load_sites(cfg):
    import pysam
    path=Path(cfg['parent'])/'integration_Optix/Optix_window.vcf.gz'
    with pysam.VariantFile(str(path)) as vf:
        if set(vf.header.samples)!=set(cfg['comparative']):
            raise RuntimeError('Comparative sample identities changed')
        return {r.pos:(r.contig,r.ref,r.alts[0]) for r in vf if biallelic(r)}


def qc(cfg, task):
    import pysam
    out=Path(cfg['output'])
    parent=Path(cfg['parent'])
    totals={name:Counter() for name in cfg['comparative']}
    records=0
    with pysam.VariantFile(str(parent/'integration_Optix/Optix_window.vcf.gz')) as vf:
        if set(vf.header.samples)!=set(totals):raise RuntimeError('Cohort mismatch')
        for rec in vf:
            if not biallelic(rec):continue
            records+=1
            for name in totals:
                call=rec.samples[name]
                c=totals[name]
                gt=call.get('GT') or (None,None)
                dp=call.get('DP') or 0
                gq=call.get('GQ')
                called=all(x is not None for x in gt)
                c['sites']+=1
                c['DP_ge_2']+=dp>=2
                c['callable_DP2_GQ20']+=called and dp>=2 and gq is not None and gq>=20
                c['covered_called_GQ_missing']+=dp>=2 and called and gq is None
                c['covered_homref']+=dp>=2 and gt==(0,0)
                c['covered_homref_GQ_missing']+=dp>=2 and gt==(0,0) and gq is None
    if not records:raise RuntimeError('No biallelic Optix SNPs')
    rows=[]
    problems=[]
    for name,c in totals.items():
        rate=c['callable_DP2_GQ20']/records
        rows.append(dict(sample=name,**{key:c[key] for key in [
            'sites','DP_ge_2','callable_DP2_GQ20','covered_called_GQ_missing',
            'covered_homref','covered_homref_GQ_missing']},call_rate=rate,
            review_low_call_rate=rate<0.8))
        if c['covered_called_GQ_missing'] or not c['callable_DP2_GQ20']:
            problems.append(name)
    with (out/'comparative_callability.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
        w.writeheader();w.writerows(rows)
    dump(out/'comparative_qc.json',dict(passed=not problems,problem_samples=problems,
        sites=records,low_call_rate_review=[r['sample'] for r in rows if r['review_low_call_rate']],
        denominator='jointly discovered biallelic SNPs in Optix window; not invariant bases',
        gate='No covered called GT with missing GQ; each sample has callable sites',
        note='80% flag is review-only, not an imposed biological exclusion'))
    if problems:raise RuntimeError('Comparative GQ/callability audit failed: '+','.join(problems))
    print('OPTIX_COMPARATIVE_QC_PASS',records,flush=True)


def counts(cfg, task):
    import pysam
    out=Path(cfg['output'])
    name=cfg['focal'][task]
    row=next(r for r in cfg['samples'] if r['sample']==name)
    sites=load_sites(cfg)
    if not sites:raise RuntimeError('No Optix sites')
    chrom=next(iter(sites.values()))[0]
    if any(v[0]!=chrom for v in sites.values()):raise RuntimeError('Multiple Optix contigs')
    limits=cfg['limits']
    # Same bcftools executable, reference, MAPQ/BQ, default BAQ and per-BAM
    # depth limit (250) as the comparative mpileup. No diploid GT is inferred.
    target=out/(name+'.pileup.bcf')
    cmd=[cfg['bcftools'],'mpileup','-Ob','-o',str(target),'-f',cfg['fasta'],
         '-r',f'{chrom}:{max(1,limits[0]-200)}-{limits[1]+200}','-a','FORMAT/DP,FORMAT/AD',
         '-q','30','-Q','20','-d','250',row['bam']]
    run(cmd)
    run([cfg['bcftools'],'index','-f','-c',str(target)])
    found={}
    skipped_indels=0
    with pysam.VariantFile(str(target)) as vf:
        if list(vf.header.samples)!=[name]:raise RuntimeError('Focal BAM SM mismatch')
        for rec in vf:
            if rec.pos not in sites:continue
            # mpileup emits separate SNP and INDEL records at the same POS.
            # Only the SNP record contains the requested nucleotide counts.
            if rec.info.get('INDEL',False):
                skipped_indels+=1
                continue
            _,ref,alt=sites[rec.pos]
            if rec.ref.upper()!=ref.upper():
                raise RuntimeError(f'Focal/reference allele mismatch: {name} {rec.contig}:{rec.pos} pileup={rec.ref} joint={ref}')
            if rec.pos in found:raise RuntimeError(f'Duplicate SNP pileup at {rec.contig}:{rec.pos}')
            c=rec.samples[name]
            ad=c.get('AD') or ()
            by_allele={allele.upper():depth for allele,depth in zip(rec.alleles,ad)}
            nr=by_allele.get(ref.upper(),0) or 0
            na=by_allele.get(alt.upper(),0) or 0
            found[rec.pos]=(nr,na,c.get('DP') or 0,
                           sum(x or 0 for x in ad)-nr-na)
    dest=out/(name+'.counts.tsv.gz')
    with gzip.open(dest,'wt',newline='') as f:
        w=csv.writer(f,delimiter='\t')
        w.writerow(['chrom','pos','ref','alt','ref_depth','alt_depth','DP','other_depth'])
        for pos,(chrom,ref,alt) in sorted(sites.items()):
            w.writerow([chrom,pos,ref,alt,*found.get(pos,(0,0,0,0))])
    dump(out/(name+'.complete.json'),dict(sample=name,command=cmd,sites=len(sites),
         counts=str(dest),pileup=str(target),no_diploid_focal_call=True,
         skipped_overlapping_indel_records=skipped_indels))
    print('OPTIX_FOCAL_COUNTS',name,len(sites),flush=True)


def fit(cfg, task):
    import pysam
    out=Path(cfg['output'])
    if not json.loads((out/'comparative_qc.json').read_text())['passed']:
        raise RuntimeError('Comparative QC did not pass')
    sites=load_sites(cfg)
    pooled={side:{pos:[0,0,0] for pos in sites} for side in ('Left','Right')}
    for name in cfg['focal']:
        row=next(r for r in cfg['samples'] if r['sample']==name)
        json.loads((out/(name+'.complete.json')).read_text())
        with gzip.open(out/(name+'.counts.tsv.gz'),'rt') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                counts=pooled[row['side']][int(r['pos'])]
                for i,key in enumerate(['ref_depth','alt_depth','DP']):counts[i]+=int(r[key])
    source=out/'uniform_Optix_model_input.tsv.gz'
    with gzip.open(source,'wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS,delimiter='\t');w.writeheader()
        with pysam.VariantFile(str(Path(cfg['parent'])/'integration_Optix/Optix_window.vcf.gz')) as vf:
            for rec in vf:
                if not biallelic(rec):continue
                common=dict(chrom=rec.contig,pos=rec.pos,ref=rec.ref,alt=rec.alts[0])
                for name in cfg['comparative']:
                    c=rec.samples[name];ad=c.get('AD') or (0,0)
                    nr,na=(ad[0] or 0,ad[1] or 0) if len(ad)>=2 else (0,0)
                    gt=c.get('GT') or (None,None)
                    dp=c.get('DP') or 0;gq=c.get('GQ')
                    valid=all(v is not None for v in gt) and dp>=2 and gq is not None and gq>=20
                    w.writerow(dict(common,sample=name,raw_GT='/'.join('.' if v is None else str(v) for v in gt),
                        ref_depth=nr,alt_depth=na,AD_total=nr+na,DP=dp,GQ='' if gq is None else gq,QC_pass=valid))
                for side in ('Left','Right'):
                    nr,na,dp=pooled[side][rec.pos]
                    w.writerow(dict(common,sample='TN_'+side,raw_GT='./.',ref_depth=nr,
                                    alt_depth=na,AD_total=nr+na,DP=dp,GQ='',QC_pass=nr+na>=2))
    model_path=out/MODEL
    if sha(model_path)!=cfg['model_sha256']:raise RuntimeError('Frozen mixture model changed')
    spec=importlib.util.spec_from_file_location('archived_reciprocal_model',model_path)
    model=importlib.util.module_from_spec(spec);spec.loader.exec_module(model)
    model.INPUT=source;model.EXPECTED_SHA256=sha(source)
    old=model.OUTDIR
    model.OUTDIR=out/'reciprocal_mixture';model.OUTDIR.mkdir(exist_ok=True)
    for key,value in list(vars(model).items()):
        if isinstance(value,Path) and value.parent==old:
            setattr(model,key,model.OUTDIR/value.name)
    model.main()
    # The archived changed-from-original diagnostics are undefined without a
    # newly called naive pooled GT. Remove them rather than report all changed.
    import pandas as pd
    for key in ('ALL_OUT','P95_OUT'):
        p=getattr(model,key)
        df=pd.read_csv(p,sep='\t')
        df.drop(columns=[c for c in df if c.startswith('changed_from_original')],inplace=True)
        df.to_csv(p,sep='\t',index=False)
    summary=pd.read_csv(model.SUMMARY_OUT,sep='\t')
    summary=summary[~summary['statistic'].str.contains('changed_from_original')]
    summary.to_csv(model.SUMMARY_OUT,sep='\t',index=False)
    products=[model.ALL_OUT,model.P95_OUT,model.CENTRES_OUT,model.SUMMARY_OUT,
              model.MANIFEST_OUT,model.OPTIMIZER_OUT]
    model.CHECKSUM_OUT.write_text(''.join(f'{sha(p)}  {p.name}\n' for p in products))
    dump(out/'corrected_calls.complete.json',dict(status='AUDITED_CALLS_AND_MIXTURE_COMPLETE',
        comparative=len(cfg['comparative']),focal_libraries=cfg['focal'],
        model_sha256=cfg['model_sha256'],input_sha256=sha(source),
        parameter_policy='Refit error, bias, dispersion and reciprocal contamination from new data',
        focal_policy='T2+T3 pooled by side; no diploid GQ filter on focal evidence',
        ascertainment='All joint biallelic Optix SNPs used here; no phenotype SNP ranking yet',
        pending=['independent non-PK/non-focal Band/Dennis/Ray SNP selection','final labelled figures'],
        original_parent_complete_marker_is_not_corrected_analysis_complete=True))
    print('OPTIX_CORRECTED_CALLS_READY',out,flush=True)


def submit(args):
    import pysam, numpy, pandas, scipy  # Fail before submission if worker dependencies are absent.
    parent=args.run.resolve(strict=True)
    cfg=json.loads((parent/'run.json').read_text())
    jobs=json.loads((parent/'jobs.json').read_text())
    if cfg.get('no_hmel25') is not True or len(cfg['comparative'])!=42 or len(cfg['focal'])!=4:
        raise RuntimeError('Not the audited Rosina 42+4 workflow')
    if Counter(r['side'] for r in cfg['samples'] if r['sample'] in cfg['focal'])!={'Left':2,'Right':2}:
        raise RuntimeError('Expected two focal libraries on each side')
    model=args.model_script or Path(__file__).resolve().parents[2]/'reciprocal_mixing'/MODEL
    if not model.is_file():raise RuntimeError('Archived reciprocal model missing: '+str(model))
    out=parent/'optix_consistent_v2'
    if out.exists():raise RuntimeError('Already prepared; inspect '+str(out/'jobs.json'))
    out.mkdir();(out/'logs').mkdir()
    frozen=out/Path(__file__).name;shutil.copy2(__file__,frozen)
    shutil.copy2(model,out/MODEL)
    cfg.update(parent=str(parent),output=str(out),model_sha256=sha(out/MODEL),
               limits=[1060000,1240000])
    dump(out/'run.json',cfg)
    deps=[] if (parent/'complete.json').is_file() else jobs['extract']
    if isinstance(deps,str):deps=[deps]
    submitted={}
    def launch(stage,array,dependency,cpus,mem,time):
        worker=out/(stage+'.sh')
        worker.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec '+shlex.join([
            sys.executable,str(frozen),'--run',str(out),'--stage',stage])+ ' --task "$SLURM_ARRAY_TASK_ID"\n')
        worker.chmod(0o755)
        command=['sbatch','--parsable','--partition=cpu','--array='+array,
            '--cpus-per-task='+str(cpus),'--mem='+mem,'--time='+time,
            '--job-name=OptixV2_'+stage,
            '--output='+str(out/'logs'/(stage+'.%A_%a.out')),
            '--error='+str(out/'logs'/(stage+'.%A_%a.err')),str(worker)]
        if dependency:command.insert(-1,'--dependency=afterok:'+':'.join(dependency))
        jid=subprocess.check_output(command,text=True).strip().split(';')[0]
        submitted[stage]=jid;dump(out/'jobs.json',submitted)
        print('SUBMITTED',stage,jid,flush=True)
        return jid
    q=launch('qc','0',deps,1,'8G','01:00:00')
    c=launch('counts','0-3',deps,2,'12G','04:00:00')
    launch('fit','0',[q,c],2,'16G','08:00:00')
    print('OUTPUT',out,flush=True)
    print('EXISTING_JOBS_UNCHANGED; phenotype ascertainment and plotting still pending',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--submit',action='store_true')
    p.add_argument('--model-script',type=Path)
    p.add_argument('--stage',choices=['qc','counts','fit'])
    p.add_argument('--task',type=int,default=0)
    a=p.parse_args()
    if a.submit:submit(a)
    elif a.stage:
        cfg=json.loads((a.run/'run.json').read_text())
        globals()[a.stage](cfg,a.task)
    else:p.error('--submit or --stage is required')

if __name__=='__main__':main()
