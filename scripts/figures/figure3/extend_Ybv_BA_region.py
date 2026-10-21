#!/usr/bin/env python3
"""Extend comparative Yb BA across the displayed span using documented existing BAMs.
Same discovery/caller/quality policies as ybv_REF_ALT_rebuild_v1. No downloads,
alignment, focal calling or mixture fitting. Seven independent genomic shards.
"""
import argparse,hashlib,importlib.util,json,os,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
BASE=ROOT/'analysis/ybv_REF_ALT_rebuild_v1'
COHORT=ROOT/'analysis/ybv_short_ALT_geography_v7'
ENGINE=ROOT/'analysis/yb_completed_reads_recovery_v1/engine.py'
CHROM='NC_140403.1';LO=1630000;HI=1692000
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')
def loadmodule(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def config(out):return json.loads((out/'run.json').read_text())

def submit(parent):
    out=parent/'Ybv_BA_broad';out.mkdir(exist_ok=True)
    (out/'logs').mkdir(exist_ok=True);(out/'parts').mkdir(exist_ok=True)
    source=json.loads((BASE/'run.json').read_text())
    if json.loads((BASE/'complete.json').read_text())['status']!='COMPLETE':raise RuntimeError('Original Yb rebuild not complete')
    paths=[Path(source['reference']),ENGINE,COHORT/'membership_geography.tsv',BASE/'SNPs.tsv',BASE/'individual_ALT_dosages.tsv.gz']
    for p in paths:
        if not p.is_file():raise RuntimeError('Missing documented source '+str(p))
    bcftools=shutil.which('bcftools')
    if not bcftools:raise RuntimeError('bcftools missing in the active environment')
    bams=list(dict.fromkeys(source['bams']))
    if not bams:raise RuntimeError('No documented BAMs')
    for p in bams:
        if not Path(p).is_file():raise RuntimeError('Documented BAM unavailable: '+p)
    shards=[dict(start=s,end=min(s+9999,HI)) for s in range(LO,HI+1,10000)]
    cfg=dict(output=str(out),reference=source['reference'],bams=bams,bcftools=bcftools,shards=shards,
        cohort_sha256=sha(COHORT/'membership_geography.tsv'),source_config_sha256=sha(BASE/'run.json'),
        engine_sha256=sha(ENGINE),region=[CHROM,LO,HI],calling='mpileup q30 Q20 d10000; call -m -A -f GQ; existing bams',
        focal_calls_changed=False,alignments=False,downloads=False,reference_fai_sha256=sha(Path(source['reference']+'.fai')))
    save(out/'run.json',cfg);(out/'bams.list').write_text('\n'.join(bams)+'\n')
    shutil.copy2(__file__,out/'worker.py');shutil.copy2(ENGINE,out/'quality_engine.py')
    shutil.copy2(parent/'ba_worker.py',out/'ba_worker.py')
    (out/'complete.json').unlink(missing_ok=True)
    jobs={}
    def sb(stage,deps=None,array=None):
        wrap=shlex.join([sys.executable,str(out/'worker.py'),'--output',str(out),'--stage',stage])
        args=['sbatch','--parsable','--partition=cpu','--job-name=YbvBA_'+stage,
              '--cpus-per-task='+('4' if stage=='joint' else '1'),'--mem='+('16G' if stage!='preflight' else '8G'),
              '--time='+('04:00:00' if stage=='joint' else '01:00:00'),
              '--output='+str(out/'logs'/(stage+'.%A_%a.out')),'--error='+str(out/'logs'/(stage+'.%A_%a.err'))]
        if deps:args+=['--dependency=afterok:'+deps,'--kill-on-invalid-dep=yes']
        if array:args+=['--array='+array]
        j=subprocess.check_output(args+['--wrap='+wrap],text=True).strip().split(';')[0]
        jobs[stage]=j;save(out/'jobs.json',jobs);print('YB_BA_'+stage.upper()+'\t'+j,flush=True);return j
    p=sb('preflight');j=sb('joint',p,'0-'+str(len(shards)-1));return sb('score',j),jobs

def preflight(out):
    import pandas as pd,pysam
    c=config(out);members=pd.read_csv(COHORT/'membership_geography.tsv',sep='\t').fillna('')
    words=members.astype(str).agg(' | '.join,axis=1).str.lower()
    if words.str.contains(r'piano|bwpk|pk_dark|pk_pale|gynandromorph|tn_left|tn_right|tn2l|tn2r|tn3l|tn3r',regex=True).any():
        raise RuntimeError('PK/focal candidate entered the frozen comparative cohort')
    wanted={s for x in members.aliases for s in x.split(';') if s};seen=set();rows=[]
    with pysam.FastaFile(c['reference']) as f:
        dictionary=dict(zip(f.references,f.lengths))
        old=pd.read_csv(BASE/'SNPs.tsv',sep='\t')
        for r in old.itertuples():
            if f.fetch(CHROM,int(r.pos1)-1,int(r.pos1)).upper()!=r.ref:raise RuntimeError('Original Yb REF/FASTA mismatch')
    for path in c['bams']:
        with pysam.AlignmentFile(path,'rb') as b:
            sm={r.get('SM') for r in b.header.to_dict().get('RG',[]) if r.get('SM')}
            d=dict(zip(b.references,b.lengths))
            if not b.has_index() or len(sm)!=1 or CHROM not in d or any(dictionary.get(k)!=v for k,v in d.items()):raise RuntimeError('BAM reference/index/sample mismatch: '+path)
            name=next(iter(sm))
            if name in seen:raise RuntimeError('Duplicate BAM sample '+name)
            seen.add(name);p=Path(path);rows.append(dict(sample=name,path=str(p.resolve()),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns))
    if not wanted<=seen:raise RuntimeError('Cohort aliases absent from BAMs: '+','.join(sorted(wanted-seen)))
    pd.DataFrame(rows).to_csv(out/'validated_BAM_manifest.tsv',sep='\t',index=False)
    save(out/'preflight.complete.json',dict(status='COMPLETE',bam_manifest_sha256=sha(out/'validated_BAM_manifest.tsv'),samples=len(seen),cohort_aliases=len(wanted)))

def joint(out,task):
    c=config(out);span=c['shards'][task];part=out/'parts'/f'{task:02d}.vcf.gz';marker=part.with_suffix('.json')
    calling={k:c[k] for k in ['reference','bams','bcftools','shards','calling','reference_fai_sha256']}
    fingerprint=hashlib.sha256(json.dumps(calling,sort_keys=True).encode()).hexdigest()+sha(out/'validated_BAM_manifest.tsv')
    if marker.exists() and part.exists():
        saved=json.loads(marker.read_text())
        if saved.get('fingerprint')==fingerprint and saved.get('vcf_sha256')==sha(part):
            if not Path(str(part)+'.csi').is_file():subprocess.run([c['bcftools'],'index','-f',str(part)],check=True)
            return
    start=max(1,span['start']-200);end=span['end']+200;tmp=Path(str(part)+'.partial')
    mp=[c['bcftools'],'mpileup','--threads','2','-f',c['reference'],'-r',f'{CHROM}:{start}-{end}',
        '-b',str(out/'bams.list'),'-q','30','-Q','20','-d','10000','-a','FORMAT/AD,FORMAT/DP','-Ou']
    call=[c['bcftools'],'call','--threads','2','-m','-A','-f','GQ','-Oz','-o',str(tmp)]
    p=subprocess.Popen(mp,stdout=subprocess.PIPE)
    try:
        r=subprocess.run(call,stdin=p.stdout);p.stdout.close();rc=p.wait()
        if rc or r.returncode:raise RuntimeError('Broad Yb joint calling failed')
    finally:
        if p.poll() is None:p.kill();p.wait()
    tmp.replace(part);subprocess.run([c['bcftools'],'index','-f',str(part)],check=True)
    save(marker,dict(status='COMPLETE',core=span,fingerprint=fingerprint,vcf_sha256=sha(part),commands=[mp,call]))

def score(out):
    import numpy as np,pandas as pd,pysam
    c=config(out)
    if sha(COHORT/'membership_geography.tsv')!=c['cohort_sha256']:raise RuntimeError('Cohort changed during job')
    engine=loadmodule(out/'quality_engine.py','quality');ba=loadmodule(out/'ba_worker.py','ba')
    members=pd.read_csv(COHORT/'membership_geography.tsv',sep='\t').fillna('');members['v']=pd.to_numeric(members.v,errors='coerce')
    members=members.loc[members.v.isin([0,1])].copy();aliases=[str(x).split(';') for x in members.aliases]
    if members['sample'].duplicated().any():raise RuntimeError('Duplicate specimen representative')
    y=members.v.to_numpy(int);sites=[];vectors=[];seen=set()
    with pysam.FastaFile(c['reference']) as fa:
        for task,span in enumerate(c['shards']):
            part=out/'parts'/f'{task:02d}.vcf.gz'
            marker=json.loads(part.with_suffix('.json').read_text())
            if marker['status']!='COMPLETE':raise RuntimeError('Incomplete shard')
            with pysam.VariantFile(str(part)) as vf:
                for rec in vf.fetch(CHROM,span['start']-1,span['end']):
                    if len(rec.ref)!=1 or len(rec.alts or ())!=1 or len(rec.alts[0])!=1 or rec.ref not in 'ACGT' or rec.alts[0] not in 'ACGT':continue
                    key=(rec.pos,rec.ref,rec.alts[0])
                    if key in seen:raise RuntimeError('Duplicate core SNP')
                    seen.add(key)
                    if fa.fetch(CHROM,rec.pos-1,rec.pos).upper()!=rec.ref:raise RuntimeError('REF mismatch')
                    vals=[]
                    for aa in aliases:
                        calls=[engine.quality(rec.samples[s],allow_pl=True)[0] for s in aa]
                        calls=[v for v in calls if v is not None]
                        vals.append(calls[0] if calls and len(set(calls))==1 else np.nan)
                    if not np.isfinite(vals).any():continue
                    sites.append(dict(chrom=CHROM,pos=rec.pos,ref=rec.ref,alt=rec.alts[0]));vectors.append(vals)
    if not sites:raise RuntimeError('No usable broad Yb SNPs')
    matrix=np.asarray(vectors).T;sites=pd.DataFrame(sites)
    if sites.pos.duplicated().any():raise RuntimeError('Multiple biallelic variants at same coordinate')
    pd.DataFrame(matrix,index=members['sample'],columns=sites.pos).to_csv(out/'individual_ALT_dosages.tsv.gz',sep='\t')
    sites.to_csv(out/'SNPs.tsv',sep='\t',index=False);members.to_csv(out/'membership.tsv',sep='\t',index=False)
    # Check fresh broad calls against the already completed 3-kb cache.
    oldsites=pd.read_csv(BASE/'SNPs.tsv',sep='\t');old=pd.read_csv(BASE/'individual_ALT_dosages.tsv.gz',sep='\t',index_col=0)
    old=old.loc[members['sample']];lookup={int(r.pos):(j,r.ref,r.alt) for j,r in enumerate(sites.itertuples())}
    compared=discordant=shared=0
    for r in oldsites.itertuples():
        hit=lookup.get(int(r.pos1))
        if hit is None or (r.ref,r.alt)!=hit[1:]:continue
        shared+=1;a=old[str(r.pos1)].to_numpy(float);b=matrix[:,hit[0]];both=np.isfinite(a)&np.isfinite(b)
        compared+=int(both.sum());discordant+=int((a[both]!=b[both]).sum())
    concordance=dict(shared_SNPs=shared,compared_calls=compared,discordant_calls=discordant,discordance_fraction=discordant/compared if compared else None)
    save(out/'old_window_concordance.json',concordance)
    if not compared or discordant/compared>.01:raise RuntimeError('Broad calls require concordance review; previous SNP-panel calls remain untouched')
    rows=[]
    for j,r in enumerate(sites.to_dict('records')):
        for rule in ['carrier','homozygote']:rows.append(dict(**r,rule=rule,**ba.score(matrix[:,j],y,rule)))
    pd.DataFrame(rows).to_csv(out/'Ybv_all_SNP_BA.tsv.gz',sep='\t',index=False)
    save(out/'complete.json',dict(status='COMPLETE',cohort_sha256=c['cohort_sha256'],source_config_sha256=c['source_config_sha256'],
        span=[LO,HI],SNPs=len(sites),individuals=len(members),positive_n=int(sum(y==1)),negative_n=int(sum(y==0)),
        concordance=concordance,existing_alignments_reused=True,downloads=False,focal_calls_changed=False))
    print('YB_BROAD_BA_COMPLETE\t'+str(len(sites))+' SNPs; '+str(len(members))+' individuals',flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--stage',choices=['preflight','joint','score'],required=True);a=p.parse_args()
    if a.stage=='joint':joint(a.output,int(os.environ['SLURM_ARRAY_TASK_ID']))
    else:globals()[a.stage](a.output)
if __name__=='__main__':main()
