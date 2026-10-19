#!/usr/bin/env python3
"""Audit meriana identity and optix genotype sharing without changing calls or selection."""
import argparse,csv,gzip,hashlib,json,re,shlex,shutil,subprocess,sys,tempfile
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
RUN=ROOT/'analysis/rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc'
SOURCES={'legacy':ROOT/'external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf','merged':RUN/'Optix_legacy_merge_v1/Optix_legacy_plus_recalled.vcf.gz'}
def table(path,rows,cols=None):
 if cols is None:cols=list(dict.fromkeys(k for r in rows for k in r))
 with (gzip.open(path,'wt',newline='') if path.suffix=='.gz' else path.open('w',newline='')) as h:
  w=csv.DictWriter(h,delimiter='\t',fieldnames=cols or ['status']);w.writeheader();w.writerows(rows)
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def canonical(gt):
 return tuple(sorted(gt)) if gt and len(gt)==2 and all(x in (0,1) for x in gt) else None
def dosage(c):
 import math
 gt=canonical(c.get('GT'));dp=c.get('DP') or 0;q=c.get('GQ');source='GQ'
 if gt is None or dp<2:return None,'missing_or_DP_below_2'
 if q is None:
  pl=c.get('PL')
  if not pl or len(pl)!=3 or any(x is None for x in pl):return None,'GQ_and_biallelic_PL_missing'
  w=[10**(-(x-min(pl))/10) for x in pl];state=sum(gt)
  if w[state]<max(w)-1e-10:return None,'GT_PL_disagree'
  q=-10*math.log10(max(1e-99,1-w[state]/sum(w)));source='PL_derived_likelihood_Q'
 return (sum(gt) if q>=20 else None),source if q>=20 else 'quality_below_20'
def do_audit(a):
 import pysam
 out=a.output;out.mkdir(parents=True,exist_ok=False)
 metadata=ROOT/'data/derived/region_panel_metadata.tsv'
 with metadata.open() as h:meta=list(csv.DictReader(h,delimiter='\t'))
 individuals={};manifest=[]
 for number,r in enumerate(meta,2):
  taxon=r.get('Taxon Name','');normal=' '.join(taxon.split())
  group='meriana_quoted' if 'meriana' in normal.lower() and re.search(r"['\"‘’“”]",normal) else 'meriana_unquoted' if 'meriana' in normal.lower() else 'besckei' if 'besckei' in normal.lower() else 'timareta' if 'timareta' in normal.lower() else ''
  if not group:continue
  names=set()
  for col in ('Accession','SequenceID','EarthCapeID'):names.update(re.findall(r'(?:ERR|SRR|DRR)\d+',r.get(col,'')))
  for s in names:
   entry=dict(sample=s,group=group,raw_taxon=taxon,normal_taxon=normal,Country=r.get('Country',''),B=r.get('B',''),D=r.get('D',''),R=r.get('R',''),metadata_line=number)
   if s in individuals and individuals[s]!=entry:raise RuntimeError('Ambiguous duplicate sample metadata: '+s)
   individuals[s]=entry;manifest.append(entry)
 table(out/'sample_identities.tsv',manifest)
 calls=[];summaries=[];absent=[];comparison={};qc=Counter()
 for label,path in SOURCES.items():
  if not path.is_file():absent.append(dict(source=label,path=str(path),reason='file_not_present'));continue
  with pysam.VariantFile(str(path)) as vf:
   samples=sorted(set(vf.header.samples)&set(individuals));counts=Counter(individuals[s]['group'] for s in samples)
   for s in set(individuals)-set(samples):absent.append(dict(source=label,sample=s,reason='sample_absent_from_VCF'))
   for rec in vf:
    if not 1060000<=rec.pos<=1240000 or len(rec.ref)!=1 or len(rec.alts or ())!=1 or len(rec.alts[0])!=1:continue
    interval='Dennis1' if 1215210<=rec.pos<=1216124 else 'Dennis2' if 1218557<=rec.pos<=1221494 else 'outside_published_D_intervals'
    key=(rec.pos,rec.ref,rec.alts[0]);groups=defaultdict(list);races=defaultdict(list)
    for s in samples:
     c=rec.samples[s];value,quality=dosage(c);qc[(label,quality)]+=1
     if label=='legacy':comparison[(s,key)]=canonical(c.get('GT'))
     elif (s,key) in comparison:
      before=comparison[(s,key)];after=canonical(c.get('GT'))
      if before is not None and after is not None:
       qc[('overlapping_raw_GT','compared')]+=1
       if before!=after:qc[('overlapping_raw_GT','discordant')]+=1
     groups[individuals[s]['group']].append(value);races[individuals[s]['normal_taxon']].append(value)
     if interval!='outside_published_D_intervals' or a.position==rec.pos:
      gt=c.get('GT') or ()
      calls.append(dict(source=label,chrom=rec.contig,pos=rec.pos,ref=rec.ref,alt=rec.alts[0],interval=interval,sample=s,group=individuals[s]['group'],taxon=individuals[s]['normal_taxon'],raw_GT='/'.join('.' if x is None else str(x) for x in gt),DP=c.get('DP'),GQ=c.get('GQ'),AD=str(c.get('AD')),PL=str(c.get('PL')),qualified_alt_dosage=value,quality_source=quality))
    stats={}
    for group in ('meriana_unquoted','meriana_quoted','besckei','timareta'):
     vals=[v for v in groups[group] if v is not None];stats[group+'_n_called']=len(vals);stats[group+'_alt_AF']=sum(vals)/(2*len(vals)) if vals else None
    true=[x for x in groups['meriana_unquoted'] if x is not None]
    same_alt=bool(true) and all(x>0 for x in true) and any(x is not None and x>0 for x in groups['besckei']) and any(x is not None and x>0 for x in groups['timareta'])
    summaries.append(dict(source=label,chrom=rec.contig,pos=rec.pos,ref=rec.ref,alt=rec.alts[0],interval=interval,all_callable_unquoted_meriana_share_ALT_with_besckei_and_timareta=same_alt,**stats))
 table(out/'Dennis_individual_raw_calls.tsv.gz',calls)
 if a.position:table(out/'requested_position_raw_calls.tsv',[r for r in calls if r['pos']==a.position])
 table(out/'Dennis_site_summary.tsv',[r for r in summaries if r['interval']!='outside_published_D_intervals'])
 candidates=[r for r in summaries if r['all_callable_unquoted_meriana_share_ALT_with_besckei_and_timareta']]
 table(out/'shared_ALT_candidate_sites.tsv',candidates)
 if a.position:table(out/'requested_position_summary.tsv',[r for r in summaries if r['pos']==a.position])
 table(out/'missing_sources_or_samples.tsv',absent)
 table(out/'QC_summary.tsv',[dict(source=k[0],metric=k[1],count=v) for k,v in sorted(qc.items())])
 summary=dict(status='AUDIT_COMPLETE_REQUIRES_IDENTITY_AND_SITE_REVIEW',metadata_sha256=hashlib.sha256(metadata.read_bytes()).hexdigest(),metadata_group_counts=dict(Counter(r['group'] for r in manifest)),requested_position=a.position,candidate_sites=len(candidates),reference_alleles_are_not_ancestral_or_causal_labels=True,unquoted_taxon_is_provisional_identity_not_independent_verification=True,selection_and_calls_unchanged=True,comparison='Full optix 1.06-1.24Mb; raw calls exported at published Dennis intervals plus optional requested position',no_causality_or_introgression_claim=True)
 save(out/'complete.json',summary)
 if a.push_repo:
  repo=a.push_repo;branch='optix-meriana-audit-'+out.name
  env=dict(__import__('os').environ,GIT_SSH_COMMAND='ssh -o BatchMode=yes')
  def git(*args,cwd=repo):subprocess.run(['git','-C',str(cwd),*map(str,args)],env=env,check=True,stdout=subprocess.DEVNULL)
  with tempfile.TemporaryDirectory(prefix='optix-meriana-publish.') as temp:
   work=Path(temp)/'checkout';git('fetch','--quiet','origin','main');git('worktree','add','--quiet','--detach',work,'origin/main')
   try:
    dest=work/'results/optix_meriana_audit'/out.name;dest.mkdir(parents=True)
    for p in out.iterdir():
     if p.is_file() and p.suffix!='.gz':shutil.copy2(p,dest/p.name)
    git('add',str(dest.relative_to(work)),cwd=work);git('commit','--quiet','-m','Audit meriana optix identity and Dennis calls',cwd=work);git('push','--quiet','git@github.com:Hanliconius/GynandromorphCode.git','HEAD:refs/heads/'+branch,cwd=work)
   finally:git('worktree','remove','--force',work)
  print('GIT_BRANCH\t'+branch)
 print('SUMMARY\t'+json.dumps(summary));print('OUTPUT\t'+str(out))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--push-repo',type=Path);p.add_argument('--position',type=int);p.add_argument('--submit',action='store_true');a=p.parse_args()
 if a.submit:
  stage=ROOT/'analysis/optix_meriana_audit_code';stage.mkdir(exist_ok=True)
  frozen=stage/('audit_'+hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:12]+'.py');shutil.copy2(__file__,frozen)
  args=[sys.executable,str(frozen),'--output',str(a.output)]
  if a.push_repo:args+=['--push-repo',str(a.push_repo)]
  if a.position:args+=['--position',str(a.position)]
  cmd=['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=6G','--time=00:30:00','--job-name=Optix_meriana','--output='+str(stage/'audit.%j.out'),'--error='+str(stage/'audit.%j.err'),'--wrap='+shlex.join(args)]
  print('JOB\t'+subprocess.check_output(cmd,text=True).strip())
 else:do_audit(a)
if __name__=='__main__':main()
