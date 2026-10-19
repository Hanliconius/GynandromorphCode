#!/usr/bin/env python3
"""Audit meriana identity and optix genotype sharing without changing calls or selection."""
import argparse,csv,gzip,hashlib,json,re,shlex,shutil,subprocess,sys,tempfile
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
RUN=ROOT/'analysis/rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc'
SOURCES={'recalled':RUN/'integration_Optix/Optix_window.vcf.gz','legacy':ROOT/'external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf','merged':RUN/'Optix_legacy_merge_v1/Optix_legacy_plus_recalled.vcf.gz'}
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
  if not group:
   if 'rosina' in normal.lower():group='rosina'
   elif re.search(r"melpomene\s+(?:ssp\.?|spp\.?)?\s*melpomene$",normal,re.I):group='melpomene_melpomene'
   elif re.search(r'piano|butterfly world|hybrid|gynandromorph|backcross',normal,re.I):continue
   else:group='other_comparative'
  if not group:continue
  names=set()
  for col in ('Accession','SequenceID','EarthCapeID'):names.update(re.findall(r'(?:ERR|SRR|DRR)\d+',r.get(col,'')))
  for s in names:
   entry=dict(sample=s,group=group,raw_taxon=taxon,normal_taxon=normal,Country=r.get('Country',''),B=r.get('B',''),D=r.get('D',''),R=r.get('R',''),metadata_line=number,Latitude=r.get('Latitude',''),Longitude=r.get('Longitude',''),SequenceID=r.get('SequenceID',''),EarthCapeID=r.get('EarthCapeID',''))
   if s in individuals and individuals[s]!=entry:raise RuntimeError('Ambiguous duplicate sample metadata: '+s)
   individuals[s]=entry;manifest.append(entry)
 table(out/'sample_identities.tsv',manifest)
 calls=[];summaries=[];absent=[];comparison={};qc=Counter();race_sites=[];reference_checks=[]
 reference=ROOT/'analysis/optix_joint_all_samples/20261014T211057Z/ref/genome.fa'
 fasta=pysam.FastaFile(str(reference)) if reference.is_file() else None
 for label,path in SOURCES.items():
  if not path.is_file():absent.append(dict(source=label,path=str(path),reason='file_not_present'));continue
  with pysam.VariantFile(str(path)) as vf:
   samples=sorted(set(vf.header.samples)&set(individuals));counts=Counter(individuals[s]['group'] for s in samples)
   for s in set(individuals)-set(samples):absent.append(dict(source=label,sample=s,reason='sample_absent_from_VCF'))
   for rec in vf:
    if not 1214000<=rec.pos<=1222500 or len(rec.ref)!=1 or len(rec.alts or ())!=1 or len(rec.alts[0])!=1:continue
    interval='Dennis1' if 1215210<=rec.pos<=1216124 else 'Dennis2' if 1218557<=rec.pos<=1221494 else 'outside_published_D_intervals'
    if fasta is not None:
     base=fasta.fetch('NC_140406.1',rec.pos-1,rec.pos).upper()
     reference_checks.append(dict(source=label,pos=rec.pos,VCF_REF=rec.ref,FASTA_REF=base,matches=base==rec.ref))
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
    for taxon,values in races.items():
     vv=[v for v in values if v is not None]
     race_sites.append(dict(source=label,pos=rec.pos,ref=rec.ref,alt=rec.alts[0],interval=interval,taxon=taxon,n_total=len(values),n_called=len(vv),ALT_AF=sum(vv)/(2*len(vv)) if vv else None))
    phenotype={1:[],0:[]}
    for sample in samples:
     try:state=int(float(individuals[sample]['D']))
     except (ValueError,TypeError):continue
     if state in phenotype:
      v,_=dosage(rec.samples[sample])
      if v is not None:phenotype[state].append(v)
    stats={}
    for state,vals in phenotype.items():
     stats['D'+str(state)+'_n_called']=len(vals);stats['D'+str(state)+'_ALT_AF']=sum(vals)/(2*len(vals)) if vals else None
    for group in ('meriana_unquoted','meriana_quoted','besckei','timareta','rosina','melpomene_melpomene'):
     vals=[v for v in groups[group] if v is not None];stats[group+'_n_called']=len(vals);stats[group+'_alt_AF']=sum(vals)/(2*len(vals)) if vals else None
    true=[x for x in groups['meriana_unquoted'] if x is not None]
    same_alt=bool(true) and all(x>0 for x in true) and any(x is not None and x>0 for x in groups['besckei']) and any(x is not None and x>0 for x in groups['timareta'])
    summaries.append(dict(source=label,chrom=rec.contig,pos=rec.pos,ref=rec.ref,alt=rec.alts[0],interval=interval,all_callable_unquoted_meriana_share_ALT_with_besckei_and_timareta=same_alt,**stats))
 table(out/'Dennis_individual_raw_calls.tsv.gz',calls)
 grouped=defaultdict(Counter)
 for r in calls:
  key=(r['source'],r['group'],r['interval'])
  grouped[key]['records']+=1
  v=r['qualified_alt_dosage'];grouped[key]['missing' if v is None else 'REF_homozygote' if v==0 else 'heterozygote' if v==1 else 'ALT_homozygote']+=1
 table(out/'Dennis_group_genotype_summary.tsv',[dict(source=k[0],group=k[1],interval=k[2],**v) for k,v in sorted(grouped.items())])

 table(out/'race_site_ALT_frequencies.tsv.gz',race_sites)
 table(out/'reference_checks.tsv',reference_checks)
 if fasta is not None:fasta.close()
 if a.position:table(out/'requested_position_raw_calls.tsv',[r for r in calls if r['pos']==a.position])
 table(out/'Dennis_site_summary.tsv',[r for r in summaries if r['interval']!='outside_published_D_intervals'])
 candidates=[r for r in summaries if r['all_callable_unquoted_meriana_share_ALT_with_besckei_and_timareta']]
 table(out/'shared_ALT_candidate_sites.tsv',candidates)
 if a.position:table(out/'requested_position_summary.tsv',[r for r in summaries if r['pos']==a.position])
 table(out/'missing_sources_or_samples.tsv',absent)
 table(out/'QC_summary.tsv',[dict(source=k[0],metric=k[1],count=v) for k,v in sorted(qc.items())])
 summary=dict(status='AUDIT_COMPLETE_REQUIRES_IDENTITY_AND_SITE_REVIEW',metadata_sha256=hashlib.sha256(metadata.read_bytes()).hexdigest(),metadata_group_counts=dict(Counter(r['group'] for r in manifest)),requested_position=a.position,candidate_sites=len(candidates),reference_alleles_are_not_ancestral_or_causal_labels=True,unquoted_taxon_is_provisional_identity_not_independent_verification=True,selection_and_calls_unchanged=True,comparison='Dennis region 1.214-1.2225Mb; raw calls exported at published Dennis intervals plus optional requested position',no_causality_or_introgression_claim=True)

 # A labelled, explicit REF/ALT view avoids ambiguous allele orientation or
 # confusion between true meriana and quoted meriana taxon labels.
 import numpy as np,matplotlib
 matplotlib.use('Agg');matplotlib.rcParams.update({'pdf.fonttype':42})
 import matplotlib.pyplot as plt
 from matplotlib.colors import LinearSegmentedColormap
 from matplotlib.backends.backend_pdf import PdfPages
 cmap=LinearSegmentedColormap.from_list('ALT',['#D7EDF5','#ECCB63','#662A83']);cmap.set_bad('#C8C8C8')
 with PdfPages(out/'Optix_Dennis_meriana_REF_ALT_audit.pdf') as pdf:
  for source in SOURCES:
   rr=[r for r in race_sites if r['source']==source and r['interval']!='outside_published_D_intervals']
   if not rr:continue
   keys=sorted({(r['pos'],r['ref'],r['alt']) for r in rr});races=sorted({r['taxon'] for r in rr});mapping={(r['taxon'],(r['pos'],r['ref'],r['alt'])):r for r in rr}
   def page(mat,labels,title):
    fig,ax=plt.subplots(figsize=(11,max(3,len(labels)*.17+1.6)));fig.subplots_adjust(left=.40,right=.97,bottom=.18,top=.88)
    im=ax.imshow(np.ma.masked_invalid(mat),aspect='auto',interpolation='nearest',vmin=0,vmax=1,cmap=cmap)
    ax.set_yticks(range(len(labels)),labels,fontsize=6);ticks=np.linspace(0,len(keys)-1,min(8,len(keys))).astype(int);ax.set_xticks(ticks,[str(keys[k][0]) for k in ticks],fontsize=6);ax.set_xlabel('Dennis1 + Dennis2 SNPs; condensed columns, Rosina coordinates',fontsize=7)
    ax.set_title(source+' | '+title,fontsize=10);cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=.035,pad=.16);cb.set_label('ALT frequency / ALT dosage: 0=REF, 1=ALT; grey=missing',fontsize=7);cb.ax.tick_params(labelsize=6)
    pdf.savefig(fig);plt.close(fig)
   mat=np.array([[mapping.get((taxon,k),{}).get('ALT_AF',np.nan) if mapping.get((taxon,k),{}).get('ALT_AF') is not None else np.nan for k in keys] for taxon in races])
   page(mat,races,'race frequencies; raw VCF allele orientation')
   focus=sorted({r['sample'] for r in calls if r['source']==source and r['group'] in ('meriana_unquoted','meriana_quoted','rosina','melpomene_melpomene')},key=lambda x:(individuals[x]['group'],x))
   lookup={(r['sample'],(r['pos'],r['ref'],r['alt'])):r for r in calls if r['source']==source}
   mat=np.array([[lookup.get((sample,k),{}).get('qualified_alt_dosage')/2 if lookup.get((sample,k),{}).get('qualified_alt_dosage') is not None else np.nan for k in keys] for sample in focus])
   labels=[individuals[x]['group']+' | '+x+' | D='+individuals[x]['D'] for x in focus]
   if focus:page(mat,labels,'individual meriana / controls; quoted taxa kept separate')

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
     if p.is_file() and p.suffix not in ('.gz','.pdf'):shutil.copy2(p,dest/p.name)
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
