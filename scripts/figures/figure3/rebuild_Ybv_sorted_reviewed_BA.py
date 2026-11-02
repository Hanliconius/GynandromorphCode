#!/usr/bin/env python3
"""Reviewed Ybv regional BA with sorted population rows and no consensus rows."""
import argparse,hashlib,importlib.util,json,re,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
SOURCE=ROOT/'analysis/Figure3_for_Arnaud'
YBV=SOURCE/'Ybv_geography_supplements'
PRIOR=ROOT/'analysis/Ybv_all_populations_20261102'
FOCAL=ROOT/'analysis/supp_genoplots_LR_cydno_historical_BDR_20261027/tables/corrected_LR_calls.tsv'
OUTPUT=ROOT/'analysis/Ybv_sorted_reviewed_BA_20261102'
KEYS=['chrom','pos','ref','alt']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def load(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def prepare(a):
 if a.output.exists():
  if not a.retry:raise RuntimeError('Output exists; use --retry for this failed run: '+str(a.output))
  if (a.output/'complete.json').exists():raise RuntimeError('Refusing to replace a completed figure')
  job=str(json.loads((a.output/'jobs.json').read_text())['plot'])
  active=subprocess.run(['squeue','-h','-j',job,'-o','%T'],capture_output=True,text=True)
  if active.returncode==0 and active.stdout.strip():raise RuntimeError('Previous job is still active: '+job)
  account=subprocess.check_output(['sacct','-X','-n','-j',job,'--format=State','--parsable2'],text=True).strip()
  states=[x.strip().split('|')[0].split()[0].rstrip('+') for x in account.splitlines() if x.strip()]
  failed={'FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE'}
  if not states or any(x not in failed for x in states):raise RuntimeError('Cannot verify previous job failed: '+account)
  prior_cfg=json.loads((a.output/'run.json').read_text())
  for path,h in prior_cfg['hashes'].items():
   if sha(path)!=h:raise RuntimeError('Original run input changed; preserve failed snapshot: '+path)
  if sha(a.output/'reviewed_specimens.tsv')!=prior_cfg['metadata_sha256']:raise RuntimeError('Reviewed snapshot changed')

 repo=Path(__file__).resolve().parents[3]
 metadata=a.metadata.resolve() if a.metadata else repo/'data/figure3_review/Ybv_reviewed_metadata_20261102.tsv'
 helper=load(YBV/'figure3_helper.py','frozen_helper')
 prior=json.loads((PRIOR/'run.json').read_text())
 paths=[metadata,PRIOR/'worker.py',PRIOR/'complete.json',PRIOR/'run.json',FOCAL,
  YBV/'worker.py',YBV/'figure3_helper.py',YBV/'geography_helper.py',YBV/'run.json',
  YBV/'all_regional_BA_and_denominators.tsv.gz',YBV/'membership.tsv',
  YBV/'Ybv_geography_regional_union_sites.tsv',helper.BROAD/'SNPs.tsv',helper.YB_ENGINE,
  SOURCE/'data/Ybv/all_plot_focus_ALT_dosages.tsv.gz',SOURCE/'tables/Ybv_focus_markers_and_updated_BA.tsv']
 paths += [Path(x['path']) for x in prior['existing_VCFs']]
 for p in paths:
  if not p.is_file():raise FileNotFoundError('Missing documented input: '+str(p))
 if json.loads((PRIOR/'complete.json').read_text()).get('status')!='COMPLETE':raise RuntimeError('All-populations source incomplete')
 a.output.mkdir(parents=True,exist_ok=a.retry);(a.output/'logs').mkdir(exist_ok=a.retry);(a.output/'tables').mkdir(exist_ok=a.retry)
 shutil.copy2(__file__,a.output/'worker.py');shutil.copy2(metadata,a.output/'reviewed_specimens.tsv')
 save(a.output/'run.json',dict(hashes={str(p):sha(p) for p in paths},existing_VCFs=prior['existing_VCFs'],
  metadata_sha256=sha(metadata),metadata_source='Google Sheet 1YcJjcnE48q_aJ1wGBIYLOudx0euTJTn9IMWDLGRES20; Figure3ab_suppGenoPlots_allSamples; read 2026-11-02',
  BA_rule='Current P eligibility, H/I phenotypes, T geography; homozygote, >=80% per-class callability; West orientation held fixed across regions with labelled East fallback',
  SNP_selection_changed=False,consensus_rows_removed=True,new_calls=False))
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=8G','--cpus-per-task=1','--job-name=Ybv_reviewBA',
  '--output='+str(a.output/'logs/plot.%j.out'),'--error='+str(a.output/'logs/plot.%j.err'),
  '--wrap='+shlex.join([sys.executable,str(a.output/'worker.py'),'--run','--output',str(a.output)])],text=True).strip().split(';')[0]
 save(a.output/'jobs.json',dict(plot=job));print('JOB\t'+job+'\nOUTPUT\t'+str(a.output))

def existing_calls(out,cfg,sites,helper):
 import numpy as np,pandas as pd,pysam
 reviewed=helper.read_metadata(out/'reviewed_specimens.tsv')
 # Restore all reviewed comparative populations, without a taxon, phenotype
 # or geography restriction. Focal tissues appear only as corrected L/R.
 focal=reviewed.specimen_id.astype(str).str.contains(r'^(?:TN(?:_|[23][LR])|NDJ(?:_|$))',case=False,regex=True)
 focal |= reviewed.taxon.astype(str).str.contains('gynandromorph',case=False,regex=False)
 reviewed['display_exclusion']=np.where(reviewed.globally_excluded,'reviewed global exclusion',np.where(focal,'focal tissues displayed separately',''))
 reviewed.to_csv(out/'tables/all_reviewed_specimen_eligibility.tsv',sep='\t',index=False)
 meta=reviewed.loc[~reviewed.globally_excluded&~focal].copy()
 meta=meta.reset_index(drop=True)
 meta['geography_display']=meta.geography_used_Fig3b.map(lambda x:str(x).strip() or 'Unresolved')
 if meta.empty:raise RuntimeError('No reviewed target taxa')
 if meta.specimen_id.duplicated().any():raise RuntimeError('Duplicate reviewed specimen identity')
 target={tuple(r[c] for c in KEYS):j for j,r in enumerate(sites.to_dict('records'))}
 if len(target)!=len(sites):raise RuntimeError('Duplicate requested Yb allele identity')
 values=np.full((len(meta),len(sites)),np.nan);sources=np.full(values.shape,'no_existing_qualified_call',dtype=object)
 decided=np.zeros(values.shape,bool);quality=load(helper.YB_ENGINE,'existing_quality_policy')
 for entry in sorted(cfg['existing_VCFs'],key=lambda x:x['priority']):
  with pysam.VariantFile(entry['path']) as vf:
   names=set(vf.header.samples);sample_map={i:sorted(helper.aliases(r)&names) for i,r in enumerate(meta.to_dict('records'))}
   if not any(sample_map.values()):continue
   needed=sites.loc[sites.pos.between(entry['start'],entry['end'])]
   if needed.empty:continue
   # Reading a completed shard is cheap and works even if its index is absent.
   records=vf.fetch(helper.SPANS['Ybv'][0],int(needed.pos.min())-1,int(needed.pos.max())) if vf.index is not None else iter(vf)
   candidates={}
   for rec in records:
    if not entry['start']<=rec.pos<=entry['end']:continue
    for alt in rec.alts or ():
     k=(helper.canon(rec.contig),rec.pos,rec.ref,alt)
     if k not in target:continue
     j=target[k]
     for i,names_i in sample_map.items():
      if decided[i,j] or not names_i:continue
      one=[]
      for sample in names_i:
       c=rec.samples[sample];gt=c.get('GT') or ();valid,_=quality.quality(c,allow_pl=True)
       bases=[rec.alleles[a] for a in gt if a is not None and a<len(rec.alleles)]
       good=valid is not None and len(gt)==2 and len(bases)==2 and all(b in (rec.ref,alt) for b in bases)
       one.append(sum(b==alt for b in bases)/2 if good else np.nan)
      candidates.setdefault((i,j),[]).extend(one)
   for (i,j),calls in candidates.items():
    good=[x for x in calls if np.isfinite(x)]
    values[i,j]=good[0] if good and len(set(good))==1 else np.nan
    sources[i,j]=entry['path'] if len(set(good))<=1 else 'conflicting_existing_calls_missing'
    decided[i,j]=True
 # Exact original-marker calls already used in earlier plots take precedence.
 cache=SOURCE/'data/Ybv/all_plot_focus_ALT_dosages.tsv.gz'
 if cache.is_file():
  focus=pd.read_csv(SOURCE/'tables/Ybv_focus_markers_and_updated_BA.tsv',sep='\t');f=pd.read_csv(cache,sep='\t',index_col=0)
  if list(map(int,f.columns))!=focus.pos.tolist():raise RuntimeError('Frozen Yb context cache order differs')
  for r in focus.itertuples():
   k=(helper.canon(r.chrom),int(r.pos),r.ref,r.alt)
   if k not in target:continue
   j=target[k]
   for i,rmeta in meta.iterrows():
    if rmeta.specimen_id in f.index:values[i,j]=float(f.loc[rmeta.specimen_id,str(r.pos)]);sources[i,j]='unchanged_existing_focus_cache'
 meta['display_taxon']=meta.population.where(meta.population.str.strip().ne(''),meta.taxon).map(lambda s:' '.join(re.sub(r'\b(?:ssp|spp|subsp)\.?(?=\s|$)','',str(s),flags=re.I).split()))
 meta['morph']=['d'+str(r.Yb_dorsal or '?')+' v'+str(r.Yb_ventral or '?') for r in meta.itertuples()]
 meta['qualified_sites']=np.isfinite(values).sum(axis=1);meta['context_only']=True
 original=pd.read_csv(YBV/'membership.tsv',sep='\t',dtype=str,keep_default_na=False)
 meta['already_in_frozen_display']=meta.specimen_id.isin(original.specimen_id)
 meta.to_csv(out/'tables/all_populations_membership.tsv',sep='\t',index=False)
 pd.DataFrame(values,index=meta.specimen_id,columns=sites.pos).to_csv(out/'tables/all_populations_ALT_dosages.tsv.gz',sep='\t')
 # Aggregate call-source accounting avoids millions of redundant text rows.
 source_audit=[]
 for i,r in meta.iterrows():
  labels,counts=np.unique(sources[i],return_counts=True)
  source_audit.extend(dict(specimen_id=r.specimen_id,source=str(label),sites=int(n)) for label,n in zip(labels,counts))
 pd.DataFrame(source_audit).to_csv(out/'tables/call_source_counts.tsv',sep='\t',index=False)
 if not np.isfinite(values).any():raise RuntimeError('Target taxa found but no existing qualified exact-allele calls; no new calling attempted')
 return meta,values



def sorted_rows(meta,values,sites,helper):
 import numpy as np,pandas as pd
 usable=np.isfinite(values).any(axis=1);m=meta.loc[usable].reset_index(drop=True).copy();v=values[usable]
 groups=list(m.groupby(['geography_display','Yb_ventral','Yb_dorsal','display_taxon'],sort=False,dropna=False))
 region_order={'East':0,'West':1,'Silvaniforms':2,'Unresolved':3};state_order={'1':0,'0':1,'?':2,'':2}
 groups.sort(key=lambda x:(region_order.get(x[0][0],3),x[0][0],state_order.get(x[0][1],2),state_order.get(x[0][2],2),x[0][3]))
 rows=[];labels=[];exports=[];blocks=[];last=None
 for (region,state,dorsal,pop),g in groups:
  if (region,state)!=last:blocks.append((len(rows),len(rows)+1,region,state));last=(region,state)
  x=v[g.index];n=np.isfinite(x).sum(axis=0);freq=np.divide(np.nansum(x,axis=0),n,out=np.full(len(sites),np.nan),where=n>0)
  morph='d'+(dorsal if dorsal in ['0','1'] else '?')+' v'+(state if state in ['0','1'] else '?')
  labels.append(region+' | '+helper.short(pop)+' | '+morph+f' | n={len(g)}');rows.append(freq)
  for j,r in enumerate(sites.to_dict('records')):exports.append(dict(row=len(rows)-1,population=pop,geography=region,morph=morph,**{k:r[k] for k in KEYS},n_total=len(g),n_called=int(n[j]),ALT_frequency=float(freq[j]) if n[j] else None,specimen_ids=';'.join(g.specimen_id)))
 return np.asarray(rows).reshape((-1,len(sites))),labels,pd.DataFrame(exports),blocks

def reconcile_inventory(broad,selected,helper,out):
 """Union exact saved alleles with the full scan; never swap REF/ALT."""
 import pandas as pd
 def normal(table):
  t=table.copy()
  t['chrom']=t.chrom.map(helper.canon)
  t['pos']=pd.to_numeric(t.pos,errors='raise').astype(int)
  for c in ['ref','alt']:t[c]=t[c].astype(str).str.strip()
  return t
 broad,selected=normal(broad),normal(selected)
 if selected.duplicated(KEYS).any():
  selected.loc[selected.duplicated(KEYS,keep=False)].to_csv(out/'tables/duplicate_selected_exact_alleles.tsv',sep='\t',index=False)
  raise RuntimeError('Saved display contains duplicate exact alleles; see duplicate_selected_exact_alleles.tsv')
 if broad.empty or broad.duplicated(KEYS).any():raise RuntimeError('Empty/duplicate broad inventory')
 known={tuple(r) for r in broad[KEYS].itertuples(index=False,name=None)}
 audit=[]
 for row in selected[KEYS].to_dict('records'):
  key=tuple(row[c] for c in KEYS)
  alternatives=broad.loc[broad.chrom.eq(row['chrom'])&broad.pos.eq(row['pos'])]
  audit.append(dict(**row,in_broad_inventory=key in known,
   action='reuse_broad_exact_allele' if key in known else 'add_saved_display_allele',
   broad_alleles_at_position=';'.join(alternatives.ref+'/'+alternatives.alt)))
 pd.DataFrame(audit).to_csv(out/'tables/selected_allele_inventory_audit.tsv',sep='\t',index=False)
 # This is an inventory reconciliation, not new variant discovery. All entries
 # were already in completed calls or in the saved displayed SNP table.
 sites=pd.concat([broad[KEYS],selected[KEYS]],ignore_index=True).drop_duplicates(KEYS).sort_values(KEYS).reset_index(drop=True)
 sites.to_csv(out/'tables/reconciled_exact_allele_inventory.tsv',sep='\t',index=False)
 return sites,selected

def run(out):
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 from matplotlib.backends.backend_pdf import PdfPages
 cfg=json.loads((out/'run.json').read_text())
 for p,h in cfg['hashes'].items():
  if sha(p)!=h:raise RuntimeError('Input changed after submission: '+p)
 if sha(out/'reviewed_specimens.tsv')!=cfg['metadata_sha256']:raise RuntimeError('Metadata snapshot changed')
 helper=load(YBV/'figure3_helper.py','frozen_helper');yb=load(YBV/'worker.py','frozen_renderer');style=load(PRIOR/'worker.py','previous_style');geo=load(YBV/'geography_helper.py','unchanged_BA_rule')
 broad=pd.read_csv(helper.BROAD/'SNPs.tsv',sep='\t');broad['chrom']=broad.chrom.map(helper.canon)
 broad=broad.loc[broad.chrom.eq(helper.SPANS['Ybv'][0])&broad.pos.between(*helper.SPANS['Ybv'][1:])].reset_index(drop=True)
 selected=pd.read_csv(YBV/'Ybv_geography_regional_union_sites.tsv',sep='\t')
 sites,selected=reconcile_inventory(broad,selected,helper,out)
 # Reuse all existing qualified calls; this does not realign, recall variants or
 # refit focal mixture estimates. Updated reviewed metadata supplies grouping.
 meta,values=existing_calls(out,cfg,sites,helper)
 meta['sample']=meta.specimen_id;meta['phenotype']=pd.to_numeric(meta.Yb_ventral,errors='coerce')
 eligible=meta.BA_Fig3b.map(helper.yes)&meta.Yb_ventral.isin(['0','1'])&meta.geography_used_Fig3b.isin(['East','West'])
 why=np.where(~meta.BA_Fig3b.map(helper.yes),'column P unchecked',np.where(~meta.Yb_ventral.isin(['0','1']),'unresolved ventral phenotype',np.where(~meta.geography_used_Fig3b.isin(['East','West']),'outside East/West or unresolved geography','')))
 meta['used_in_BA']=eligible;meta['BA_exclusion_reason']=why
 meta['all_sites_missing']=~np.isfinite(values).any(axis=1)
 meta.to_csv(out/'tables/specimen_accounting.tsv',sep='\t',index=False)
 members=meta.loc[eligible].copy()
 if members.empty:raise RuntimeError('No eligible BA specimens')
 geo.ensure_allowed(members[['sample','taxon','population']])
 members['phenotype']=members.Yb_ventral.astype(int)
 regional=yb.regional_scores(sites,values[eligible],members,geo)
 regional['pos']=regional.pos.astype(int)
 lookup={helper.key(r.chrom,r.pos,r.ref,r.alt):j for j,r in enumerate(sites.itertuples())}
 keys=[helper.key(r.chrom,r.pos,r.ref,r.alt) for r in selected.itertuples()]
 if len(keys)!=len(set(keys)) or any(k not in lookup for k in keys):raise RuntimeError('Reconciled exact-allele mapping failed; see selected_allele_inventory_audit.tsv')
 indices=[lookup[k] for k in keys]
 regional['original_marker']=[k in set(helper.key(r.chrom,r.pos,r.ref,r.alt) for r in selected.loc[selected.original_marker].itertuples()) for k in (helper.key(r.chrom,r.pos,r.ref,r.alt) for r in regional.itertuples())]
 # Keep the saved SNP columns; colour the regional BA tracks by exact membership.
 regional['displayed_in_genoplot']=[helper.key(r.chrom,r.pos,r.ref,r.alt) in set(keys) for r in regional.itertuples()]
 regional.to_csv(out/'tables/recalculated_regional_BA_all_SNPs.tsv.gz',sep='\t',index=False)
 selected_new=regional.iloc[indices].reset_index(drop=True)
 if selected[KEYS].to_dict('records')!=selected_new[KEYS].to_dict('records'):raise RuntimeError('Genoplot SNP order changed')
 # Supply the current BA values to the fixed-site display, never select again.
 for c in selected.columns:
  if c not in selected_new:selected_new[c]=selected[c].to_numpy()
 selected_new.to_csv(out/'tables/fixed_selected_SNPs_updated_BA.tsv',sep='\t',index=False)
 counts=members.groupby(['geography_used_Fig3b','Yb_ventral']).size().reset_index(name='eligible_specimens')
 counts.to_csv(out/'tables/BA_cohort_counts.tsv',sep='\t',index=False)
 old=pd.read_csv(YBV/'all_regional_BA_and_denominators.tsv.gz',sep='\t');old['chrom']=old.chrom.map(helper.canon)
 cols=['West_best_BA','East_best_BA','West_fixed_direction_BA','East_fixed_direction_BA']
 shifts=regional[KEYS+cols].merge(old[KEYS+cols],on=KEYS,how='left',suffixes=('_new','_old'),validate='one_to_one')
 for c in cols:shifts[c+'_change']=shifts[c+'_new']-shifts[c+'_old']
 shifts.to_csv(out/'tables/BA_changes.tsv.gz',sep='\t',index=False)
 rows,labels,exports,blocks=sorted_rows(meta,values[:,indices],selected_new,helper)
 if not len(rows):raise RuntimeError('No populations have qualified selected-site calls')
 focal=pd.read_csv(FOCAL,sep='\t');focal=focal.loc[focal.locus.eq('Ybv')].copy();focal['chrom']=focal.chrom.map(helper.canon)
 for _,g in focal.groupby(KEYS):
  if len(g[['Left_ALT_dosage','Right_ALT_dosage']].drop_duplicates())!=1:raise RuntimeError('Conflicting corrected focal values')
 f=selected_new[KEYS].merge(focal.drop_duplicates(KEYS),on=KEYS,how='left',validate='one_to_one',indicator=True)
 if not f['_merge'].eq('both').all():raise RuntimeError('Missing exact saved focal allele')
 lr=f[['Left_ALT_dosage','Right_ALT_dosage']].to_numpy(float).T
 saved=np.vstack([rows,lr]),labels+['Left (corrected)','Right (corrected)'],exports,blocks
 yb.grouped_rows=lambda *args:saved
 original=PdfPages.savefig
 def finish(self,figure=None,**kwargs):
  style.remove_original_marker_squares(figure);style.add_phenotype_table(figure,'Ybv')
  for t in figure.texts:
   if 'West v+, West v−' in t.get_text():
    t.set_text('Rows: reviewed column T, then v+ / v− / unknown; all populations; no pooled consensus rows.\nBA: current column P eligibility, homozygote rule, >=80% per-class callability; West-best direction fixed across regions, with labelled East fallback.\nColumns: unchanged previous regional-union SNP set; BA recalculated, SNPs not reselected. Blue points identify displayed sites.\nALT frequency: 0=REF, 1=ALT; grey=missing. n=specimens; per-site denominators exported. Corrected L/R unchanged.')
  return original(self,figure,**kwargs)
 PdfPages.savefig=finish
 pdf=out/'Ybv_all_populations_geography_phenotype_updated_BA.pdf'
 try:_,pages=yb.draw(pdf,selected_new,regional,None,None,helper,'Yb(v): populations by geography and phenotype',json.loads((YBV/'run.json').read_text())['cutoff'])
 finally:PdfPages.savefig=original
 exports.to_csv(out/'tables/population_frequencies_and_membership.tsv.gz',sep='\t',index=False)
 # Verify genotype preservation by exact site identity; all grouping derives
 # from the current reviewed metadata, never from the observed genotype pattern.
 save(out/'complete.json',dict(status='COMPLETE',PDF=str(pdf),pages=pages,SNPs=len(selected_new),population_rows=len(rows),displayed_specimens=int(np.isfinite(values[:,indices]).any(axis=1).sum()),BA_eligible_specimens=len(members),BA_eligible_with_no_calls=int(members.all_sites_missing.sum()),consensus_rows=0,BA_recalculated=True,SNP_selection_changed=False,corrected_focal_unchanged=True,metadata_sha256=cfg['metadata_sha256']))
 (out/'README.txt').write_text('Current reviewed Google Sheet snapshot supplies geography, phenotype codes and BA eligibility. Population display uses all eligible comparative specimens with qualified calls, irrespective of P; BA uses only P=true, known v phenotype and East/West geography. Global exclusions and focal separation remain enforced. Silvaniforms and Unresolved retain separate display blocks and do not enter East/West BA. Homozygote BA uses >=80% callability separately per phenotype class; both allele orientations are exported, West-best orientation is held fixed across regional tracks, with labelled East fallback. The full original SNP inventory plus exact saved display alleles is rescored; added inventory entries are audited without changing REF/ALT or genoplot columns; the regional-union genoplot SNP list is NOT reselected. Four pooled consensus rows removed. Corrected focal values copied exactly from the prior completed output. No downloads, alignment, variant calling or mixture refitting. Eligible specimens with no calls remain explicit in denominators and accounting.\n')
 print('PDF\t'+str(pdf)+'\nSTATUS\tCOMPLETE')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--submit',action='store_true');p.add_argument('--retry',action='store_true');p.add_argument('--run',action='store_true');p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--metadata',type=Path,help='Explicit reviewed TSV snapshot; otherwise use original frozen snapshot');a=p.parse_args()
 if a.submit:prepare(a)
 elif a.run:run(a.output)
 else:p.error('Use --submit')
if __name__=='__main__':main()
