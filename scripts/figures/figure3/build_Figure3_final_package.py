#!/usr/bin/env python3
"""Rebuild Figure 3 and supplements from the author-reviewed specimen snapshot.

Uses existing calls only. One specimen contributes once; concordant qualified
run calls are united, conflicts missing. Main focus marker identities stay fixed.
"""
import argparse,csv,getpass,hashlib,importlib.util,json,re,shlex,shutil,subprocess,sys,zipfile
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
RUN=ROOT/'analysis/rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc'
PARENT=ROOT/'analysis/Optix_Ybv_integrated_figure'
DEN=RUN/'Optix_Arnaud_review_nonPK_v2/merged'
DF=ROOT/'analysis/Dennis_figure_set_v1'
YF=ROOT/'analysis/Ybv_figure_set_v1'
HO=ROOT/'analysis/LR_window_heterozygosity_v1'
OLD=ROOT/'analysis/yb_completed_reads_recovery_v1/frozen/HmRos_Ybv_unchanged_focal_calls.tsv.gz'
VCF=RUN/'Optix_legacy_merge_v1/Optix_legacy_plus_recalled.vcf.gz'
BROAD=PARENT/'Ybv_BA_broad'
YB_FULL=ROOT/'analysis/ybv_REF_ALT_rebuild_v1/full_window.vcf.gz'
YB_ENGINE=ROOT/'analysis/yb_completed_reads_recovery_v1/engine.py'
GFF=ROOT/'analysis/genomewide_mixing_input/reference/GCF_056956725.1/GCF_056956725.1_ASM5695672v1_genomic.gff.gz'
SPANS={'Dennis':('NC_140406.1',1060000,1240000),'Ybv':('NC_140403.1',1630000,1692000)}
DENY={'ERR1303880','ERR1303882','ERR1303884','ERR1303885','ERR260303','ERR260304'}
BLUE='#0072D5';YELLOW='#FFD000'
SUPP=['#D7EDF5','#ECCB63','#662A83']
TRAITS=['B','D','R','Yb_dorsal','Yb_ventral','N','Sb','Br']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def load(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def yes(v):return str(v).lower() in ('true','1','1.0')
def normal(s):return ' '.join(str(s).split())
def aliases(r):return set(re.findall(r'\b(?:ERR|SRR|DRR|SAMN|SAMEA)\d+\b',' '.join(str(r.get(k,'')) for k in ['specimen_id','run_accessions','BioSample_accessions'])))
def blocked(r):return yes(r.get('exclude_from_all_calculations','')) or bool(re.search('numata|pardalinus',str(r.get('taxon',''))+' '+str(r.get('population','')),re.I)) or bool(aliases(r)&DENY)
def canon(chrom):return {'Hmr18':'NC_140406.1','18':'NC_140406.1','Hmr15':'NC_140403.1','15':'NC_140403.1'}.get(str(chrom),str(chrom))
def key(chrom,pos,ref,alt):return (canon(chrom),int(pos),str(ref),str(alt))
def read_metadata(p):
 import pandas as pd
 m=pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False).fillna('')
 if m.specimen_id.duplicated().any():raise RuntimeError('Duplicate stable specimen ID')
 for c in ['BA_Fig3a','BA_Fig3b','exclude_from_all_calculations']:
  if not m[c].str.lower().isin(['true','false']).all():raise RuntimeError('Invalid Boolean in '+c)
 for c in TRAITS:
  if not m[c].isin(['0','1','?','']).all():raise RuntimeError('Invalid phenotype code in '+c)
 # Cythera spelling variants denote the same population, not distinct cohorts.
 cythera=m.population.str.fullmatch(r'Heliconius\s+melpomene\s+(?:(?:ssp|spp)\.?\s+)?cythera',case=False)
 m.loc[cythera,'population']='Heliconius melpomene cythera'
 m['globally_excluded']=[blocked(r) for r in m.to_dict('records')]
 # This snapshot is the input authority. No genotype-dependent exclusions.
 return m
def selected(m,locus):
 flag,trait=('BA_Fig3a','D') if locus=='Dennis' else ('BA_Fig3b','Yb_ventral')
 s=m.loc[m[flag].map(yes)&~m.globally_excluded].copy()
 if not len(s) or not s[trait].isin(['0','1']).all():raise RuntimeError('Selected cohort has unresolved phenotype: '+locus)
 s['phenotype']=s[trait].astype(int)
 if locus=='Ybv' and not s.geography_used_Fig3b.isin(['West','East']).all():raise RuntimeError('Yb cohort has unresolved geography')
 return s
def union(values):
 import numpy as np
 x=np.asarray(values,float);called=np.isfinite(x);n=called.sum(axis=0)
 lo=np.min(np.where(called,x,np.inf),axis=0);hi=np.max(np.where(called,x,-np.inf),axis=0)
 conflict=(n>0)&(lo!=hi)
 return np.where((n>0)&~conflict,lo,np.nan),int(conflict.sum())
def map_sources(meta,source_records):
 lookup={}
 for i,r in enumerate(meta.to_dict('records')):
  for a in aliases(r):
   if a in lookup and lookup[a]!=i:raise RuntimeError('Run/BioSample links multiple specimens: '+a)
   lookup[a]=i
 mapping={};used={}
 for sample,ids in source_records:
  hit={lookup[a] for a in ids if a in lookup}
  if len(hit)>1:raise RuntimeError('Source aliases link different reviewed specimens: '+sample)
  if hit:
   i=next(iter(hit));mapping.setdefault(i,[]).append(sample);used[sample]=i
 missing=[meta.iloc[i].specimen_id for i in range(len(meta)) if i not in mapping]
 if missing:raise RuntimeError('Selected specimen lacks existing calls: '+', '.join(missing))
 return mapping
def prepare(a):
 meta=a.metadata;out=a.output;repo=Path(__file__).resolve().parents[3]
 restriction=out/'Ybv_BA_rule_comparison/lower_panel'
 restricted=getattr(a,'yb_cohort','west')=='lower-panel'
 paths=[meta,VCF,DF/'frozen_markers_and_BA.tsv',DEN/'Dennis_BA90_intervals_frequencies.tsv',
 BROAD/'complete.json',BROAD/'membership.tsv',BROAD/'SNPs.tsv',BROAD/'individual_ALT_dosages.tsv.gz',
 YF/'individual_ALT_dosages.tsv.gz',OLD,PARENT/'worker.py',PARENT/'Ybv_gene_features.json',HO/'worker.py',
 HO/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz',RUN/'optix_consistent_v2/reciprocal_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz',
  repo/'scripts/figures/figure3/compare_Ybv_geography_BA.py',YB_FULL,YB_ENGINE,GFF]
 if restricted:paths += [restriction/n for n in ['complete.json','membership.tsv','selected_population_rows.tsv','Ybv_27_SNP_BA_comparison.tsv']]
 for p in paths:
  if not p.is_file():raise RuntimeError('Missing documented input: '+str(p))
 if json.loads((BROAD/'complete.json').read_text())['status']!='COMPLETE':raise RuntimeError('Yb calls incomplete')
 m=read_metadata(meta);selected(m,'Dennis');selected(m,'Ybv')
 if restricted:
  done=json.loads((restriction/'complete.json').read_text())
  if done.get('status')!='COMPLETE' or done.get('cohort')!='lower-panel':raise RuntimeError('Completed exact lower-panel diagnostic required')
 if out.exists():
  if not a.overwrite:raise RuntimeError('Use --overwrite to replace this package in place')
  prior=json.loads((out/'run.json').read_text())
  if prior.get('workflow')!='Figure3_reviewed_package':raise RuntimeError('Output owned by another workflow')
  jobs=json.loads((out/'jobs.json').read_text()) if (out/'jobs.json').exists() else {}
  active=set(subprocess.check_output(['squeue','-h','-u',getpass.getuser(),'-o','%A'],text=True).split())
  if active&set(map(str,jobs.values())):raise RuntimeError('Package jobs still active')
 out.mkdir(parents=True,exist_ok=True)
 for d in ['logs','code','data','tables']:(out/d).mkdir(exist_ok=True)
 # Only generated files belonging to this workflow are replaced.
 for p in [out/'complete.json',out/'Figure3_for_Arnaud.zip',*out.glob('*.pdf')]:p.unlink(missing_ok=True)
 for d in ['Dennis','Ybv']:(out/'data'/d).mkdir(exist_ok=True);(out/'data'/d/'complete.json').unlink(missing_ok=True)
 if restricted:
  for name in ['membership.tsv','selected_population_rows.tsv','Ybv_27_SNP_BA_comparison.tsv']:shutil.copy2(restriction/name,out/'code'/('restricted_'+name))
 shutil.copy2(meta,out/'reviewed_specimens.tsv');shutil.copy2(__file__,out/'code/worker.py')
 for src,dst in [(PARENT/'worker.py','panel_helper.py'),(HO/'worker.py','Ho_helper.py'),
  (PARENT/'Ybv_gene_features.json','Ybv_gene_features.json'),(repo/'scripts/figures/figure3/compare_Ybv_geography_BA.py','geography_helper.py'),(YB_ENGINE,'quality_engine.py')]:shutil.copy2(src,out/'code'/dst)
 save(out/'run.json',dict(workflow='Figure3_reviewed_package',metadata_sha256=sha(meta),
  metadata_source='Author-reviewed Google Sheet snapshot, 2026-10-22',metadata_url='https://docs.google.com/spreadsheets/d/1YcJjcnE48q_aJ1wGBIYLOudx0euTJTn9IMWDLGRES20/edit',
  input_hashes={str(p):sha(p) for p in paths},code_sha256=sha(__file__),no_new_calls=True,
  spans=SPANS,focus_markers_retained={'Dennis':24,'Ybv':27},Ybv_cohort='lower-panel' if restricted else 'west',
  Ybv_BA='Exact frozen lower-panel specimens; homozygote rule; allele direction selected within this cohort' if restricted else 'West only; homozygote rule; final pooled orientation held fixed as in geography diagnostic',
  Ybv_lower_panel='Same exact specimens as main BA; eight recorded rows including elevatus and merged cythera' if restricted else 'All final eligible comparison specimens; geography/race examples with elevatus replacing linaresi',
  population_aliases={'Heliconius melpomene ssp. cythera':'Heliconius melpomene cythera'},
  Dennis_BA='Carrier rule; final reviewed cohort',duplicate_runs='Concordant union per specimen; conflicts missing',
  main_palette={'allele_1':YELLOW,'allele_2':BLUE},supplement_palette=SUPP))
 jobs={}
 for stage,array,deps in [('stats','0-1',None),('render',None,'stats')]:
  args=['sbatch','--parsable','--partition=nano','--time=00:30:00','--cpus-per-task=1','--mem=8G',
   '--job-name=Fig3_'+stage,'--output='+str(out/'logs'/(stage+'.%A_%a.out')),'--error='+str(out/'logs'/(stage+'.%A_%a.err'))]
  if array:args+=['--array='+array]
  if deps:args+=['--dependency=afterok:'+jobs[deps],'--kill-on-invalid-dep=yes']
  args+=['--wrap='+shlex.join([sys.executable,str(out/'code/worker.py'),'--stage',stage,'--output',str(out)])]
  jobs[stage]=subprocess.check_output(args,text=True).strip().split(';')[0];save(out/'jobs.json',jobs)
  print(stage.upper()+'\t'+jobs[stage])
 print('OUTPUT\t'+str(out)+'\nMAIN_PDF\t'+str(out/'01_Figure3_main_BA_dots.pdf'))
def restrict_membership(meta,path):
 import pandas as pd
 frozen=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False)
 if frozen.specimen_id.duplicated().any() or not set(frozen.specimen_id)<=set(meta.specimen_id):raise RuntimeError('Restricted specimen identity mismatch')
 selected=meta.set_index('specimen_id').loc[frozen.specimen_id].reset_index()
 for c in ['phenotype','population','geography_used_Fig3b']:
  a=selected[c].astype(str).to_numpy();b=frozen[c].astype(str).to_numpy()
  if not (a==b).all():raise RuntimeError('Reviewed metadata differs from completed restricted diagnostic: '+c)
 if not len(selected) or set(selected.phenotype)!={0,1}:raise RuntimeError('Both restricted phenotype classes required')
 return selected

def use_focus_calls(sites,matrix,focus,narrow,out):
 """Use the same qualified focal-marker calls as the completed cohort diagnostic."""
 import numpy as np,pandas as pd
 lookup={key(r.chrom,r.pos,r.ref,r.alt):i for i,r in enumerate(sites.itertuples())};audit=[]
 for j,r in enumerate(focus.itertuples()):
  k=key(r.chrom,r.pos,r.ref,r.alt);i=lookup.get(k)
  if i is None:
   lookup[k]=len(sites);sites=pd.concat([sites,focus.iloc[[j]][['chrom','pos','ref','alt']]],ignore_index=True)
   matrix=np.column_stack([matrix,narrow[:,j]]);status='retained_marker_added_from_existing_focus_cache';missing_diff=0
  else:
   a=matrix[:,i];b=narrow[:,j];both=np.isfinite(a)&np.isfinite(b)
   if not np.allclose(a[both],b[both],atol=1e-12,rtol=0):raise RuntimeError('Broad/focus genotype disagreement at '+str(r.pos))
   missing_diff=int((np.isfinite(a)!=np.isfinite(b)).sum());matrix[:,i]=b;status='concordant_calls_use_exact_focus_cache'
  audit.append(dict(pos=r.pos,ref=r.ref,alt=r.alt,status=status,different_missingness=missing_diff))
 pd.DataFrame(audit).to_csv(out/'tables/Ybv_focus_call_source_alignment.tsv',sep='\t',index=False)
 return sites,matrix

def prepare_render_only(a):
 out=a.output
 if not a.overwrite:raise RuntimeError('Use --overwrite for an in-place redraw')
 cfg=json.loads((out/'run.json').read_text())
 if cfg.get('workflow')!='Figure3_reviewed_package' or cfg.get('Ybv_cohort')!='lower-panel':raise RuntimeError('Completed restricted Figure3 package required')
 for locus in ['Dennis','Ybv']:
  done=json.loads((out/'data'/locus/'complete.json').read_text())
  if done.get('status')!='COMPLETE':raise RuntimeError('Statistics incomplete: '+locus)
 if sha(out/'reviewed_specimens.tsv')!=cfg['metadata_sha256']:raise RuntimeError('Frozen metadata changed')
 jobs=json.loads((out/'jobs.json').read_text())
 active=set(subprocess.check_output(['squeue','-h','-u',getpass.getuser(),'-o','%A'],text=True).split())
 if active&set(map(str,jobs.values())):raise RuntimeError('Package jobs still active')
 shutil.copy2(__file__,out/'code/worker.py')
 cfg['render_code_sha256']=sha(__file__);cfg['render_only']=True
 save(out/'run.json',cfg)
 (out/'complete.json').unlink(missing_ok=True)
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--cpus-per-task=1','--mem=8G',
  '--job-name=Fig3_render','--output='+str(out/'logs/render.%A_%a.out'),'--error='+str(out/'logs/render.%A_%a.err'),
  '--wrap='+shlex.join([sys.executable,str(out/'code/worker.py'),'--stage','render','--output',str(out)])],text=True).strip().split(';')[0]
 jobs['render']=job;save(out/'jobs.json',jobs)
 print('RENDER\t'+job+'\nMAIN_PDF\t'+str(out/'01_Figure3_main_BA_dots.pdf'))

def statistics(out,locus):
 import numpy as np,pandas as pd
 folder=out/'data'/locus;cfg=json.loads((out/'run.json').read_text())
 if sha(out/'reviewed_specimens.tsv')!=cfg['metadata_sha256']:raise RuntimeError('Metadata snapshot changed')
 meta=selected(read_metadata(out/'reviewed_specimens.tsv'),locus).reset_index(drop=True)
 restricted=locus=='Ybv' and cfg.get('Ybv_cohort')=='lower-panel'
 if restricted:meta=restrict_membership(meta,out/'code/restricted_membership.tsv')
 for p,h in cfg['input_hashes'].items():
  if sha(p)!=h:raise RuntimeError('Input changed since submission: '+p)
 geo=load(out/'code/geography_helper.py','geo');conflicts=[]
 if locus=='Dennis':
  import pysam
  focus=pd.read_csv(DF/'frozen_markers_and_BA.tsv',sep='\t');focus['chrom']=focus.chrom.map(canon)
  if len(focus)!=24:raise RuntimeError('Dennis focus count differs')
  with pysam.VariantFile(str(VCF)) as vf:
   names=list(vf.header.samples);mapping=map_sources(meta,[(s,{s}) for s in names]);needed=sorted({s for ss in mapping.values() for s in ss});indices={s:i for i,s in enumerate(needed)}
   recs=[];vectors=[];seen=set();chrom,lo,hi=SPANS[locus]
   for rec in vf:
    if canon(rec.contig)!=chrom or not lo<=rec.pos<=hi or len(rec.ref)!=1 or len(rec.alts or ())!=1 or len(rec.alts[0])!=1 or rec.ref not in 'ACGT' or rec.alts[0] not in 'ACGT':continue
    kk=key(rec.contig,rec.pos,rec.ref,rec.alts[0])
    if kk in seen:raise RuntimeError('Duplicate Dennis SNP identity')
    seen.add(kk);vals=[]
    for s in needed:
     c=rec.samples[s];gt=c.get('GT') or ();dp=c.get('DP');gq=c.get('GQ')
     good=len(gt)==2 and all(v in (0,1) for v in gt) and dp is not None and dp>=2 and gq is not None and gq>=20
     vals.append(sum(gt)/2 if good else np.nan)
    recs.append(dict(zip(['chrom','pos','ref','alt'],kk)));vectors.append(vals)
  if not recs:raise RuntimeError('No Dennis SNP records in documented VCF')
  sites=pd.DataFrame(recs);raw=np.asarray(vectors).T;rows=[]
  for i in range(len(meta)):
   vals,n=union(raw[[indices[s] for s in mapping[i]]]);rows.append(vals)
   conflicts.append(dict(specimen_id=meta.iloc[i].specimen_id,locus=locus,source_samples='; '.join(mapping[i]),conflicting_sites=n))
  matrix=np.asarray(rows)
  lookup={key(r.chrom,r.pos,r.ref,r.alt):j for j,r in enumerate(sites.itertuples())}
  try:fi=[lookup[key(r.chrom,r.pos,r.ref,r.alt)] for r in focus.itertuples()]
  except KeyError as e:raise RuntimeError('Frozen Dennis focus marker missing: '+str(e))
  narrow=matrix[:,fi];oldfreq=pd.read_csv(DEN/'Dennis_BA90_intervals_frequencies.tsv',sep='\t',index_col=0)
  kk={key(*col.split(':')):col for col in oldfreq.columns}
  if oldfreq.shape!=(4,24) or len(kk)!=24:raise RuntimeError('Frozen corrected Dennis table shape differs')
  columns=[kk[key(r.chrom,r.pos,r.ref,r.alt)] for r in focus.itertuples()]
  focal=oldfreq.loc[:,columns].iloc[[1,2]].to_numpy(float)
  orient=np.asarray([yes(x) for x in focus.pattern_is_alt]);focal=np.where(orient,focal,1-focal)
 else:
  sm=pd.read_csv(BROAD/'membership.tsv',sep='\t').fillna('');sites=pd.read_csv(BROAD/'SNPs.tsv',sep='\t')
  raw=pd.read_csv(BROAD/'individual_ALT_dosages.tsv.gz',sep='\t',index_col=0)
  if list(map(int,raw.columns))!=sites.pos.astype(int).tolist():raise RuntimeError('Broad Yb matrix/site order mismatch')
  source_records=[(r['sample'],set(re.split(r'[;, ]+',r.get('aliases','')))|{r['sample']}) for r in sm.to_dict('records')]
  mapping=map_sources(meta,source_records);rows=[];small=[]
  cached=pd.read_csv(YF/'individual_ALT_dosages.tsv.gz',sep='\t',index_col=0)
  old=pd.read_csv(OLD,sep='\t').sort_values('pos1');focus=old[['chrom','pos1','ref','alt']].rename(columns={'pos1':'pos'}) if 'chrom' in old else old[['pos1','ref','alt']].rename(columns={'pos1':'pos'}).assign(chrom=SPANS[locus][0])
  focus['chrom']=focus.chrom.map(canon);focus=focus[['chrom','pos','ref','alt']]
  if len(focus)!=27 or list(map(int,cached.columns))!=focus.pos.astype(int).tolist():raise RuntimeError('Frozen Yb focus identity mismatch')
  for i in range(len(meta)):
   samples=mapping[i];vals,n=union(raw.loc[samples].to_numpy(float));short,nshort=union(cached.loc[samples].to_numpy(float));rows.append(vals);small.append(short)
   conflicts.append(dict(specimen_id=meta.iloc[i].specimen_id,locus=locus,source_samples='; '.join(samples),conflicting_sites=n,focus_conflicting_sites=nshort))
  matrix=np.asarray(rows);narrow=np.asarray(small);focal=old[['left_rosina_dosage','right_rosina_dosage']].to_numpy(float).T
  flip=old.rosina_base.eq(old.ref).to_numpy();focal[:,flip]=1-focal[:,flip]
 meta['source_samples']=['; '.join(mapping[i]) for i in range(len(meta))]
 if restricted:sites,matrix=use_focus_calls(sites,matrix,focus,narrow,out)
 y=meta.phenotype.to_numpy(int);rule='carrier' if locus=='Dennis' else 'homozygote'
 geography=np.asarray(meta.geography_used_Fig3b) if locus=='Ybv' else np.array(['West']*len(meta))
 raw_scores,comparisons,orientation,_=geo.calculate(matrix,y,geography,rule)
 scores=sites.copy();scores['rule']=rule;scores['pooled_orientation']=orientation
 scores['balanced_accuracy']=comparisons['pooled'] if locus=='Dennis' or restricted else comparisons['West_only']
 scores['pooled_BA']=comparisons['pooled']
 regional=raw_scores['pooled' if locus=='Dennis' or restricted else 'West_only']
 for k in ['n_called_0','n_called_1','rate_0','rate_1']:scores[k]=regional[k]
 scores['best_within_cohort_BA']=np.maximum(regional['ALT'],regional['REF'])
 scores['status']=np.where(np.isfinite(scores.balanced_accuracy),'scored','tie_or_insufficient_callability')
 trait_scores=None;trait_summary=None
 if locus=='Dennis':
  # Reuse the reviewed Figure 3a specimen set for all three Optix traits.
  # This keeps cohort membership constant while allowing phenotype counts to vary.
  trait_frames=[];summary_rows=[];geography=np.asarray(['West']*len(meta))
  for trait,label in [('B','Band'),('D','Dennis'),('R','Ray')]:
   if not meta[trait].isin(['0','1']).all():raise RuntimeError('Optix trait unresolved in reviewed Figure 3a cohort: '+trait)
   trait_y=meta[trait].astype(int).to_numpy()
   trait_raw,trait_comparisons,trait_orientation,_=geo.calculate(matrix,trait_y,geography,'carrier')
   trait_ba=trait_comparisons['pooled']
   if trait=='D' and not np.allclose(trait_ba,scores.balanced_accuracy.to_numpy(float),equal_nan=True,atol=1e-12,rtol=0):
    raise RuntimeError('Recomputed Dennis BA differs from the frozen main calculation')
   counts={str(v):int((trait_y==v).sum()) for v in [0,1]}
   summary_rows.append(dict(trait=trait,label=label,rule='carrier',cohort='BA_Fig3a reviewed specimens',positive_n=counts['1'],negative_n=counts['0'],SNPs=len(sites),BA_ge_0_90=int(np.nansum(trait_ba>=.9))))
   trait_frames.append(pd.DataFrame(dict(chrom=sites.chrom,pos=sites.pos,ref=sites.ref,alt=sites.alt,trait=trait,label=label,rule='carrier',
    balanced_accuracy=trait_ba,pooled_orientation=trait_orientation,n_called_0=trait_raw['pooled']['n_called_0'],n_called_1=trait_raw['pooled']['n_called_1'],
    rate_0=trait_raw['pooled']['rate_0'],rate_1=trait_raw['pooled']['rate_1'],positive_n=counts['1'],negative_n=counts['0'])))
  trait_scores=pd.concat(trait_frames,ignore_index=True);trait_summary=pd.DataFrame(summary_rows)
 focus=focus.merge(scores,on=['chrom','pos','ref','alt'],how='left',validate='one_to_one',suffixes=('_original',''))
 focus['passes_updated_BA90']=focus.balanced_accuracy.ge(.9)
 if restricted:
  expected=pd.read_csv(out/'code/restricted_Ybv_27_SNP_BA_comparison.tsv',sep='\t').set_index(['chrom','pos','ref','alt'])
  values=expected.reindex(pd.MultiIndex.from_frame(focus[['chrom','pos','ref','alt']])).pooled_homozygote_BA.to_numpy(float)
  if not np.allclose(values,focus.balanced_accuracy.to_numpy(float),equal_nan=True,atol=1e-10,rtol=0):raise RuntimeError('Restricted focus BA differs from completed diagnostic')
 meta['used_in_main_BA']=True if locus=='Dennis' or restricted else meta.geography_used_Fig3b.eq('West')
 # All output matrices are independent specimens in reviewed-sheet order.
 meta.to_csv(folder/'membership.tsv',sep='\t',index=False);sites.to_csv(folder/'all_sites.tsv',sep='\t',index=False)
 pd.DataFrame(matrix,index=meta.specimen_id,columns=sites.pos).to_csv(folder/'all_ALT_dosages.tsv.gz',sep='\t')
 pd.DataFrame(narrow,index=meta.specimen_id,columns=focus.pos).to_csv(folder/'focus_ALT_dosages.tsv.gz',sep='\t')
 pd.DataFrame(focal,index=['Left (corrected)','Right (corrected)'],columns=focus.pos).to_csv(folder/'corrected_focal_ALT_dosages.tsv',sep='\t')
 scores.to_csv(folder/'all_SNP_BA.tsv.gz',sep='\t',index=False);focus.to_csv(out/'tables'/(locus+'_focus_markers_and_updated_BA.tsv'),sep='\t',index=False)
 if trait_scores is not None:
  trait_scores.to_csv(out/'tables/Optix_B_D_R_full_span_BA.tsv.gz',sep='\t',index=False,compression='gzip')
  trait_summary.to_csv(out/'tables/Optix_B_D_R_BA_cohorts.tsv',sep='\t',index=False)
 pd.DataFrame(conflicts).to_csv(out/'tables'/(locus+'_run_alias_conflicts.tsv'),sep='\t',index=False)
 selection=meta.loc[meta.used_in_main_BA];counts={str(state):int(sum(selection.phenotype.eq(state))) for state in [0,1]}
 all_individual_data(out,locus,focus)
 save(folder/'complete.json',dict(status='COMPLETE',locus=locus,individuals=len(meta),main_BA_individuals=len(selection),phenotype_counts=counts,
  SNPs=len(sites),focus_SNPs=len(focus),focus_SNPs_passing_updated_BA90=int(focus.passes_updated_BA90.sum()),
  duplicated_run_groups=sum(len(mapping[i])>1 for i in mapping),conflicting_calls=sum(r['conflicting_sites'] for r in conflicts),
  metadata_sha256=cfg['metadata_sha256'],rule=rule,cohort='lower-panel' if restricted else 'reviewed' if locus=='Dennis' else 'west',
  orientation='Best orientation within restricted cohort' if restricted else 'Final eligible pooled orientation held fixed; Ybv evaluated within West',
  calls_changed=False,new_SNP_selection=False))
 print('STATS\t'+locus+'\t'+json.dumps(counts)+'\nCOMPLETE\t'+str(folder/'complete.json'))

def all_individual_data(out,locus,focus):
 """Display available specimens, including BA-excluded references, without pooling them."""
 import numpy as np,pandas as pd,pysam
 meta=read_metadata(out/'reviewed_specimens.tsv');meta=meta.loc[~meta.globally_excluded].copy()
 path=VCF if locus=='Dennis' else YB_FULL
 engine=load(out/'code/quality_engine.py','quality') if locus=='Ybv' else None
 with pysam.VariantFile(str(path)) as vf:
  source_names=set(vf.header.samples)
  meta=meta.loc[[bool(aliases(r)&source_names) for r in meta.to_dict('records')]].reset_index(drop=True)
  mapping=map_sources(meta,[(s,{s}) for s in sorted(source_names)])
  needed=sorted({s for ss in mapping.values() for s in ss});lookup={s:i for i,s in enumerate(needed)}
  target={key(r.chrom,r.pos,r.ref,r.alt):j for j,r in enumerate(focus.itertuples())};bypos={int(r.pos):r for r in focus.itertuples()}
  vectors={j:[] for j in range(len(focus))};incompatible=np.zeros((len(needed),len(focus)),bool)
  for rec in vf:
   if rec.pos not in bypos or canon(rec.contig)!=SPANS[locus][0]:continue
   if len(rec.ref)!=1 or any(len(a)!=1 or a not in 'ACGT' for a in (rec.alts or ())):continue
   r=bypos[rec.pos];j=target[key(r.chrom,r.pos,r.ref,r.alt)]
   if rec.ref!=r.ref:raise RuntimeError('Focus reference mismatch')
   if locus=='Dennis' and tuple(rec.alts or ())!=(r.alt,):continue
   vals=[]
   for s in needed:
    c=rec.samples[s];gt=c.get('GT') or ()
    if locus=='Dennis':
     good=len(gt)==2 and all(a in (0,1) for a in gt) and c.get('DP') is not None and c['DP']>=2 and c.get('GQ') is not None and c['GQ']>=20
     value=sum(gt)/2 if good else np.nan
    else:
     dose,_=engine.quality(c,allow_pl=True)
     bases=[rec.alleles[a] for a in gt if a is not None and a<len(rec.alleles)]
     compatible=len(gt)==2 and len(bases)==2 and all(b in (r.ref,r.alt) for b in bases)
     if dose is not None and not compatible:incompatible[lookup[s],j]=True
     value=sum(b==r.alt for b in bases)/2 if dose is not None and compatible else np.nan
    vals.append(value)
   vectors[j].append(vals)
 raw=np.full((len(needed),len(focus)),np.nan)
 for j,vs in vectors.items():
  if vs:raw[:,j]=union(vs)[0]
 raw[incompatible]=np.nan
 matrix=[]
 for i in range(len(meta)):matrix.append(union(raw[[lookup[s] for s in mapping[i]]])[0])
 matrix=np.asarray(matrix);folder=out/'data'/locus
 existing=pd.read_csv(folder/'focus_ALT_dosages.tsv.gz',sep='\t',index_col=0)
 test=pd.DataFrame(matrix,index=meta.specimen_id,columns=focus.pos.astype(str))
 if not set(existing.index)<=set(test.index) or not np.allclose(existing.to_numpy(float),test.loc[existing.index].to_numpy(float),equal_nan=True,atol=1e-12,rtol=0):
  raise RuntimeError('All-individual display calls disagree with analysis focus cache: '+locus)
 selected_ids=set(pd.read_csv(folder/'membership.tsv',sep='\t').specimen_id)
 meta['included_comparison']=meta.specimen_id.isin(selected_ids)
 restricted=locus=='Ybv' and json.loads((out/'run.json').read_text()).get('Ybv_cohort')=='lower-panel'
 meta['used_in_main_BA']=meta.included_comparison & (True if locus=='Dennis' or restricted else meta.geography_used_Fig3b.eq('West'))
 meta.to_csv(folder/'all_plot_membership.tsv',sep='\t',index=False)
 test.to_csv(folder/'all_plot_focus_ALT_dosages.tsv.gz',sep='\t')

def get_data(out,locus):
 import numpy as np,pandas as pd
 folder=out/'data'/locus;done=json.loads((folder/'complete.json').read_text())
 if done['status']!='COMPLETE':raise RuntimeError('Statistics incomplete')
 meta=pd.read_csv(folder/'membership.tsv',sep='\t',dtype=str,keep_default_na=False)
 meta['phenotype']=meta.phenotype.astype(int);meta['used_in_main_BA']=meta.used_in_main_BA.map(yes)
 focus=pd.read_csv(out/'tables'/(locus+'_focus_markers_and_updated_BA.tsv'),sep='\t')
 dos=pd.read_csv(folder/'focus_ALT_dosages.tsv.gz',sep='\t',index_col=0).loc[meta.specimen_id]
 focal=pd.read_csv(folder/'corrected_focal_ALT_dosages.tsv',sep='\t',index_col=0).to_numpy(float)
 scores=pd.read_csv(folder/'all_SNP_BA.tsv.gz',sep='\t')
 if any(blocked(r) for r in meta.to_dict('records')):raise RuntimeError('Excluded taxon in calculation')
 rows=None
 if locus=='Ybv' and done.get('cohort')=='lower-panel':rows=pd.read_csv(out/'code/restricted_selected_population_rows.tsv',sep='\t',dtype=str,keep_default_na=False)
 trait_scores=pd.read_csv(out/'tables/Optix_B_D_R_full_span_BA.tsv.gz',sep='\t') if locus=='Dennis' else None
 display_meta=display_dos=None
 if locus=='Ybv':
  display_meta=pd.read_csv(folder/'all_plot_membership.tsv',sep='\t',dtype=str,keep_default_na=False)
  display_dos=pd.read_csv(folder/'all_plot_focus_ALT_dosages.tsv.gz',sep='\t',index_col=0)
 return dict(population_rows=rows,meta=meta,focus=focus,dos=dos,focal=focal,scores=scores,trait_scores=trait_scores,
  display_meta=display_meta,display_dos=display_dos,complete=done)
def pool(x,min_rate=0):
 import numpy as np
 x=np.asarray(x,float);n=np.isfinite(x).sum(axis=0)
 return np.divide(np.nansum(x,axis=0),n,out=np.full(x.shape[1],np.nan),where=(n>0)&(n>=len(x)*min_rate)),n
def short(s):return normal(s).replace('Heliconius ','H. ')
def lower_cohort(data,locus):
 # Restricted Yb metadata already contains exactly the frozen lower-panel cohort.
 m=data['meta']
 return m.loc[m.used_in_main_BA] if locus=='Dennis' else m
def main_rows(data,locus):
 import numpy as np,pandas as pd
 m=data['meta'];focus=data['focus'];dos=data['dos'].copy();focal=data['focal'].copy()
 if locus=='Dennis':
  original=np.asarray([yes(x) for x in focus.pattern_is_alt]);dos.iloc[:,:]=np.where(original,dos.to_numpy(float),1-dos.to_numpy(float));focal=np.where(original,focal,1-focal)
  allele1=np.where(original,focus.alt,focus.ref).tolist();allele0=np.where(original,focus.ref,focus.alt).tolist()
 else:allele1=focus.alt.tolist();allele0=focus.ref.tolist()
 selected=lower_cohort(data,locus);groups={};bread={};member_rows=[]
 for state in [1,0]:
  mm=selected.loc[selected.phenotype.eq(state)];ids=mm.specimen_id.tolist();values,n=pool(dos.loc[ids],.8 if locus=='Dennis' else 0)
  label=('Dennis '+('present' if state else 'absent') if locus=='Dennis' else 'Yb(v)'+('+' if state else '-'))+f' consensus (n={len(ids)})'
  bread[state]=(label,values,'consensus');groups[state]=[]
  cols=['population'] if locus=='Dennis' else ['geography_used_Fig3b','population']
  ranks=mm.groupby(cols,sort=True).size().sort_values(ascending=False,kind='stable')
  picks=list(ranks.index[:4])
  if locus=='Ybv' and data.get('population_rows') is not None:
   frozen=data['population_rows'];frozen=frozen.loc[pd.to_numeric(frozen.phenotype).eq(state)]
   picks=list(frozen[['geography','population']].drop_duplicates().itertuples(index=False,name=None))
   if len(picks)!=4:raise RuntimeError('Expected four frozen population rows per phenotype')
  if locus=='Ybv':
   for i,(geo,race) in enumerate(list(picks)):
    if not re.search(r'linaresi|linearesi',race,re.I):continue
    candidates=[k for k in ranks.index if re.search('elevatus',k[1],re.I) and k not in picks]
    if not candidates:raise RuntimeError('No same-phenotype elevatus example to restore')
    picks[i]=next((k for k in candidates if k[0]==geo),candidates[0])
  for group in picks:
   geography,population=('',group) if locus=='Dennis' else group
   take=mm.population.eq(population)
   if locus=='Ybv':take &= mm.geography_used_Fig3b.eq(geography)
   ss=mm.loc[take,'specimen_id'];v,_=pool(dos.loc[ss],.8 if locus=='Dennis' else 0)
   groups[state].append((short(population)+f' | {"D" if locus=="Dennis" else "v"}{state} | n={len(ss)}',v,'frequency'))
   member_rows.extend(dict(locus=locus,geography=geography,population=population,phenotype=state,specimen_id=s) for s in ss)
 if locus=='Ybv' and data.get('population_rows') is not None:
  ids=[r['specimen_id'] for r in member_rows]
  if len(ids)!=len(set(ids)) or set(ids)!=set(selected.specimen_id):raise RuntimeError('Displayed Ybv rows do not exactly cover BA cohort')
 rows=groups[1]+[bread[1],('Left (corrected)',focal[0],'genotype'),('Right (corrected)',focal[1],'genotype'),bread[0]]+groups[0]
 if locus=='Ybv' and data.get('population_rows') is not None:
  def ordered(state,names):
   result=[]
   for name in names:
    hits=[r for r in groups[state] if name in r[0].lower()]
    if len(hits)!=1:raise RuntimeError('Expected exactly one displayed Ybv population row: '+name)
    result.extend(hits)
   if len(result)!=len(groups[state]):raise RuntimeError('Ybv row reordering omitted a population')
   return result
  negative=ordered(0,['malleti','nov. ecu','heurippa','elevatus'])
  positive=ordered(1,['cythera','vulcanus','thelxinoe','besckei'])
  rows=negative+[bread[0],('Left (corrected)',focal[0],'genotype'),('Right (corrected)',focal[1],'genotype'),bread[1]]+positive
  rows.extend(western_cydno_pachinus_rows(data))
 return dict(sites=focus.pos.to_numpy(),a=allele1,b=allele0,rows=rows,examples=member_rows)

def western_cydno_pachinus_rows(data):
 """Western cydno/pachinus context rows; never members of the Ybv BA cohort."""
 import numpy as np,pandas as pd
 meta=data.get('display_meta');matrix=data.get('display_dos')
 if meta is None or matrix is None:return []
 country=meta.get('country',pd.Series('',index=meta.index)).str.contains(r'Colombia|Panam|Costa Rica',case=False,regex=True)
 longitude=pd.to_numeric(meta.get('longitude',pd.Series('',index=meta.index)),errors='coerce')
 west=country&(longitude<=-70)
 taxa=meta.taxon.str.contains(r'Heliconius\s+(?:cydno|pachinus)',case=False,regex=True)
 state=meta.Yb_ventral.isin(['0','1'])
 candidates=meta.loc[west&taxa&state&~meta.globally_excluded.map(yes)].copy()
 if len(candidates) and candidates.included_comparison.map(yes).any():raise RuntimeError('Western cydno/pachinus unexpectedly entered the reviewed comparison cohort')
 if candidates.empty:return []
 candidates['_v']=candidates.Yb_ventral.astype(int)
 out=[];audit=[]
 for (taxon,v),group in candidates.groupby(['taxon','_v'],sort=True):
  ids=[sid for sid in group.specimen_id if sid in matrix.index]
  if not ids:continue
  x=matrix.loc[ids].to_numpy(float);callable_fraction=np.isfinite(x).mean(axis=1)
  keep=np.asarray(callable_fraction>=.8)
  ids=[sid for sid,ok in zip(ids,keep) if ok]
  if len(ids)<3:continue
  values,n=pool(matrix.loc[ids].to_numpy(float))
  label=re.sub(r'\bssp\.\s*','',short(taxon))+f' | West v{v} | n={len(ids)}'
  out.append((label,values,'frequency'))
  cohort=dict(taxon=taxon,geography='West (coordinates: Colombia/Panama/Costa Rica, longitude <= -70)',Yb_ventral=v,
   individuals=len(ids),BA_cohort=False,sample_ids=';'.join(ids),callability_rule='individual >=80% of 27 focus SNPs')
  for pos,ncalled,value in zip(map(int,matrix.columns),n,values):
   audit.append(dict(**cohort,pos=pos,n_called=int(ncalled),ALT_frequency=float(value) if np.isfinite(value) else ''))
 data['_western_cydno_audit']=audit
 return out
def draw_main_letters(ax,data):
 import numpy as np
 from matplotlib.patches import Rectangle
 from matplotlib.colors import LinearSegmentedColormap
 import matplotlib.patheffects as pe
 cmap=LinearSegmentedColormap.from_list('main_blue_yellow',[YELLOW,'#F5F7F4',BLUE],N=257)
 rows=data['rows'];heights=np.asarray([2.7 if r[2]=='genotype' else 1 for r in rows]);centres=np.cumsum(heights)-heights/2-.5
 ax.set_xlim(-.5,len(data['sites'])-.5);ax.set_ylim(heights.sum()-.5,-.5);ax.set_xticks([])
 ax.set_yticks(centres,[r[0] for r in rows],fontsize=7);ax.tick_params(length=0,pad=7)
 for i,(label,vals,kind) in enumerate(rows):
  cy,rh=centres[i],heights[i];bold=kind=='genotype'
  for j,x in enumerate(vals):
   if not np.isfinite(x):ax.text(j,cy,'-',ha='center',va='center',color='#888888');continue
   a,b=data['a'][j],data['b'][j]
   if bold and x not in (0,.5,1):raise RuntimeError('Non-discrete focal genotype')
   if bold:
    if x==.5:
     for yy,base,col in [(cy-.43*rh,a,BLUE),(cy,b,YELLOW)]:
      ax.add_patch(Rectangle((j-.47,yy),.94,.43*rh,facecolor=col,edgecolor='none'))
      ax.text(j,yy+.215*rh,base,ha='center',va='center',color='white' if col==BLUE else '#17202A',fontsize=13,weight='bold')
    else:
     ax.add_patch(Rectangle((j-.47,cy-.43*rh),.94,.86*rh,facecolor=BLUE if x==1 else YELLOW,lw=0))
     ax.text(j,cy,a if x==1 else b,ha='center',va='center',color='white' if x==1 else '#17202A',fontsize=13,weight='bold')
    ax.add_patch(Rectangle((j-.47,cy-.43*rh),.94,.86*rh,facecolor='none',edgecolor='#555555',lw=.6))
   else:
    if kind=='frequency':ax.add_patch(Rectangle((j-.47,cy-.43),.94,.86,facecolor=cmap(x),edgecolor='white',lw=.3))
    col=BLUE if x>.5 else YELLOW if x<.5 else '#333333';base=a if x>.5 else b if x<.5 else '='
    text=ax.text(j,cy,base,ha='center',va='center',color=col if kind=='consensus' else 'white' if x>.72 else '#17202A',fontsize=10 if kind=='consensus' else 8,weight='bold')
    if kind=='consensus' and x<.5:text.set_path_effects([pe.withStroke(linewidth=.65,foreground='#8C6D00')])
  if bold:ax.get_yticklabels()[i].set_weight('bold');ax.get_yticklabels()[i].set_fontsize(9)
 for sp in ax.spines.values():sp.set_visible(False)
def annotation(ax,helper,features,native,locus,dots=None):
 import numpy as np
 from matplotlib.patches import Rectangle
 chrom,lo,hi=SPANS[locus];ax.set_xlim(lo,hi);ax.set_ylim(-.2,5.2);ax.axhline(0,c='#999999',lw=.5)
 for label,a,b in helper.INTERVALS['Dennis' if locus=='Dennis' else 'Yb(v)']:
  ax.add_patch(Rectangle((a,.1),b-a,.24,facecolor='#555555'));ax.text((a+b)/2,.5 if label!='D2' else .9,label,ha='center',fontsize=7)
 gs=[g for g in native if g['chrom']==chrom]
 if locus=='Ybv':
  c=features['cortex'];gs=[g for g in gs if not (g['start']==c['start'] and g['end']==c['end'])]
  gs+=[dict(chrom=chrom,label='cortex',**{k:c[k] for k in ['start','end','strand','exons']})]
  iv=features['ivory'];gs+=[dict(chrom=chrom,start=min(x[0] for x in iv['exons']),end=max(x[1] for x in iv['exons']),strand='+',label='ivory (continues)',exons=iv['exons'])]
 lanes=[lo-1]*3
 for g in sorted(gs,key=lambda z:z['start']):
  lane=next((i for i,end in enumerate(lanes) if g['start']>end),int(np.argmin(lanes)));lanes[lane]=g['end']+(hi-lo)*.04
  y=1.5+lane*.95;a=max(lo,g['start']);b=min(hi,g['end'])
  if b<a:continue
  col='#267B82' if g['label'].startswith('ivory') else '#444444';ax.plot([a,b],[y,y],c=col,lw=.65)
  for a1,b1 in g['exons']:
   if b1<lo or a1>hi:continue
   ax.add_patch(Rectangle((max(lo,a1),y-.11),min(hi,b1)-max(lo,a1)+1,.22,facecolor=col,lw=0))
  x1,x2=(a,b) if g['strand']=='+' else (b,a);ax.annotate('',xy=(x2,y),xytext=(x1,y),arrowprops=dict(arrowstyle='->',color=col,lw=.6))
  ax.text((a+b)/2,y+.2,g['label'],fontsize=6.5,ha='center',color=col)
 if dots is not None:
  s=dots.loc[dots.balanced_accuracy.ge(.9)];ax.scatter(s.pos,np.full(len(s),4.45),s=9,c='#D7191C',lw=0,zorder=8)
  ax.text(-.01,4.45,('D BA >=0.90' if locus=='Dennis' else 'Ybv BA >=0.90'),transform=ax.get_yaxis_transform(),ha='right',va='center',color='#D7191C',fontsize=7)
 ticks=np.linspace(lo,hi,5);ax.set_xticks(ticks,[f'{v/1e6:.3f}' for v in ticks]);ax.set_yticks([]);ax.set_xlabel(('Chromosome 18' if locus=='Dennis' else 'Chromosome 15')+' (Rosina, Mb)',fontsize=7)
 for sp in ax.spines.values():sp.set_visible(False)

def coordinate_audit(out,data,features,helper):
 import pandas as pd
 focus=data['focus'];sites=data['scores'].loc[data['scores'].balanced_accuracy.ge(.9)]
 exons=features['ivory']['exons'];first=min(exons,key=lambda e:e[0])
 intervals=helper.INTERVALS['Yb(v)'];target=next(((a,b) for label,a,b in intervals if label=='Yb(v)'),(1658079,1659072))
 keys={key(SPANS['Ybv'][0],r.pos,r.ref,r.alt) for r in focus.itertuples()}
 rows=[]
 for r in sites.itertuples():
  nearest=focus.iloc[(focus.pos-int(r.pos)).abs().argmin()]
  rows.append(dict(pos=int(r.pos),ref=r.ref,alt=r.alt,balanced_accuracy=r.balanced_accuracy,
   displayed_exact_alleles=key(SPANS['Ybv'][0],r.pos,r.ref,r.alt) in keys,
   nearest_displayed_pos=int(nearest.pos),offset_from_nearest_displayed_bp=int(r.pos)-int(nearest.pos),
   Ybv_start=target[0],Ybv_end=target[1],right_of_Ybv=int(r.pos)>target[1],
   first_ivory_exon_start=first[0],first_ivory_exon_end=first[1],right_of_first_ivory_exon=int(r.pos)>first[1]))
 columns=['pos','ref','alt','balanced_accuracy','displayed_exact_alleles','nearest_displayed_pos','offset_from_nearest_displayed_bp','Ybv_start','Ybv_end','right_of_Ybv','first_ivory_exon_start','first_ivory_exon_end','right_of_first_ivory_exon']
 pd.DataFrame(rows,columns=columns).to_csv(out/'tables/Ybv_BA90_coordinate_audit.tsv',sep='\t',index=False)
 print('YBV_BA90_COORDINATES\t'+json.dumps(rows,allow_nan=False))

def main_figures(out,data,helper,features,native,windows):
 import numpy as np,matplotlib.pyplot as plt
 from matplotlib.patches import ConnectionPatch
 top=max(float(w[['left_Ho','right_Ho']].max().max()) for w in windows.values());ylim=min(1.,max(.3,float(np.ceil(top*1.1*10)/10)))
 all_examples=[]
 for dots,filename in [(True,'01_Figure3_main_BA_dots.pdf'),(False,'02_Supplement_full_BA.pdf')]:
  fig=plt.figure(figsize=(14,17.5));audit=[]
  for i,locus in enumerate(['Dennis','Ybv']):
   d=data[locus];letters=main_rows(d,locus);base=.515 if i==0 else .035;left=.35;width=.62
   fig.text(.025,base+.46,'A  Optix / Dennis' if i==0 else 'B  Ivory / Yb(v)',fontsize=13,weight='bold')
   cohort=d['complete']['phenotype_counts'];description=('Carrier BA' if locus=='Dennis' else ('Restricted-population homozygote BA' if d['complete'].get('cohort')=='lower-panel' else 'Western homozygote BA'))+f' | + n={cohort["1"]}, - n={cohort["0"]}'
   fig.text(left,base+.45,description,fontsize=7)
   if not dots:
    ba=fig.add_axes([left,base+.37,width,.075]);s=d['scores'];ba.scatter(s.pos,s.balanced_accuracy,s=6,c='#333333',alpha=.65,lw=0)
    ba.axhline(.5,c='#aaaaaa',ls=':',lw=.5);ba.axhline(.9,c='#D7191C',ls='--',lw=.6);ba.set_xlim(*SPANS[locus][1:]);ba.set_ylim(0,1.02);ba.set_ylabel('Phenotype BA');ba.tick_params(axis='x',labelbottom=False);ba.spines[['top','right']].set_visible(False)
   ho=fig.add_axes([left,base+(.325 if dots else .265),width,.12 if dots else .09])
   for q,ls in [(.95,'-'),(.99,'--')]:
    w=windows[(locus,q)]
    for side,col in [('left','#1764AB'),('right','#D04B35')]:ho.plot(w.mid,w[side+'_Ho'],c=col,ls=ls,lw=1.25 if q==.95 else .8,label=side.title()+f' | P >= {q}')
   ho.set_xlim(*SPANS[locus][1:]);ho.set_xticks(np.linspace(*SPANS[locus][1:],5));ho.set_ylim(0,ylim);ho.set_yticks(np.round(np.arange(0,ylim+1e-8,.05),2));ho.set_ylabel('Heterozygous fraction\namong callable SNPs');ho.tick_params(axis='x',labelbottom=False);ho.legend(fontsize=7,ncol=2,loc='upper right');ho.spines[['top','right']].set_visible(False)
   ann=fig.add_axes([left,base+(.225 if dots else .185),width,.085 if dots else .065]);annotation(ann,helper,features,native,locus,d['scores'] if dots else None)
   heat=fig.add_axes([left,base+.010,width,.19 if dots else .16]);draw_main_letters(heat,letters)
   for j,pos in enumerate(letters['sites']):fig.add_artist(ConnectionPatch(xyA=(int(pos),0),coordsA=ann.transData,xyB=(j,-.5),coordsB=heat.transData,color='#BFC3C7',lw=.35,alpha=.5,clip_on=False))
   fig.text(left,base-.004,'Blue: original Dennis-associated allele; yellow: other allele' if locus=='Dennis' else 'Blue: ALT; yellow: REF. Consensus and BA use the same restricted population rows.',fontsize=7)
   audit.append(dict(locus=locus,cohort=cohort,rows=[r[0] for r in letters['rows']],genes=[g['label'] for g in native if g['chrom']==SPANS[locus][0]],focus_SNPs=len(letters['sites']),xlim_bp=list(SPANS[locus][1:]),x_ticks_bp=[int(x) for x in np.linspace(*SPANS[locus][1:],5)],Ho_y_tick_step=.05,first_window_midpoint_bp=float(windows[(locus,.95)].mid.iloc[0]),lower_panel_counts={str(v):int(sum(lower_cohort(d,locus).phenotype.eq(v))) for v in [0,1]}))
   if dots:all_examples.extend(letters['examples'])
  fig.text(.025,.004,'Reviewed specimen metadata; one contribution per specimen. Red dots: updated Dennis/Ybv BA >=0.90. Western cydno/pachinus rows are display-only, outside BA and consensus. Existing 24 Dennis / 27 Yb(v) focus SNP identities retained.\nHo: unchanged corrected calls, 10-kb windows / 2.5-kb slide, common callable denominator. No new calling or mixture fit.\nTwo-letter cells are unphased heterozygotes; = means a pooled 50:50 tie. Ivory continues beyond the fixed genomic view.',fontsize=6.5)
  fig.savefig(out/filename);plt.close(fig);save(out/'tables'/filename.replace('.pdf','.json'),dict(panels=audit,Ho_ylim=[0,ylim]))
 import pandas as pd
 pd.DataFrame(all_examples).to_csv(out/'tables/main_example_membership.tsv',sep='\t',index=False)

def specimen_label(r):
 codes=' '.join(k+str(r.get(k,'?')) for k in ['B','D','R'])+' d'+str(r.get('Yb_dorsal','?'))+' v'+str(r.get('Yb_ventral','?'))
 return normal(r['specimen_name'])+' | '+short(r['population'])+' | '+normal(r['run_accessions'])+' | '+codes+' | BA='+('yes' if yes(r.get('used_in_main_BA')) else 'no')

def heat_page(pdf,rows,positions,title,locus,frequency=False,page=None):
 import numpy as np,matplotlib.pyplot as plt
 from matplotlib.colors import LinearSegmentedColormap,Normalize
 from matplotlib.patches import ConnectionPatch,Rectangle
 labels=[r[0] for r in rows];matrix=np.asarray([r[1] for r in rows],float);nr=len(rows)
 label_width=max(4.8,max(map(len,labels))*.055+0.3);width=max(14,label_width+7.1);height=max(5.2,3+nr*.24);left=label_width/width
 fig=plt.figure(figsize=(width,height));heat=fig.add_axes([left,.15,.98-left,.67]);loc=fig.add_axes([left,.88,.98-left,.055])
 cmap=LinearSegmentedColormap.from_list('previous_genotype_palette',SUPP,N=257);cmap.set_bad('#C8C8C8')
 heat.pcolormesh(np.arange(len(positions)+1)-.5,np.arange(nr+1)-.5,np.ma.masked_invalid(matrix),vmin=0,vmax=1,cmap=cmap,shading='flat',rasterized=False)
 heat.set_ylim(nr-.5,-.5);heat.set_xlim(-.5,len(positions)-.5);heat.set_yticks(np.arange(nr),labels,fontsize=7);heat.set_xticks([]);heat.tick_params(length=0,pad=8)
 for i,row in enumerate(rows):
  if row[2]=='genotype':
   heat.axhline(i-.5,c='#666666',lw=.5);heat.axhline(i+.5,c='#666666',lw=.5);heat.get_yticklabels()[i].set_weight('bold')
 lo,hi=int(min(positions))-180,int(max(positions))+180;loc.set_xlim(lo,hi);loc.set_ylim(0,1.3);loc.axhline(.1,c='#999999',lw=.5)
 intervals=[('D1',1215210,1216124),('D2',1218557,1221494)] if locus=='Dennis' else [('Yb(v)',1658079,1659072)]
 for label,a,b in intervals:loc.add_patch(Rectangle((a,.2),b-a,.2,color='#666666'));loc.text((a+b)/2,.65,label,ha='center',fontsize=8)
 loc.set_xticks(np.linspace(lo,hi,4),[f'{x/1e6:.4f}' for x in np.linspace(lo,hi,4)],fontsize=7);loc.set_yticks([])
 for j,pos in enumerate(positions):fig.add_artist(ConnectionPatch(xyA=(int(pos),.1),coordsA=loc.transData,xyB=(j,-.5),coordsB=heat.transData,color='#BFC3C7',lw=.3,alpha=.35,clip_on=False))
 for ax in [heat,loc]:
  for sp in ax.spines.values():sp.set_visible(False)
 fig.text(.025,.975,title+(f' | page {page}' if page else ''),fontsize=12,weight='bold')
 bar=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=cmap),cax=fig.add_axes([max(left,.65),.085,min(.27,.98-max(left,.65)),.014]),orientation='horizontal',ticks=[0,.5,1])
 bar.ax.tick_params(labelsize=7);bar.set_label('ALT frequency' if frequency else 'ALT dosage: REF/REF = 0; REF/ALT = 0.5; ALT/ALT = 1',fontsize=7)
 fig.text(.025,.025,'Previous palette: pale blue / gold / deep purple. Grey = missing. All specimen labels are a single line.\n'+('Each frequency pools only reviewed comparison members; n_called is exported for every marker. A pooled frequency of 0.5 need not represent heterozygous individuals.' if frequency else 'BA=no references are display-only and contribute to no phenotype statistics. Corrected Left/Right shown separately; tissues are not independent specimens.'),fontsize=7)
 pdf.savefig(fig);plt.close(fig)

def supplements(out,data):
 import numpy as np,pandas as pd
 from matplotlib.backends.backend_pdf import PdfPages
 freq_exports=[];all_counts={}
 for locus,prefix in [('Dennis','03_Dennis_all_individuals'),('Ybv','05_Ybv_all_individuals')]:
  d=data[locus];folder=out/'data'/locus;meta=pd.read_csv(folder/'all_plot_membership.tsv',sep='\t',dtype=str,keep_default_na=False)
  matrix=pd.read_csv(folder/'all_plot_focus_ALT_dosages.tsv.gz',sep='\t',index_col=0);meta=meta.sort_values(['taxon','population','specimen_name','specimen_id'])
  focal=[('Left (corrected)',d['focal'][0],'genotype'),('Right (corrected)',d['focal'][1],'genotype')];individuals=[(specimen_label(r),matrix.loc[r['specimen_id']].to_numpy(float),'individual') for r in meta.to_dict('records')]
  with PdfPages(out/(prefix+'.pdf')) as pdf:
   for i in range(0,len(individuals),55):heat_page(pdf,individuals[i:i+55]+focal,d['focus'].pos,locus+' | all available reviewed individuals',locus,page=i//55+1)
  all_counts[locus]=len(meta)
  rows=[];m=d['meta'];cols=['population','phenotype'] if locus=='Dennis' else ['geography_used_Fig3b','population','phenotype']
  for group,g in m.groupby(cols,sort=True):
   values,n=pool(d['dos'].loc[g.specimen_id],.8 if locus=='Dennis' else 0);geo='' if locus=='Dennis' else str(group[0])+' | '
   pop=group[0] if locus=='Dennis' else group[1];state=group[-1]
   rows.append((geo+short(pop)+f' | {"D" if locus=="Dennis" else "v"}{state} | n={len(g)}',values,'frequency'))
   for j,r in enumerate(d['focus'].itertuples()):freq_exports.append(dict(locus=locus,population=pop,geography='' if locus=='Dennis' else group[0],phenotype=state,pos=r.pos,ref=r.ref,alt=r.alt,n_total=len(g),n_called=int(n[j]),ALT_frequency=values[j]))
  filename='04_Dennis_race_frequencies.pdf' if locus=='Dennis' else '06_Ybv_race_frequencies.pdf'
  with PdfPages(out/filename) as pdf:
   for i in range(0,len(rows),50):heat_page(pdf,rows[i:i+50]+focal,d['focus'].pos,locus+' | reviewed comparison frequencies by population',locus,frequency=True,page=i//50+1)
 with PdfPages(out/'07_BA_selection_individuals.pdf') as pdf:
  for locus,d in data.items():
   m=d['meta'].loc[d['meta'].used_in_main_BA].sort_values(['population','phenotype','specimen_name']);rows=[]
   for r in m.to_dict('records'):rows.append((specimen_label(r),d['dos'].loc[r['specimen_id']].to_numpy(float),'individual'))
   focal=[('Left (corrected)',d['focal'][0],'genotype'),('Right (corrected)',d['focal'][1],'genotype')]
   for i in range(0,len(rows),55):heat_page(pdf,rows[i:i+55]+focal,d['focus'].pos,locus+' | individuals used in main BA',locus,page=i//55+1)
 with PdfPages(out/'08_Phenotype_frequencies.pdf') as pdf:
  for locus,d in data.items():
   m=d['meta'].loc[d['meta'].used_in_main_BA];rows=[]
   for state in [1,0]:
    ids=m.loc[m.phenotype.eq(state),'specimen_id'];values,n=pool(d['dos'].loc[ids],.8 if locus=='Dennis' else 0)
    rows.append((('Dennis' if locus=='Dennis' else ('Restricted Yb(v)' if d['complete'].get('cohort')=='lower-panel' else 'Western Yb(v)'))+f' | phenotype {state} | n={len(ids)}',values,'frequency'))
   rows += [('Left (corrected)',d['focal'][0],'genotype'),('Right (corrected)',d['focal'][1],'genotype')]
   heat_page(pdf,rows,d['focus'].pos,locus+' | frequencies pooled by phenotype',locus,frequency=True)
 pd.DataFrame(freq_exports).to_csv(out/'tables/race_frequencies_and_denominators.tsv.gz',sep='\t',index=False)
 return all_counts

def optix_trait_BA_figure(out,data):
 import numpy as np,matplotlib.pyplot as plt,pandas as pd
 from matplotlib.backends.backend_pdf import PdfPages
 scores=data['Dennis']['trait_scores']
 if scores is None:raise RuntimeError('Optix Band/Dennis/Ray BA data missing')
 palette={'B':'#1F77B4','D':'#4D4D4D','R':'#009E73'};labels={'B':'Band','D':'Dennis','R':'Ray'}
 summary=pd.read_csv(out/'tables/Optix_B_D_R_BA_cohorts.tsv',sep='\t')
 path=out/'09_Optix_B_D_R_full_span_BA.pdf'
 fig=plt.figure(figsize=(12,7.2))
 for i,trait in enumerate(['B','D','R']):
  ax=fig.add_axes([.11,.70-i*.215,.84,.16]);s=scores.loc[scores.trait.eq(trait)]
  ax.scatter(s.pos,s.balanced_accuracy,s=2.2,c=palette[trait],alpha=.48,lw=0,rasterized=True)
  ax.axhline(.5,c='#999999',ls=':',lw=.6);ax.axhline(.9,c='#D7191C',ls='--',lw=.65)
  ax.set_xlim(*SPANS['Dennis'][1:]);ax.set_ylim(0,1.02);ax.set_ylabel('BA',fontsize=8)
  row=summary.loc[summary.trait.eq(trait)].iloc[0]
  ax.set_title(f'{labels[trait]} | carrier rule | + n={int(row.positive_n)}, − n={int(row.negative_n)}',loc='left',fontsize=9,pad=3)
  ax.tick_params(axis='x',labelsize=7);ax.tick_params(axis='y',labelsize=7);ax.spines[['top','right']].set_visible(False)
  if trait!='R':ax.tick_params(axis='x',labelbottom=False)
  else:
   ticks=np.linspace(*SPANS['Dennis'][1:],5);ax.set_xticks(ticks,[f'{x/1e6:.3f}' for x in ticks]);ax.set_xlabel('Rosina chromosome 18 (Mb)',fontsize=8)
 fig.text(.04,.965,'Optix / Dennis region | Band, Dennis, and Ray phenotype BA',fontsize=13,weight='bold')
 fig.text(.11,.925,'Same reviewed Figure 3a individuals for all traits; one pooled carrier orientation per SNP; DP ≥2 and GQ ≥20; ≥80% callability per phenotype class.',fontsize=8)
 fig.text(.11,.035,'Red dashed line: BA=0.90. Existing marker selection and the Dennis main-figure BA dots are unchanged. Full-span calculations use existing calls only.',fontsize=7)
 fig.savefig(path);plt.close(fig)
 return path

def render(out):
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams.update({'pdf.fonttype':42,'font.size':8,'font.family':'DejaVu Sans'})
 cfg=json.loads((out/'run.json').read_text())
 if sha(out/'reviewed_specimens.tsv')!=cfg['metadata_sha256']:raise RuntimeError('Metadata snapshot changed')
 data={locus:get_data(out,locus) for locus in ['Dennis','Ybv']}
 helper=load(out/'code/panel_helper.py','panel');old=load(out/'code/Ho_helper.py','ho');features=json.loads((out/'code/Ybv_gene_features.json').read_text());native=helper.genes();windows={}
 for locus in data:
  path=RUN/'optix_consistent_v2/reciprocal_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz' if locus=='Dennis' else HO/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz'
  table=old.normalize(pd.read_csv(path,sep='\t'),'Optix' if locus=='Dennis' else 'Ivory_Yb')
  for q in [.95,.99]:
   w=old.windows(table,*SPANS[locus][1:],'bp',10000,q);windows[(locus,q)]=w
   if not (w.end-w.start+1).eq(10000).all() or not w.start.diff().dropna().eq(2500).all():raise RuntimeError('Ho geometry changed')
   w.to_csv(out/'tables'/f'{locus}_Ho_10kb_step2500_P{q}.tsv.gz',sep='\t',index=False)
 coordinate_audit(out,data['Ybv'],features,helper)
 main_figures(out,data,helper,features,native,windows);all_counts=supplements(out,data)
 optix_trait_BA_figure(out,data)
 cydno_audit=data['Ybv'].get('_western_cydno_audit',[])
 if cydno_audit:pd.DataFrame(cydno_audit).to_csv(out/'tables/Ybv_western_cydno_pachinus_display_only.tsv',sep='\t',index=False)
 audit=[]
 reviewed=read_metadata(out/'reviewed_specimens.tsv')
 raw=pd.read_csv(out/'reviewed_specimens.tsv',sep='\t',dtype=str,keep_default_na=False)
 changed=raw.population.ne(reviewed.population)
 pd.DataFrame(dict(specimen_id=raw.loc[changed,'specimen_id'],source_population=raw.loc[changed,'population'],plot_population=reviewed.loc[changed,'population'])).to_csv(out/'tables/population_label_normalizations.tsv',sep='\t',index=False)
 for r in reviewed.to_dict('records'):
  audit.append(dict(specimen_id=r['specimen_id'],specimen_name=r['specimen_name'],taxon=r['taxon'],population=r['population'],
   globally_excluded=r['globally_excluded'],Dennis_used=r['specimen_id'] in set(data['Dennis']['meta'].specimen_id),
   Ybv_used_in_regional_context=r['specimen_id'] in set(data['Ybv']['meta'].specimen_id),
   Ybv_used_in_main_BA=r['specimen_id'] in set(data['Ybv']['meta'].loc[data['Ybv']['meta'].used_in_main_BA,'specimen_id'])))
 pd.DataFrame(audit).to_csv(out/'tables/final_specimen_membership.tsv',sep='\t',index=False)
 summary={locus:d['complete'] for locus,d in data.items()}
 consensus=[]
 for locus,d in data.items():
  shutil.copy2(out/'data'/locus/'all_SNP_BA.tsv.gz',out/'tables'/(locus+'_all_SNP_BA.tsv.gz'))
  mm=lower_cohort(d,locus)
  for state in [1,0]:
   ss=mm.loc[mm.phenotype.eq(state),'specimen_id'];values,n=pool(d['dos'].loc[ss],.8 if locus=='Dennis' else 0)
   for j,r in enumerate(d['focus'].itertuples()):consensus.append(dict(locus=locus,phenotype=state,pos=r.pos,ref=r.ref,alt=r.alt,n_total=len(ss),n_called=int(n[j]),ALT_frequency=values[j]))
 pd.DataFrame(consensus).to_csv(out/'tables/main_consensus_frequencies.tsv',sep='\t',index=False)
 pd.DataFrame([dict(locus=k,negative_n=v['phenotype_counts']['0'],positive_n=v['phenotype_counts']['1'],all_genoplot_individuals=all_counts[k],focus_SNPs=v['focus_SNPs'],focus_SNPs_passing_current_BA90=v['focus_SNPs_passing_updated_BA90']) for k,v in summary.items()]).to_csv(out/'tables/cohort_and_marker_summary.tsv',sep='\t',index=False)
 readme='Figure 3 review package for Arnaud\n\n01 Main figure: red Dennis/Ybv BA>=0.90 dots, unchanged Ho, full gene annotations and nucleotide focus panels. Bright blue/yellow.\n02 Supplement: identical layout with full BA scatter plots for the main Dennis and Ybv comparisons.\n03/05 All available reviewed individual genotypes, paginated, single-line names; corrected gynandromorph Left/Right states are at the bottom of each page. BA=no references are display-only.\n04/06 Population frequencies from the reviewed comparison cohorts; corrected gynandromorph Left/Right states are at the bottom of each page. Yb(v) retains East/West labels.\n07 Individuals actually contributing to main BA, with Left/Right at the bottom of each page. 08 Frequencies pooled by phenotype with Left/Right at the bottom.\n09 Full-span Optix BA for Band, Dennis, and Ray. Traits use the same reviewed Figure 3a specimen membership and carrier rule; trait-specific positive/negative counts are shown.\n\nThe Ivory/Yb(v) main lower genotype panel adds western cydno/pachinus population frequencies as display-only context, outside BA and consensus cohorts. Species, ventral phenotype, per-individual callability, and per-SNP denominators are recorded in tables/Ybv_western_cydno_pachinus_display_only.tsv.\nSupplement genotype palette: pale blue REF/REF, gold heterozygote, deep purple ALT/ALT; grey missing. Frequency midpoint need not mean every specimen is heterozygous.\nMain colours: Dennis blue=original associated allele/yellow=other; Yb(v) blue=ALT/yellow=REF. Two letters within a box are unphased heterozygotes. Consensus is a majority base per site, not a phased haplotype.\n\nBA uses >=80% callability in each phenotype class. Dennis main: carrier rule, reviewed cohort. Yb(v): homozygote rule in the explicitly frozen restricted lower-panel cohort; allele direction is chosen within that cohort. Allele orientation and per-class denominators are recorded for every score.\nOne biological specimen contributes once. Qualified calls across multiple accessions are united when concordant; conflicts are missing and exported. Genotypes are read from completed call sets; no calling or mixture fitting. Legacy Optix and newer call sets retain their documented calling histories. Optix quality: DP>=2, reported GQ>=20. Yb: existing GQ/PL quality engine, unchanged.\nExisting 24 Dennis/27 Yb focus marker identities are retained without reselection. Tables report updated BA and whether each still exceeds 0.90; red dots show all currently passing SNPs across the full region.\nMain Yb lower consensus, example populations and main BA use the exact same frozen restricted membership. Example selection restores the prior geography/race grouping with elevatus replacing linaresi. Cythera ssp/spp spelling variants are merged. Ho y-tick spacing is explicitly 0.05. Exact red-dot positions and first ivory exon relationships are exported in tables/Ybv_BA90_coordinate_audit.tsv.\nFull spans: Optix 1,060,000-1,240,000; Yb 1,630,000-1,692,000. Ho: 10kb windows/2.5kb slide; same corrected focal calls and common denominator. Ivory annotation continues beyond the view.\n\nReviewed metadata snapshot and all membership/denominator tables are included. Numata and pardalinus are excluded from all plots and calculations (retained only in metadata/provenance). Existing references excluded from BA may appear only in all-individual displays, labelled BA=no. No metadata-only specimen is fabricated as a sequenced individual.\n'
 (out/'README.txt').write_text(readme.replace('excluded from every output and calculation','excluded from all plots and calculations (retained only in metadata/provenance)'))
 pdfs=sorted(out.glob('*.pdf'))
 if len(pdfs)!=9 or any(p.stat().st_size==0 for p in pdfs):raise RuntimeError('Expected nine non-empty PDFs')
 pdf_pages={}
 try:
  from pypdf import PdfReader
  pdf_pages={p.name:len(PdfReader(p).pages) for p in pdfs}
 except ImportError:pass
 completed=dict(status='COMPLETE',PDFs=[str(p) for p in pdfs],pdf_pages=pdf_pages,cohorts=summary,
  all_individual_counts=all_counts,metadata_sha256=cfg['metadata_sha256'],new_calls=False,focal_calls_changed=False,
  SNPs_reselected=False,main_palette=[YELLOW,BLUE],supplement_palette=SUPP,spans=SPANS)
 archive=out/'Figure3_for_Arnaud.zip.partial'
 with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
  for p in [*pdfs,out/'README.txt',out/'reviewed_specimens.tsv',*sorted((out/'tables').glob('*'))]:
   if p.is_file():z.write(p,str(p.relative_to(out)))
  z.writestr('complete.json',json.dumps(completed,indent=2)+'\n')
 archive.replace(out/'Figure3_for_Arnaud.zip');save(out/'complete.json',completed)
 print('STATUS\tCOMPLETE\nMAIN_PDF\t'+str(out/'01_Figure3_main_BA_dots.pdf')+'\nSHARE_ZIP\t'+str(out/'Figure3_for_Arnaud.zip'))
 print('COHORTS\t'+json.dumps({k:v['phenotype_counts'] for k,v in summary.items()}))
def main():
 import os
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--metadata',type=Path,default=Path(__file__).resolve().parents[3]/'data/figure3_review/final_specimens.tsv')
 p.add_argument('--render-only',action='store_true');p.add_argument('--yb-cohort',choices=['lower-panel','west'],default='lower-panel');p.add_argument('--output',type=Path,default=ROOT/'analysis/Figure3_for_Arnaud_final_20261027');p.add_argument('--submit',action='store_true');p.add_argument('--overwrite',action='store_true');p.add_argument('--stage',choices=['stats','render']);a=p.parse_args()
 if a.submit and a.render_only:prepare_render_only(a)
 elif a.submit:prepare(a)
 elif a.stage=='stats':statistics(a.output,'Dennis' if int(os.environ['SLURM_ARRAY_TASK_ID'])==0 else 'Ybv')
 elif a.stage=='render':render(a.output)
 else:p.error('Use --submit')
if __name__=='__main__':main()
