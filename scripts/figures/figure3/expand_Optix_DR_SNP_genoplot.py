#!/usr/bin/env python3
"""More inclusive Optix columns with the historical D/R BA display preserved.

Selection = previous 176 exact alleles union archived D or R BA >= 0.93.
Archived highlights remain BA > 0.99 AND historical focal-compatible pattern.
No BA recalculation, allele flipping, new calling, or genotype-based exclusions.
"""
import argparse,hashlib,importlib.util,json,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
SOURCE=ROOT/'analysis/supp_genoplots_LR_cydno_historical_BDR_20261027'
OUTPUT=ROOT/'analysis/Optix_DR_inclusive_SNPs_20261102'
LEGACY=ROOT/'external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf'
KEYS=['chrom','pos','ref','alt']
def load(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def keys(t):
 return [(str(c) if str(c) not in ['Hmr18','18'] else 'NC_140406.1',int(p),str(r),str(a)) for c,p,r,a in t[KEYS].itertuples(index=False,name=None)]
def flags(s):
 if s.dtype==bool:return s
 t=s.astype(str).str.lower().map({'true':True,'false':False})
 if t.isna().any():raise RuntimeError('Invalid archived highlight flag')
 return t

def select_sites(previous,scores,cutoff):
 import pandas as pd,numpy as np
 old=previous.copy();old['chrom']=[k[0] for k in keys(old)]
 if old.duplicated(KEYS).any():raise RuntimeError('Duplicate old genoplot alleles')
 dr=scores.loc[scores.trait.isin(['D','R'])].copy();dr['chrom']=[k[0] for k in keys(dr)]
 if dr.duplicated(['trait']+KEYS).any():raise RuntimeError('Duplicate archived trait/allele')
 selected=set(keys(old))|set(keys(dr.loc[dr.balanced_accuracy.ge(cutoff)]))
 sites=pd.DataFrame(sorted(selected),columns=KEYS);prior=dict(zip(keys(old),old.balanced_accuracy))
 sites['previously_displayed']=[k in prior for k in keys(sites)]
 sites['updated_Dennis_selection_BA']=[prior.get(k,np.nan) for k in keys(sites)]
 for trait in ['D','R']:
  g=dr.loc[dr.trait.eq(trait)];lookup=dict(zip(keys(g),g.balanced_accuracy));hi=dict(zip(keys(g),flags(g.historical_highlight)))
  sites['archived_'+trait+'_BA']=[lookup.get(k,np.nan) for k in keys(sites)]
  sites['archived_'+trait+'_highlight']=[bool(hi.get(k,False)) for k in keys(sites)]
 # Compatibility only: rendering uses an explicit union label, not this field.
 sites['balanced_accuracy']=sites[['updated_Dennis_selection_BA','archived_D_BA','archived_R_BA']].max(axis=1)
 if not set(keys(dr.loc[flags(dr.historical_highlight)]))<=selected:raise RuntimeError('Historical highlight outside declared selection')
 return sites

def prepare(a):
 repo=Path(__file__).resolve().parents[3]
 library=repo/'scripts/figures/figure3/remake_Dennis_DR_supplement.py'
 old=load(a.source/'worker.py','old_supp');helper=load(old.MORPH/'figure3_helper.py','input_helper');base=load(a.source/'append_helper.py','input_focal')
 paths=[a.source/n for n in ['complete.json','worker.py','append_helper.py','run.json',
  'tables/Dennis_BA093_selected_SNPs.tsv','tables/Dennis_BA093_frequencies.tsv.gz',
  'tables/corrected_LR_calls.tsv','tables/Optix_September_B_D_R_BA.tsv.gz']]
 paths += [old.MORPH/'worker.py',old.MORPH/'figure3_helper.py',base.SOURCE/'reviewed_specimens.tsv',helper.VCF,LEGACY,library,
  base.OPTIX,base.HO/'worker.py',base.SOURCE/'tables/Dennis_focus_markers_and_updated_BA.tsv',base.SOURCE/'data/Dennis/corrected_focal_ALT_dosages.tsv']
 for p in paths:
  if not p.is_file():raise FileNotFoundError('Missing documented input: '+str(p))
 if json.loads((a.source/'complete.json').read_text()).get('status')!='COMPLETE':raise RuntimeError('Source incomplete')
 if a.output.exists():raise RuntimeError('Preserving previous output: '+str(a.output))
 a.output.mkdir(parents=True);(a.output/'logs').mkdir();(a.output/'tables').mkdir()
 shutil.copy2(__file__,a.output/'worker.py');shutil.copy2(library,a.output/'track_helper.py')
 # Large VCFs: record stat identity, avoid hashing gigabytes on the login node.
 large={str(p):dict(size=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in [helper.VCF,LEGACY]}
 save(a.output/'run.json',dict(source=str(a.source),cutoff=.93,hashes={str(p):sha(p) for p in paths if p not in [helper.VCF,LEGACY]},
  VCFs=large,selection='previous 176 union every archived D or R BA >=0.93; no focal gate; no cap',BA_recalculated=False))
 cmd=[sys.executable,str(a.output/'worker.py'),'--run','--output',str(a.output)]
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=12G','--cpus-per-task=1','--job-name=Optix_DRwide',
  '--output='+str(a.output/'logs/plot.%j.out'),'--error='+str(a.output/'logs/plot.%j.err'),'--wrap='+shlex.join(cmd)],text=True).strip().split(';')[0]
 save(a.output/'jobs.json',dict(plot=job));print('JOB\t'+job+'\nOUTPUT\t'+str(a.output))

def read_calls(sites,meta,helper,out):
 import numpy as np,pandas as pd,pysam
 target={k:j for j,k in enumerate(keys(sites))};available=set();found={};samples={};vectors={}
 # The merged calls have priority. Legacy records are used only for exact alleles
 # absent from the merged VCF, never to patch a low-quality merged genotype.
 for priority,path in enumerate([helper.VCF,LEGACY]):
  with pysam.VariantFile(str(path)) as vf:
   names=list(vf.header.samples);available.update(names);samples[priority]=names
   for rec in vf:
    if len(rec.alts or ())!=1:continue
    k=(helper.canon(rec.contig),rec.pos,rec.ref,rec.alts[0])
    if k not in target:continue
    if k in found:
     if found[k]==priority:raise RuntimeError('Duplicate exact allele in VCF: '+str(k))
     continue
    found[k]=priority;vals=[]
    for sample in names:
     c=rec.samples[sample];gt=c.get('GT') or ();dp=c.get('DP');gq=c.get('GQ')
     good=len(gt)==2 and all(a in (0,1) for a in gt) and dp is not None and dp>=2 and gq is not None and gq>=20
     vals.append(sum(gt)/2 if good else np.nan)
    vectors[k]=dict(zip(names,vals))
 missing=[k for k in target if k not in found]
 pd.DataFrame([dict(zip(KEYS,k)) for k in missing],columns=KEYS).to_csv(out/'tables/absent_exact_alleles.tsv',sep='\t',index=False)
 if missing:raise RuntimeError(f'{len(missing)} selected exact alleles absent from both documented VCFs; see tables/absent_exact_alleles.tsv')
 mapping={i:sorted(helper.aliases(r)&available) for i,r in enumerate(meta.to_dict('records'))}
 matrix=np.full((len(meta),len(sites)),np.nan);audit=[]
 for k,j in target.items():
  for i in range(len(meta)):
   vals=[vectors[k].get(s,np.nan) for s in mapping[i]]
   if vals:matrix[i,j]=helper.union(np.asarray(vals).reshape(-1,1))[0][0]
  audit.append(dict(zip(KEYS,k),source=str([helper.VCF,LEGACY][found[k]]),called_specimens=int(np.isfinite(matrix[:,j]).sum())))
 pd.DataFrame(audit).to_csv(out/'tables/exact_allele_call_sources.tsv',sep='\t',index=False)
 return matrix,mapping

def run(out):
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 from matplotlib.backends.backend_pdf import PdfPages
 cfg=json.loads((out/'run.json').read_text());src=Path(cfg['source'])
 for p,h in cfg['hashes'].items():
  if sha(Path(p))!=h:raise RuntimeError('Source changed: '+p)
 for p,identity in cfg['VCFs'].items():
  st=Path(p).stat()
  if dict(size=st.st_size,mtime_ns=st.st_mtime_ns)!=identity:raise RuntimeError('VCF changed: '+p)
 old=load(src/'worker.py','source_style');base=load(src/'append_helper.py','focal_source');lib=load(out/'track_helper.py','track_style')
 morph=load(old.MORPH/'worker.py','morph_renderer');helper=load(old.MORPH/'figure3_helper.py','figure_helper')
 previous=pd.read_csv(src/'tables/Dennis_BA093_selected_SNPs.tsv',sep='\t')
 scores=pd.read_csv(src/'tables/Optix_September_B_D_R_BA.tsv.gz',sep='\t');scores['historical_highlight']=flags(scores.historical_highlight)
 sites=select_sites(previous,scores,cfg['cutoff']);sites.to_csv(out/'tables/selected_SNPs.tsv',sep='\t',index=False)
 exp=pd.read_csv(src/'tables/Dennis_BA093_frequencies.tsv.gz',sep='\t',keep_default_na=False,na_values=[''])
 reviewed=helper.read_metadata(base.SOURCE/'reviewed_specimens.tsv');meta=helper.selected(reviewed,'Dennis').reset_index(drop=True)
 meta['normal_population']=meta.population.map(morph.population_name)
 meta['morph']=[' '.join(t+str(r[t]) for t in ['B','D','R']) for r in meta.to_dict('records')]
 groups=[];membership=[]
 for _,g in exp.groupby(['population','morph','geography','phenotype'],sort=False,dropna=False):
  r=g.iloc[0];members=meta.loc[meta.normal_population.eq(morph.population_name(r.population))&meta.morph.eq(r.morph)&meta.D.eq(str(int(float(r.phenotype))))]
  if len(members)!=int(r.n_total):raise RuntimeError('Reviewed specimen count differs from saved population: '+str(r.population))
  groups.append((r,g,members.index.tolist()))
  for i in members.index:membership.append(dict(row=len(groups)-1,specimen_id=meta.loc[i,'specimen_id'],population=r.population,morph=r.morph))
 ids=sorted({i for _,_,group in groups for i in group});meta=meta.loc[ids].reset_index(drop=True);remap={old:i for i,old in enumerate(ids)}
 values,mapping=read_calls(sites,meta,helper,out);meta.to_csv(out/'tables/display_specimens.tsv',sep='\t',index=False)
 pd.DataFrame(membership).to_csv(out/'tables/population_membership.tsv',sep='\t',index=False)
 pd.DataFrame(values,index=meta.specimen_id,columns=sites.pos).to_csv(out/'tables/individual_ALT_dosages.tsv.gz',sep='\t')
 newkeys=keys(sites);original=set(keys(previous));rows=[];labels=[];exports=[]
 for row,(r,g,indices) in enumerate(groups):
  x=values[[remap[i] for i in indices]];n=np.isfinite(x).sum(axis=0)
  freq=np.divide(np.nansum(x,axis=0),n,out=np.full(len(sites),np.nan),where=n>0)
  saved=dict(zip(keys(g),g.ALT_frequency));saved_n=dict(zip(keys(g),g.n_called))
  for j,k in enumerate(newkeys):
   if k in original:freq[j]=saved[k];n[j]=int(saved_n[k]) # retain the original plotted cells and denominators verbatim
  rows.append(freq);labels.append(helper.short(r.population)+f' | {r.morph} | n={int(r.n_total)}')
  for j,k in enumerate(newkeys):exports.append(dict(row=row,population=r.population,morph=r.morph,**dict(zip(KEYS,k)),n_total=len(x),n_called=int(n[j]),ALT_frequency=float(freq[j]) if np.isfinite(freq[j]) else None,previously_displayed=k in original))
 pd.DataFrame(exports).to_csv(out/'tables/population_frequencies.tsv.gz',sep='\t',index=False)
 lr,audit=base.focal(sites,'Dennis')
 frozen=pd.read_csv(src/'tables/corrected_LR_calls.tsv',sep='\t');frozen=frozen.loc[frozen.locus.eq('Dennis')]
 for _,g in frozen.groupby(KEYS):
  if len(g[['Left_ALT_dosage','Right_ALT_dosage']].drop_duplicates())!=1:raise RuntimeError('Conflicting old focal values')
 lookup={k:r for k,r in zip(keys(frozen.drop_duplicates(KEYS)),frozen.drop_duplicates(KEYS).to_dict('records'))}
 for j,k in enumerate(newkeys):
  if k in original:
   if k not in lookup:raise RuntimeError('Previous focal allele missing')
   lr[:,j]=[lookup[k]['Left_ALT_dosage'],lookup[k]['Right_ALT_dosage']]
 pd.DataFrame([dict(**dict(zip(KEYS,k)),Left_ALT_dosage=lr[0,j],Right_ALT_dosage=lr[1,j]) for j,k in enumerate(newkeys)]).to_csv(out/'tables/corrected_LR_calls.tsv',sep='\t',index=False)
 selected=set(newkeys);scores['displayed_in_genoplot']=[k in selected for k in keys(scores)]
 highlighted=scores.trait.isin(['D','R'])&scores.historical_highlight
 if not scores.loc[highlighted,'displayed_in_genoplot'].all():raise RuntimeError('Highlighted BA SNP lacks genoplot column')
 scores.to_csv(out/'tables/archived_BA_display_membership.tsv.gz',sep='\t',index=False)
 # Reuse the concise D/R renderer, with its colour overlay driven by the
 # original highlight flags, exactly as in the requested archived BA version.
 scores['displayed_in_genoplot']=scores.historical_highlight
 pdf=out/'Optix_DR_inclusive_genoplot.pdf';original_save=PdfPages.savefig;connector_counts=[]
 def finish(self,figure=None,**kwargs):
  lib.add_tracks(figure,scores)
  for ax in figure.axes:
   if ax.get_ylabel() in ['D','R']:
    ax.collections[0].set_facecolor('#C8CED6');ax.collections[0].set_alpha(.52)
    ax.collections[1].set_sizes([10]);ax.collections[1].set_alpha(.98)
  for t in figure.texts:
   if t.get_text().startswith('Genoplot:'):t.set_text(f'Genoplot: previous 176 SNPs + archived D or R BA ≥ 0.93 | {len(sites)} SNPs | no cap')
   elif t.get_text().startswith('Archived D/R BA'):t.set_text('D/R BA unchanged · coloured = historical BA > 0.99 + focal pattern · dashed = 0.8')
  old.add_phenotype_table(figure,'Dennis');connector_counts.append(lib.simplify_objects(figure))
  return original_save(self,figure,dpi=600,**kwargs)
 PdfPages.savefig=finish
 try:morph.render(pdf,sites,np.vstack([rows,lr]),labels+['Left (corrected)','Right (corrected)'],'Optix','Dennis',.93,len(sites),helper.SPANS['Dennis'],helper)
 finally:PdfPages.savefig=original_save
 save(out/'complete.json',dict(status='COMPLETE',PDF=str(pdf),SNPs=len(sites),previous_SNPs=len(previous),added_SNPs=len(sites)-len(previous),
  BA_recalculated=False,BA_colouring='unchanged historical highlights',all_highlighted_alleles_in_genoplot=True,previous_cells_unchanged=True,
  population_rows=len(labels),connectors=connector_counts,new_variant_calls=False,new_mixture_fit=False))
 print(f'PDF\t{pdf}\nSNPs\t{len(sites)}\nADDED\t{len(sites)-len(previous)}\nSTATUS\tCOMPLETE')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=SOURCE);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--submit',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
 if a.submit:prepare(a)
 elif a.run:run(a.output)
 else:p.error('Use --submit')
if __name__=='__main__':main()
