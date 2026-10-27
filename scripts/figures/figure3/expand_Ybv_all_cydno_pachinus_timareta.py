#!/usr/bin/env python3
"""New final supplements: corrected L/R, West cydno/pachinus, historical B/D/R BA.

Uses frozen comparative exports and existing qualified calls. No SNP reselection,
downloads, alignment, variant calling, or changes to previous output packages.
"""
import argparse,hashlib,importlib.util,json,re,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
SOURCE=ROOT/'analysis/Figure3_for_Arnaud'
MORPH=SOURCE/'final_morph_supplements'
YBV=SOURCE/'Ybv_geography_supplements'
HISTORICAL_SCORES={
 'B':ROOT/'analysis/B_western_postman_vs_western_Amazon_DR_reanalysis/checkpoints/B_western_BA_and_adjusted_LR.tsv',
 'D':ROOT/'analysis/drop_problematic_reanalysis/checkpoints/longDennis_original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.tsv',
 'R':ROOT/'analysis/drop_problematic_reanalysis/checkpoints/longRay_original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.tsv'}
HISTORICAL_DIR=ROOT/'analysis/figures/final_B_D_R_combined/science_condensation_trials'
HISTORICAL_HIGH=HISTORICAL_DIR/'Science_trials_high_correct_BDR_BAgt099_sites.tsv'
HISTORICAL_COUNTS=HISTORICAL_DIR/'Science_trials_site_count_audit.tsv'
HISTORICAL_MANIFESTS={
 'B':ROOT/'analysis/B_western_postman_vs_western_Amazon_DR_reanalysis/checkpoints/B_western_retained_sample_manifest.tsv',
 'DR':ROOT/'analysis/drop_problematic_reanalysis/checkpoints/original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv'}
PREVIOUS=ROOT/'analysis/supp_genoplots_LR_cydno_historical_BDR_20261027'
OUTPUT=ROOT/'analysis/Ybv_all_cydno_pachinus_timareta_20261027'
TAXA=r'\b(?:cydno|pachinus|timareta)\b'
KEYS=['chrom','pos','ref','alt']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def load(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def prepare(a):
 if a.output.exists():raise RuntimeError('Preserving existing output: '+str(a.output))
 helper=load(YBV/'figure3_helper.py','frozen_helper')
 prior=json.loads((PREVIOUS/'run.json').read_text())
 inputs=[PREVIOUS/'complete.json',PREVIOUS/'run.json',PREVIOUS/'append_helper.py',
  PREVIOUS/'tables/corrected_LR_calls.tsv',SOURCE/'reviewed_specimens.tsv',
  YBV/'worker.py',YBV/'figure3_helper.py',YBV/'complete.json',YBV/'run.json',YBV/'membership.tsv',
  YBV/'all_regional_BA_and_denominators.tsv.gz',helper.YB_ENGINE,
  SOURCE/'data/Ybv/all_plot_focus_ALT_dosages.tsv.gz',SOURCE/'tables/Ybv_focus_markers_and_updated_BA.tsv']
 for stem in ['Ybv_geography_27_original','Ybv_geography_regional_union']:
  inputs += [YBV/(stem+tail) for tail in ['_sites.tsv','_frequencies_and_callability.tsv.gz']]
 inputs += [Path(r['path']) for r in prior['existing_VCFs']]
 for p in inputs:
  if not p.is_file():raise FileNotFoundError('Missing documented input: '+str(p))
 if json.loads((PREVIOUS/'complete.json').read_text()).get('status')!='COMPLETE':raise RuntimeError('Previous supplement incomplete')
 a.output.mkdir(parents=True);(a.output/'logs').mkdir();(a.output/'tables').mkdir()
 shutil.copy2(__file__,a.output/'worker.py')
 save(a.output/'run.json',dict(hashes={str(p):sha(p) for p in inputs},existing_VCFs=prior['existing_VCFs'],
  all_geographies=True,unknown_phenotypes_displayed=True,BA_recalculated=False))
 cmd=[sys.executable,str(a.output/'worker.py'),'--run','--output',str(a.output)]
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=8G','--cpus-per-task=1','--job-name=Ybv_alltaxa',
  '--output='+str(a.output/'logs/plot.%j.out'),'--error='+str(a.output/'logs/plot.%j.err'),'--wrap='+shlex.join(cmd)],text=True).strip().split(';')[0]
 save(a.output/'jobs.json',dict(plot=job));print('JOB\t'+job+'\nOUTPUT\t'+str(a.output))

def phenotype_codes(label,locus):
 traits=['B','D','R'] if locus=='Dennis' else ['d','v']
 if label.strip().startswith(('Left (corrected)','Right (corrected)')):return {t:None for t in traits}
 result={}
 for trait in traits:
  hit=re.search(r'(?<!\w)'+trait+r'([01?])(?=\W|$)',label)
  if hit:result[trait]=hit.group(1);continue
  sign=re.search(r'(?<!\w)v([+\-−])',label) if trait=='v' else None
  result[trait]=('1' if sign.group(1)=='+' else '0') if sign else ('mixed' if label.strip().startswith('All ') else '?')
 return result

def add_phenotype_table(fig,locus):
 """Insert a narrow phenotype table without changing the genotype panel width."""
 from matplotlib.patches import Rectangle
 heat=fig.axes[0];labels=[x.get_text() for x in heat.get_yticklabels()]
 traits=['B','D','R'] if locus=='Dennis' else ['d','v']
 oldw,h=fig.get_size_inches();extra=.95;neww=oldw+extra
 fig.set_size_inches(neww,h)
 old_axes=list(fig.axes)
 oldleft=heat.get_position().x0*oldw
 for ax in old_axes:
  p=ax.get_position();ax.set_position([(p.x0*oldw+extra)/neww,p.y0,p.width*oldw/neww,p.height])
 for text in fig.texts:
  x,y=text.get_position();text.set_position((x*oldw/neww,y))
 # Row labels remain where they were before the additional table was inserted.
 oldpad=heat.yaxis.majorTicks[0].get_pad() if heat.yaxis.majorTicks else 6
 heat.tick_params(axis='y',pad=oldpad+extra*72)
 hp=heat.get_position();table=fig.add_axes([(oldleft+.06)/neww,hp.y0,.80/neww,hp.height])
 table.set(xlim=(-.5,len(traits)-.5),ylim=heat.get_ylim(),xticks=[],yticks=[])
 for sp in table.spines.values():sp.set_visible(False)
 table.patch.set_alpha(0)
 for j,t in enumerate(traits):table.text(j,1.015,t if locus=='Dennis' else 'Yb('+t+')',transform=table.get_xaxis_transform(),ha='center',va='bottom',fontsize=7,weight='bold')
 audit=[]
 for i,label in enumerate(labels):
  codes=phenotype_codes(label,locus);audit.append(dict(locus=locus,row=i,label=label,**codes))
  for j,t in enumerate(traits):
   value=codes[t]
   if value is None:continue
   fill='#333333' if value=='1' else '#FFFFFF' if value=='0' else '#D4D4D4'
   table.add_patch(Rectangle((j-.34,i-.31),.68,.62,facecolor=fill,edgecolor='#555555',lw=.45))
   if value in ['?','mixed']:table.text(j,i,'?' if value=='?' else '−',ha='center',va='center',fontsize=6,color='#333333')
 fig.text(.025,.007,'Phenotype boxes: filled = present; open = absent; ? = unknown; − = pooled/mixed. Focal rows remain unscored where no phenotype code was supplied.',fontsize=6)
 return audit

def remove_original_marker_squares(fig):
 import numpy as np
 heat=fig.axes[0];removed=0
 for line in list(heat.lines):
  if line.get_marker()=='s' and np.all(np.asarray(line.get_ydata(),float)==-.5):line.remove();removed+=1
 return removed


def cydno_calls(out,cfg,sites,helper):
 import numpy as np,pandas as pd,pysam
 reviewed=helper.read_metadata(SOURCE/'reviewed_specimens.tsv')
 taxa=reviewed.taxon.str.contains(TAXA,case=False,regex=True)
 meta=reviewed.loc[taxa&~reviewed.globally_excluded].copy()
 meta=meta.reset_index(drop=True)
 meta['geography_display']=meta.geography_used_Fig3b.map(lambda x:x if x in ['West','East'] else 'Unresolved')
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
 meta['display_taxon']=meta.taxon.map(lambda s:' '.join(re.sub(r'\b(?:ssp|spp|subsp)\.?(?=\s|$)','',str(s),flags=re.I).split()))
 meta['morph']=['d'+str(r.Yb_dorsal or '?')+' v'+str(r.Yb_ventral or '?') for r in meta.itertuples()]
 meta['qualified_sites']=np.isfinite(values).sum(axis=1);meta['context_only']=True
 original=pd.read_csv(YBV/'membership.tsv',sep='\t',dtype=str,keep_default_na=False)
 meta['already_in_frozen_display']=meta.specimen_id.isin(original.specimen_id)
 meta.to_csv(out/'tables/all_cydno_pachinus_timareta_membership.tsv',sep='\t',index=False)
 pd.DataFrame(values,index=meta.specimen_id,columns=sites.pos).to_csv(out/'tables/all_cydno_pachinus_timareta_ALT_dosages.tsv.gz',sep='\t')
 audit=[]
 for i,r in meta.iterrows():
  for j,site in enumerate(sites.to_dict('records')):audit.append(dict(specimen_id=r.specimen_id,**{c:site[c] for c in KEYS},source=sources[i,j]))
 pd.DataFrame(audit).to_csv(out/'tables/all_target_call_sources.tsv.gz',sep='\t',index=False)
 if not np.isfinite(values).any():raise RuntimeError('Target taxa found but no existing qualified exact-allele calls; no new calling attempted')
 return meta,values


def context_rows(meta,values,sites,helper):
 import numpy as np,pandas as pd
 usable=np.isfinite(values).any(axis=1);meta=meta.loc[usable].copy();values=values[usable];meta=meta.reset_index(drop=True)
 rows=[];labels=[];exports=[]
 for (geography,taxon,morph,state),g in meta.groupby(['geography_display','display_taxon','morph','Yb_ventral'],sort=True):
  x=values[g.index];n=np.isfinite(x).sum(axis=0);v=np.divide(np.nansum(x,axis=0),n,out=np.full(len(sites),np.nan),where=n>0)
  rows.append(v);labels.append(geography+' | '+helper.short(taxon)+' | '+morph+f' | n={len(g)}')
  for j,r in enumerate(sites.to_dict('records')):exports.append(dict(taxon=taxon,morph=morph,geography=geography,phenotype=int(state) if state in ['0','1'] else None,**{k:r[k] for k in KEYS},n_total=len(g),n_called=int(n[j]),ALT_frequency=float(v[j]) if n[j] else None,specimen_ids=';'.join(g.specimen_id),context_only=True))
 return np.asarray(rows).reshape((-1,len(sites))),labels,pd.DataFrame(exports)


def run(out):
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 from matplotlib.backends.backend_pdf import PdfPages
 cfg=json.loads((out/'run.json').read_text())
 for p,h in cfg['hashes'].items():
  if sha(p)!=h:raise RuntimeError('Input changed: '+p)
 base=load(PREVIOUS/'append_helper.py','frozen_rows')
 yb=load(YBV/'worker.py','frozen_yb_renderer');helper=load(YBV/'figure3_helper.py','frozen_helper')
 regional=pd.read_csv(YBV/'all_regional_BA_and_denominators.tsv.gz',sep='\t');cutoff=json.loads((YBV/'run.json').read_text())['cutoff']
 union=pd.read_csv(YBV/'Ybv_geography_regional_union_sites.tsv',sep='\t')
 meta,calls=cydno_calls(out,cfg,union,helper)
 lookup={tuple(r[c] for c in KEYS):j for j,r in enumerate(union.to_dict('records'))}
 saved_focal=pd.read_csv(PREVIOUS/'tables/corrected_LR_calls.tsv',sep='\t');saved_focal=saved_focal.loc[saved_focal.locus.eq('Ybv')]
 for _,g in saved_focal.groupby(KEYS):
  if len(g[['Left_ALT_dosage','Right_ALT_dosage']].drop_duplicates())!=1:raise RuntimeError('Inconsistent saved focal calls')
 saved_focal=saved_focal.drop_duplicates(KEYS)
 outputs=[]
 for stem in ['Ybv_geography_27_original','Ybv_geography_regional_union']:
  sites=pd.read_csv(YBV/(stem+'_sites.tsv'),sep='\t');exp=pd.read_csv(YBV/(stem+'_frequencies_and_callability.tsv.gz'),sep='\t',keep_default_na=False,na_values=[''])
  # Replace, rather than duplicate, target population rows. Frozen regional
  # aggregate rows and BA still describe their original comparative cohorts.
  target_rows=exp.row_type.eq('population_morph') & exp.population.fillna('').str.contains(TAXA,case=False,regex=True)
  removed=exp.loc[target_rows].copy();exp=exp.loc[~target_rows].copy()
  def label(r):
   text=f'All {r.geography} v'+('+' if int(r.phenotype)==1 else '−') if r.row_type=='region_phenotype' else '    '+helper.short(r.population)+' | '+r.morph
   return text+f' | n={int(r.n_total)}'
  rows,labels=base.frozen_rows(exp,sites,['row'],label);blocks=[]
  for i,(_,g) in enumerate(exp.groupby('row',sort=False)):
   if g.iloc[0].row_type=='region_phenotype':blocks.append((i,i+1,g.iloc[0].geography,int(g.iloc[0].phenotype)))
  idx=[lookup[tuple(r[c] for c in KEYS)] for r in sites.to_dict('records')]
  added,names,context=context_rows(meta,calls[:,idx],sites,helper)
  if not len(added):raise RuntimeError('No qualified target calls for '+stem)
  joined=sites[KEYS].merge(saved_focal,on=KEYS,how='left',validate='one_to_one',indicator=True)
  if not joined['_merge'].eq('both').all():raise RuntimeError('Missing saved corrected focal allele')
  lr=joined[['Left_ALT_dosage','Right_ALT_dosage']].to_numpy(float).T
  blocks.append((len(rows),len(rows)+len(added),'All target taxa',0))
  saved=np.vstack([rows,added,lr]),labels+names+['Left (corrected)','Right (corrected)'],exp,blocks
  yb.grouped_rows=lambda *args,saved=saved:saved
  original=PdfPages.savefig
  def finish(self,figure=None,**kwargs):
   remove_original_marker_squares(figure);add_phenotype_table(figure,'Ybv')
   figure.text(.025,.025,'Added taxa: all recorded regions; ? = unresolved phenotype. Regional summaries and BA retain their original cohorts.',fontsize=6)
   return original(self,figure,**kwargs)
  PdfPages.savefig=finish
  path=out/(stem+'_all_cydno_pachinus_timareta_LR.pdf')
  try:_,pages=yb.draw(path,sites,regional,None,None,helper,'Yb(v): all cydno, pachinus and timareta',cutoff)
  finally:PdfPages.savefig=original
  context.to_csv(out/'tables'/(stem+'_target_frequencies.tsv.gz'),sep='\t',index=False)
  removed.to_csv(out/'tables'/(stem+'_replaced_population_rows.tsv.gz'),sep='\t',index=False)
  audit=meta.copy();audit['qualified_sites_this_panel']=np.isfinite(calls[:,idx]).sum(axis=1)
  audit['displayed_this_panel']=audit.qualified_sites_this_panel.gt(0)
  audit['not_displayed_reason']=np.where(audit.displayed_this_panel,'','no qualified exact-allele calls in selected SNPs')
  audit.to_csv(out/'tables'/(stem+'_specimen_accounting.tsv'),sep='\t',index=False)
  outputs.append(dict(PDF=str(path),pages=pages,SNPs=len(sites),target_specimens=int(audit.displayed_this_panel.sum()),target_population_rows=len(names),unavailable_specimens=int((~audit.displayed_this_panel).sum())))
 (out/'README.txt').write_text('All reviewed non-excluded cydno, pachinus and timareta with available qualified exact-allele calls are displayed, including unresolved phenotypes and geography. ssp/spp spellings are normalized before pooling. Existing target population rows are replaced to avoid duplication. Other population rows, regional summaries, BA, selected SNPs and corrected L/R remain frozen. New target context is display-only; regional summaries are not sums of all displayed context rows. Specimen accounting lists zero-call individuals. No downloads, alignments or new calls.\n')
 save(out/'complete.json',dict(status='COMPLETE',outputs=outputs,BA_changed=False,new_calls=False,downloads=False,all_geographies=True))
 for r in outputs:print('PDF\t'+r['PDF']+'\tTARGET_SPECIMENS='+str(r['target_specimens'])+'\tUNAVAILABLE='+str(r['unavailable_specimens']))
 print('STATUS\tCOMPLETE')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--submit',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
 if a.submit:prepare(a)
 elif a.run:run(a.output)
 else:p.error('Use --submit')
if __name__=='__main__':main()
