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
OUTPUT=ROOT/'analysis/supp_genoplots_LR_cydno_historical_BDR_20261027'
KEYS=['chrom','pos','ref','alt']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def load(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def prepare(a):
 if a.output.exists():raise RuntimeError('Refusing to replace an existing output directory: '+str(a.output))
 repo=Path(__file__).resolve().parents[3]
 helper=load(MORPH/'figure3_helper.py','documented_source_helper')
 inputs=[MORPH/'worker.py',MORPH/'figure3_helper.py',MORPH/'Dennis_complete.json',
  YBV/'worker.py',YBV/'figure3_helper.py',YBV/'run.json',YBV/'complete.json',YBV/'membership.tsv',YBV/'all_regional_BA_and_denominators.tsv.gz',
  SOURCE/'reviewed_specimens.tsv',SOURCE/'run.json',
  repo/'scripts/figures/figure3/append_LR_to_final_supp_genoplots.py',helper.YB_ENGINE,
  *HISTORICAL_SCORES.values(),*HISTORICAL_MANIFESTS.values(),HISTORICAL_HIGH,HISTORICAL_COUNTS]
 for locus in ['Dennis','Ybv']:
  inputs += [SOURCE/'tables'/f'{locus}_focus_markers_and_updated_BA.tsv',SOURCE/'data'/locus/'corrected_focal_ALT_dosages.tsv']
 inputs += [helper.HO/'worker.py',helper.HO/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz',
  helper.RUN/'optix_consistent_v2/reciprocal_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz']
 for suffix in ['BA093','compact']:
  inputs += [MORPH/f'Dennis_{suffix}_{tail}' for tail in ['selected_SNPs.tsv','frequencies.tsv.gz','morph_frequencies.pdf']]
 for stem in ['Ybv_geography_27_original','Ybv_geography_regional_union']:
  inputs += [YBV/(stem+tail) for tail in ['_sites.tsv','_frequencies_and_callability.tsv.gz','.pdf']]
 for p in inputs:
  if not p.is_file():raise FileNotFoundError('Missing documented source: '+str(p))
 for p in [MORPH/'Dennis_complete.json',YBV/'complete.json']:
  if json.loads(p.read_text()).get('status')!='COMPLETE':raise RuntimeError('Source incomplete: '+str(p))
 vcfs=[]
 cfg_path=helper.BROAD/'run.json'
 if cfg_path.is_file():
  bc=json.loads(cfg_path.read_text());inputs.append(cfg_path)
  for i,span in enumerate(bc['shards']):
   p=helper.BROAD/'parts'/f'{i:02d}.vcf.gz'
   if p.is_file():vcfs.append(dict(path=str(p),start=span['start'],end=span['end'],priority=0))
 if helper.YB_FULL.is_file():vcfs.append(dict(path=str(helper.YB_FULL),start=1657079,end=1660072,priority=1))
 if not vcfs:raise RuntimeError('No existing Yb VCF source available for western context; no new calling submitted')
 inputs += [Path(r['path']) for r in vcfs]
 for p in [SOURCE/'data/Ybv/all_plot_membership.tsv',SOURCE/'data/Ybv/all_plot_focus_ALT_dosages.tsv.gz']:
  if p.is_file():inputs.append(p)
 out=a.output;out.mkdir(parents=True);(out/'logs').mkdir();(out/'tables').mkdir()
 for p,name in [(Path(__file__),'worker.py'),(repo/'scripts/figures/figure3/append_LR_to_final_supp_genoplots.py','append_helper.py')]:shutil.copy2(p,out/name)
 save(out/'run.json',dict(workflow='final_supp_LR_west_cydno_historical_BDR_BA',hashes={str(p):sha(p) for p in inputs},
  existing_VCFs=vcfs,historical_Optix_BA=True,BA_recalculated=False,
  historical_figure='Science_trials_1_2_pubregions_cartoons.comparison.pdf; Analyses 158/159',
  historical_code_commit='cbc9e087942d91776ecfa41ef3dfdd66c773051e',
  previous_outputs_modified=False,SNPs_reselected=False,new_calls=False))
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=8G','--cpus-per-task=1',
  '--job-name=Supp_BDR','--output='+str(out/'logs/plot.%j.out'),'--error='+str(out/'logs/plot.%j.err'),
  '--wrap='+shlex.join([sys.executable,str(out/'worker.py'),'--run','--output',str(out)])],text=True).strip().split(';')[0]
 save(out/'jobs.json',dict(plot=job));print('JOB\t'+job+'\nOUTPUT\t'+str(out))

def optix_BA(out,cfg):
 """The actual Analysis-158/159 archived score clouds and highlight flags."""
 import numpy as np,pandas as pd
 # Keep the author's standing numata/pardalinus exclusion auditable even when
 # reusing an old aggregate score; never silently display a forbidden cohort.
 for label,path in HISTORICAL_MANIFESTS.items():
  members=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False)
  if 'cohort' in members:members=members.loc[members.cohort.eq('WGS_plus_CS')]
  if 'keep' in members:members=members.loc[members['keep'].str.lower().isin(['true','t','1','yes','y'])]
  taxon=next((c for c in ['taxon_name','taxon','Taxon Name'] if c in members),None)
  if members.empty or taxon is None:raise RuntimeError('Cannot validate archived comparative membership: '+str(path))
  if members[taxon].str.contains('numata|pardalinus',case=False,regex=True).any():raise RuntimeError('Archived BA retains numata/pardalinus; refusing to substitute this analysis or change its cohort silently')
  members.to_csv(out/'tables'/('September_'+label+'_retained_membership.tsv'),sep='\t',index=False)
 frames=[]
 for trait,path in HISTORICAL_SCORES.items():
  raw=pd.read_csv(path,sep='\t',keep_default_na=False,na_values=['','NA','NaN'])
  needed=set(KEYS+['analysis_cohort','scheme','balanced_accuracy'])
  if not needed<=set(raw.columns):raise RuntimeError('Archived score table lacks required schema: '+str(path))
  # Match 115b/142 and then 158 normalise_track verbatim: WGS_plus_CS,
  # DP2_GQ20, finite BA, fixed span; do NOT newly gate score_eligible.
  raw=raw.loc[raw.analysis_cohort.eq('WGS_plus_CS')&raw.scheme.eq('DP2_GQ20')].copy()
  raw['pos']=pd.to_numeric(raw.pos,errors='coerce');raw['balanced_accuracy']=pd.to_numeric(raw.balanced_accuracy,errors='coerce')
  if trait=='B':raw=raw.drop_duplicates(KEYS,keep='first')
  raw=raw.loc[np.isfinite(raw.pos)&np.isfinite(raw.balanced_accuracy)&raw.pos.between(1060000,1240000)].copy()
  raw['pos']=raw.pos.astype(int)
  if raw.empty or raw.duplicated(['chrom','pos']).any():raise RuntimeError('Empty or ambiguous archived '+trait+' BA track')
  if not raw.chrom.astype(str).isin(['Hmr18','18','NC_140406.1']).all():raise RuntimeError('Historical BA is not on the documented Rosina/Hmr18 chromosome')
  if not raw.balanced_accuracy.between(0,1).all():raise RuntimeError('Invalid archived BA')
  raw['trait']=trait;raw['source']=str(path);frames.append(raw)
 scores=pd.concat(frames,ignore_index=True)
 high=pd.read_csv(HISTORICAL_HIGH,sep='\t')
 if not {'trait','pos','BA'}<=set(high):raise RuntimeError('Archived Analysis-158 highlight export has unexpected columns')
 high['pos']=pd.to_numeric(high.pos,errors='raise').astype(int)
 if high.duplicated(['trait','pos']).any() or not high.BA.gt(.99).all() or not high.trait.isin(['B','D','R']).all():raise RuntimeError('Invalid historical focal-compatible highlight set')
 # These counts identify the screenshot, not a new SNP-selection target.
 observed=dict(B=int(high.loc[high.trait.eq('B'),'pos'].nunique()),DR_union=int(high.loc[high.trait.isin(['D','R']),'pos'].nunique()))
 if observed!={'B':10,'DR_union':45}:raise RuntimeError('Archived highlights differ from the screenshot (expected B=10, D/R union=45): '+str(observed))
 counts=pd.read_csv(HISTORICAL_COUNTS,sep='\t').set_index('category')['n']
 for category in ['Band focal-compatible','D/R focal-compatible']:
  if category not in counts or int(counts[category])!=367:raise RuntimeError('Archived count audit does not match the screenshot: '+category)
 lookup={(r.trait,int(r.pos)):float(r.balanced_accuracy) for r in scores.itertuples()}
 for r in high.itertuples():
  if (r.trait,r.pos) not in lookup or not np.isclose(lookup[(r.trait,r.pos)],r.BA,atol=1e-10,rtol=0):raise RuntimeError('Historical high-point BA disagrees with archived cloud at '+str((r.trait,r.pos)))
 highlights=set(zip(high.trait,high.pos))
 scores['historical_highlight']=[(r.trait,int(r.pos)) in highlights for r in scores.itertuples()]
 scores.to_csv(out/'tables/Optix_September_B_D_R_BA.tsv.gz',sep='\t',index=False)
 shutil.copy2(HISTORICAL_HIGH,out/'tables'/HISTORICAL_HIGH.name);shutil.copy2(HISTORICAL_COUNTS,out/'tables'/HISTORICAL_COUNTS.name)
 provenance=dict(source='Exact archived Analysis-158/159 B/D/R BA; no recalculation or modern-cohort substitution',
  score_paths={k:str(v) for k,v in HISTORICAL_SCORES.items()},filter='WGS_plus_CS; DP2_GQ20; finite scores; 1,060,000–1,240,000',
  historical_highlights='BA strictly >0.99 and historical focal-compatible genotype pattern; archived flags, not a current validation set',
  expected_screenshot_counts=observed,modern_genoplot_unchanged=True,
  historical_cohorts_are_not_current_genoplot_cohorts=True,BA_recalculated=False)
 save(out/'tables/Optix_September_BDR_BA_provenance.json',provenance)
 return scores,provenance['source']


def add_BA_to_figure(fig,sites,scores):
 """Reinstate the three quantitative tracks from the supplied Analysis-159 image."""
 import numpy as np
 from matplotlib.patches import Rectangle
 w,h=fig.get_size_inches();newh=h+4.;fig.set_size_inches(w,newh)
 for ax in list(fig.axes):
  p=ax.get_position();ax.set_position([p.x0,p.y0*h/newh,p.width,p.height*h/newh])
 for text in list(fig.texts):
  x,y=text.get_position()
  text.set_position((x,(newh-.25)/newh if y>.97 else (newh-.6)/newh if y>.94 else y*h/newh))
  if .94<y<.97:text.set_text('Genoplot: '+text.get_text())
 palette={'B':'#0072B2','D':'#E69F00','R':'#CC79A7'}
 for trait,bottom in [('B',h+1.65),('D',h+.80),('R',h-.05)]:
  s=scores.loc[scores.trait.eq(trait)];ax=fig.add_axes([.33,bottom/newh,.65,.62/newh])
  # Original track transform clamps BA below .5 to the baseline.
  y=np.clip(s.balanced_accuracy.to_numpy(float),.5,1.)
  ax.scatter(s.pos,y,s=2.2,c='#C8CED6',alpha=.52,lw=0)
  selected=s.historical_highlight.to_numpy(bool)
  ax.scatter(s.loc[selected,'pos'],y[selected],s=10,c=palette[trait],alpha=.98,lw=0)
  ax.axhline(.5,c='#AEB6BF',lw=.5);ax.axhline(.8,c='#666666',ls=(0,(2,2)),lw=.6)
  ax.set(xlim=(1060000,1240000),ylim=(.495,1.015),yticks=[],xticks=[])
  ax.text(-.013,.50,trait,transform=ax.transAxes,ha='right',va='center',weight='bold',fontsize=10,color='#242424')
  for sp in ax.spines.values():sp.set_visible(False)
 ann=fig.add_axes([.33,(h+2.50)/newh,.65,.40/newh]);ann.set(xlim=(1060000,1240000),ylim=(0,1.5),yticks=[],xticks=[])
 ann.add_patch(Rectangle((1088602,.98),1097846-1088602,.18,facecolor='#3A7CA5',lw=0));ann.text((1088602+1097846)/2,1.22,'optix',ha='center',fontsize=8,color='#3A7CA5',style='italic')
 for name,a,b,trait,y in [('B',1112651,1114314,'B',.40),('R1',1156746,1160220,'R',.40),('R2',1174547,1187833,'R',.40),('D1',1215210,1216124,'D',.40),('D2',1218557,1221494,'D',.10)]:
  ann.add_patch(Rectangle((a,y),b-a,.16,facecolor=palette[trait],lw=0));ann.text((a+b)/2,y+.19,name,ha='center',fontsize=7,color=palette[trait])
 for sp in ann.spines.values():sp.set_visible(False)
 fig.text(.33,(h+3.04)/newh,'September archived B/D/R BA | scale 0.5–1.0; dashed = 0.8 | coloured: historical focal-compatible BA >0.99.',fontsize=7)
 fig.text(.33,(h+2.29)/newh,'Historical BA cohorts/highlights are separate from the current genoplot selection; genotype rows and SNP columns below are unchanged.',fontsize=6.5)


def west_assignment(r):
 """Only reviewed West or specific documented western cydno/pachinus localities."""
 region=str(r.get('geography_used_Fig3b','')).strip()
 if region=='East':return ''
 if region=='West':return 'reviewed West assignment'
 country=str(r.get('country','')).lower()
 if 'panam' in country or 'costa rica' in country:return 'recorded Panama/Costa Rica locality'
 try:longitude=float(r.get('longitude',''))
 except (TypeError,ValueError):return ''
 if 'colombia' in country and longitude<=-75:return 'recorded western Colombian locality (longitude <= -75)'
 return ''

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
 taxa=reviewed.taxon.str.contains(r'Heliconius\s+(?:cydno|pachinus)\b',case=False,regex=True)
 meta=reviewed.loc[taxa&~reviewed.globally_excluded&reviewed.Yb_ventral.isin(['0','1'])].copy()
 meta['western_evidence']=[west_assignment(r) for r in meta.to_dict('records')];meta=meta.loc[meta.western_evidence.ne('')].reset_index(drop=True)
 if meta.empty:raise RuntimeError('No independently documented West cydno/pachinus metadata')
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
 meta['morph']=['d'+str(r.Yb_dorsal or '?')+' v'+str(r.Yb_ventral) for r in meta.itertuples()]
 meta['qualified_sites']=np.isfinite(values).sum(axis=1);meta['context_only']=True
 original=pd.read_csv(YBV/'membership.tsv',sep='\t',dtype=str,keep_default_na=False)
 meta['already_in_frozen_display']=meta.specimen_id.isin(original.specimen_id)
 meta.to_csv(out/'tables/western_cydno_pachinus_membership.tsv',sep='\t',index=False)
 pd.DataFrame(values,index=meta.specimen_id,columns=sites.pos).to_csv(out/'tables/western_cydno_pachinus_ALT_dosages.tsv.gz',sep='\t')
 audit=[]
 for i,r in meta.iterrows():
  for j,site in enumerate(sites.to_dict('records')):audit.append(dict(specimen_id=r.specimen_id,**{c:site[c] for c in KEYS},source=sources[i,j]))
 pd.DataFrame(audit).to_csv(out/'tables/western_cydno_call_sources.tsv.gz',sep='\t',index=False)
 if not np.isfinite(values).any():raise RuntimeError('West candidates found but no existing qualified exact-allele calls; no new calling attempted')
 return meta,values


def context_rows(meta,values,sites,helper):
 import numpy as np,pandas as pd
 usable=np.isfinite(values).any(axis=1);meta=meta.loc[usable].copy();values=values[usable];meta=meta.reset_index(drop=True)
 rows=[];labels=[];exports=[]
 for (state,taxon,morph),g in meta.groupby(['Yb_ventral','display_taxon','morph'],sort=True):
  x=values[g.index];n=np.isfinite(x).sum(axis=0);v=np.divide(np.nansum(x,axis=0),n,out=np.full(len(sites),np.nan),where=n>0)
  rows.append(v);labels.append('West context | '+helper.short(taxon)+' | '+morph+f' | n={len(g)}')
  for j,r in enumerate(sites.to_dict('records')):exports.append(dict(taxon=taxon,morph=morph,phenotype=int(state),**{k:r[k] for k in KEYS},n_total=len(g),n_called=int(n[j]),ALT_frequency=float(v[j]) if n[j] else None,specimen_ids=';'.join(g.specimen_id),context_only=True))
 return np.asarray(rows).reshape((-1,len(sites))),labels,pd.DataFrame(exports)


def run(out):
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 from matplotlib.backends.backend_pdf import PdfPages
 cfg=json.loads((out/'run.json').read_text())
 for p,h in cfg['hashes'].items():
  if sha(p)!=h:raise RuntimeError('Input changed: '+p)
 base=load(out/'append_helper.py','frozen_LR_functions');morph=load(MORPH/'worker.py','frozen_morph_renderer');helper=load(MORPH/'figure3_helper.py','frozen_morph_helper')
 scores,BAorigin=optix_BA(out,cfg);outputs=[];focal_audits=[];phenotype_audit=[]
 versions=json.loads((MORPH/'Dennis_complete.json').read_text())['versions']
 for suffix in ['BA093','compact']:
  stem='Dennis_'+suffix;sites=pd.read_csv(MORPH/(stem+'_selected_SNPs.tsv'),sep='\t');exp=pd.read_csv(MORPH/(stem+'_frequencies.tsv.gz'),sep='\t',keep_default_na=False,na_values=[''])
  rows,labels=base.frozen_rows(exp,sites,['population','morph','geography','phenotype'],lambda r:helper.short(r.population)+f' | {r.morph} | n={int(r.n_total)}')
  lr,audit=base.focal(sites,'Dennis');focal_audits.extend(audit)
  original_save=PdfPages.savefig
  def save_with_BA(self,figure=None,**kwargs):
   add_BA_to_figure(figure,sites,scores)
   phenotype_audit.extend(dict(PDF=path.name,**r) for r in add_phenotype_table(figure,'Dennis'))
   return original_save(self,figure,**kwargs)
  PdfPages.savefig=save_with_BA
  path=out/(stem+'_LR_September_BDR_BA.pdf');version=next(v for v in versions if v['version']==suffix)
  try:morph.render(path,sites,np.vstack([rows,lr]),labels+['Left (corrected)','Right (corrected)'],'Optix: archived B/D/R association and current morph genotypes','Dennis',version['threshold'],len(sites),helper.SPANS['Dennis'],helper)
  finally:PdfPages.savefig=original_save
  for tail in ['_selected_SNPs.tsv','_frequencies.tsv.gz']:shutil.copy2(MORPH/(stem+tail),out/'tables'/(stem+tail))
  outputs.append(dict(PDF=str(path),SNPs=len(sites),BA_tracks=['Band','Dennis','Ray'],historical_September_BA=True,original_rows_unchanged=True))
 yb=load(YBV/'worker.py','frozen_ybv_renderer');helper=load(YBV/'figure3_helper.py','frozen_ybv_helper')
 regional=pd.read_csv(YBV/'all_regional_BA_and_denominators.tsv.gz',sep='\t');cutoff=json.loads((YBV/'run.json').read_text())['cutoff']
 union=pd.read_csv(YBV/'Ybv_geography_regional_union_sites.tsv',sep='\t');meta,calls=cydno_calls(out,cfg,union,helper);lookup={tuple(r[c] for c in KEYS):j for j,r in enumerate(union.to_dict('records'))}
 for stem,title in [('Ybv_geography_27_original','Yb(v): original 27 SNPs by geography and morph'),('Ybv_geography_regional_union',f'Yb(v): original markers + regional BA ≥{cutoff:.2f} union')]:
  sites=pd.read_csv(YBV/(stem+'_sites.tsv'),sep='\t');exp=pd.read_csv(YBV/(stem+'_frequencies_and_callability.tsv.gz'),sep='\t',keep_default_na=False,na_values=[''])
  def label(r):
   text=f'All {r.geography} v'+('+' if int(r.phenotype)==1 else '−') if r.row_type=='region_phenotype' else '    '+helper.short(r.population)+' | '+r.morph
   return text+f' | n={int(r.n_total)}'
  rows,labels=base.frozen_rows(exp,sites,['row'],label);blocks=[]
  for i,g in exp.groupby('row',sort=False):
   r=g.iloc[0]
   if r.row_type=='region_phenotype':blocks.append((int(i),int(i)+1,r.geography,int(r.phenotype)))
  idx=[lookup[tuple(r[c] for c in KEYS)] for r in sites.to_dict('records')];added,names,context=context_rows(meta,calls[:,idx],sites,helper)
  if not len(added):raise RuntimeError('No West context rows have qualified calls in '+stem)
  blocks.append((len(rows),len(rows)+len(added),'West context',0))
  lr,audit=base.focal(sites,'Ybv');focal_audits.extend(audit)
  saved=np.vstack([rows,added,lr]),labels+names+['Left (corrected)','Right (corrected)'],exp,blocks
  yb.grouped_rows=lambda *args,saved=saved:saved
  original_save=PdfPages.savefig
  def save_context_note(self,figure=None,**kwargs):
   remove_original_marker_squares(figure)
   phenotype_audit.extend(dict(PDF=path.name,**r) for r in add_phenotype_table(figure,'Ybv'))
   for text in figure.texts:
    if 'West v+' in text.get_text() and 'n = specimens' in text.get_text():text.set_text(text.get_text().replace('Black squares mark original sites. ','')+'\nWest cydno/pachinus context rows are display-only subsets; frozen regional summaries and BA cohorts are unchanged. Corrected L/R use the same ALT scale.')
   return original_save(self,figure,**kwargs)
  PdfPages.savefig=save_context_note
  path=out/(stem+'_LR_western_cydno_pachinus.pdf')
  try:_,pages=yb.draw(path,sites,regional,None,None,helper,title,cutoff)
  finally:PdfPages.savefig=original_save
  for tail in ['_sites.tsv','_frequencies_and_callability.tsv.gz']:shutil.copy2(YBV/(stem+tail),out/'tables'/(stem+tail))
  context.to_csv(out/'tables'/(stem+'_added_western_context.tsv.gz'),sep='\t',index=False)
  outputs.append(dict(PDF=str(path),SNPs=len(sites),pages=pages,added_context_rows=names,original_rows_unchanged=True,Ybv_BA_changed=False))
 pd.DataFrame(focal_audits).to_csv(out/'tables/corrected_LR_calls.tsv',sep='\t',index=False)
 pd.DataFrame(phenotype_audit).to_csv(out/'tables/phenotype_display_codes.tsv',sep='\t',index=False)
 (out/'README.txt').write_text('Final supplementary genoplots: additions only\n\nTwo Optix versions retain all current comparative frequency rows and SNP columns, with corrected Left/Right at the bottom and the archived September Analysis-158/159 B/D/R BA tracks above. Historical cohorts and coloured focal-compatible BA>0.99 flags are reused without recalculation; they are not the current genotype display cohort or an independent validation.\nTwo Yb geography versions retain the original regional BA, summaries and SNP columns, add independently documented western cydno/pachinus context rows, and retain corrected Left/Right at the bottom. These added rows are display-only subsets, not new BA individuals; overlaps with the original summaries are exported. Original-marker black squares have been removed; connector lines remain.\nPhenotype boxes: filled=present, open=absent, ?=unknown, minus=pooled/mixed; focal phenotype boxes are blank unless supplied. Codes and all added per-site denominators are exported. Original comparative exports are byte-identical copies. No download, alignment, variant calling or SNP reselection.\n\nOptix BA provenance: '+BAorigin+'\n')
 save(out/'complete.json',dict(status='COMPLETE',outputs=outputs,Optix_BA_source=BAorigin,original_SNP_selection_unchanged=True,original_population_rows_unchanged=True,
  Western_context_in_main_BA=False,Western_context_in_regional_BA=False,phenotype_boxes=True,Yb_connector_squares_removed=True,BA_recalculated=False,new_calls=False,downloads=False,alignments=False))
 for r in outputs:print('PDF\t'+r['PDF']+'\tSNPs='+str(r['SNPs']))
 print('STATUS\tCOMPLETE')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--submit',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
 if a.submit:prepare(a)
 elif a.run:run(a.output)
 else:p.error('Use --submit')
if __name__=='__main__':main()
