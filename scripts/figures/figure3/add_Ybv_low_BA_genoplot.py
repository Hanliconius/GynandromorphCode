#!/usr/bin/env python3
"""Add symmetric low-BA sites to the completed reviewed Ybv genoplot."""
import argparse,csv,gzip,hashlib,importlib.util,json,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis')
SOURCE=ROOT/'Ybv_sorted_reviewed_BA_geography_v2_20261102'
OUTPUT=ROOT/'Ybv_sorted_reviewed_BA_high_low_20261102'
YBV=ROOT/'Figure3_for_Arnaud/Ybv_geography_supplements'
APPEND=ROOT/'supp_genoplots_LR_cydno_historical_BDR_20261027/append_helper.py'
KEYS=['chrom','pos','ref','alt']
def load(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def prepare(a):
 if a.output.exists():raise RuntimeError('Preserving existing output: '+str(a.output))
 base=load(APPEND,'focal_inputs');h=load(YBV/'figure3_helper.py','style_inputs')
 paths=[a.source/n for n in ['complete.json','run.json','worker.py','reviewed_specimens.tsv','tables/all_populations_membership.tsv',
  'tables/all_populations_ALT_dosages.tsv.gz','tables/reconciled_exact_allele_inventory.tsv',
  'tables/recalculated_regional_BA_all_SNPs.tsv.gz','tables/fixed_selected_SNPs_updated_BA.tsv']]
 paths += [YBV/'worker.py',YBV/'figure3_helper.py',YBV/'run.json',APPEND,base.HO/'worker.py',base.HO/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz',
  base.SOURCE/'tables/Ybv_focus_markers_and_updated_BA.tsv',base.SOURCE/'data/Ybv/corrected_focal_ALT_dosages.tsv',
  ROOT/'supp_genoplots_LR_cydno_historical_BDR_20261027/tables/corrected_LR_calls.tsv',ROOT/'Ybv_all_populations_20261102/worker.py']
 for p in paths:
  if not p.is_file():raise FileNotFoundError('Missing documented completed source: '+str(p))
 if json.loads((a.source/'complete.json').read_text()).get('status')!='COMPLETE':raise RuntimeError('Reviewed BA source incomplete')
 cutoff=a.cutoff if a.cutoff is not None else float(json.loads((YBV/'run.json').read_text())['cutoff'])
 if not .5<cutoff<=1:raise ValueError('Cutoff must be >0.5 and <=1')
 a.output.mkdir(parents=True);(a.output/'logs').mkdir();(a.output/'tables').mkdir()
 shutil.copy2(__file__,a.output/'worker.py')
 save(a.output/'run.json',dict(source=str(a.source),cutoff=cutoff,low_cutoff=round(1-cutoff,12),hashes={str(p):sha(p) for p in paths},
  selection='previous displayed SNPs union all current West/East fixed-direction BA >=t or <=1-t; no SNP cap',row_order=a.order,BA_recalculated=False,new_calls=False))
 cmd=[sys.executable,str(a.output/'worker.py'),'--run','--output',str(a.output)]
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:15:00','--mem=8G','--cpus-per-task=1','--job-name=Ybv_BAtails',
  '--output='+str(a.output/'logs/plot.%j.out'),'--error='+str(a.output/'logs/plot.%j.err'),'--wrap='+shlex.join(cmd)],text=True).strip().split(';')[0]
 save(a.output/'jobs.json',dict(plot=job));print(f'JOB\t{job}\nHIGH_BA\t{cutoff}\nLOW_BA\t{round(1-cutoff,12)}\nOUTPUT\t{a.output}')
def choose(scores,previous,cutoff):
 import numpy as np
 s=scores.copy();old=set(map(tuple,previous[KEYS].itertuples(index=False,name=None)))
 if s.duplicated(KEYS).any() or previous.duplicated(KEYS).any():raise RuntimeError('Duplicate exact alleles')
 present=set(map(tuple,s[KEYS].itertuples(index=False,name=None)))
 if not old<=present:raise RuntimeError('Previous displayed alleles absent from completed score inventory')
 s['previously_displayed']=[tuple(r) in old for r in s[KEYS].itertuples(index=False,name=None)]
 for region in ['West','East']:
  b=s[region+'_fixed_direction_BA'];s[region+'_high_tail']=b.ge(cutoff);s[region+'_low_tail']=b.le(round(1-cutoff,12))
  s[region+'_tail']=np.where(s[region+'_high_tail'],'high',np.where(s[region+'_low_tail'],'low',np.where(b.notna(),'middle','unscored')))
 mask=s.previously_displayed|s.West_high_tail|s.West_low_tail|s.East_high_tail|s.East_low_tail
 return s.loc[mask].sort_values(KEYS).reset_index(drop=True),s

def reorder_population_rows(rows,labels,exports,order):
 """Change display order only; preserve every frequency and specimen group."""
 import re
 if order=='geography-phenotype':return rows,labels,exports,None
 first=exports.groupby('row',sort=False).first()
 if set(first.index)!=set(range(len(rows))):raise RuntimeError('Unexpected population row identifiers')
 regions={'East':0,'West':1,'Silvaniforms':2,'Unresolved':3};states={'1':0,'0':1,'?':2}
 def keys(i):
  r=first.loc[i];m=re.search(r'v([01?])',str(r.morph));state=m.group(1) if m else '?'
  return (states[state],regions.get(r.geography,3),str(r.geography),str(r.population),str(r.morph))
 idx=sorted(range(len(rows)),key=keys)
 mapping={old:new for new,old in enumerate(idx)}
 exported=exports.copy();exported['row']=exported.row.map(mapping)
 exported=exported.sort_values(['row']+KEYS,kind='stable').reset_index(drop=True)
 blocks=[];last=None
 for j,i in enumerate(idx):
  group=keys(i)[:3]
  if group!=last:blocks.append((j,j+1,str(first.loc[i].geography),group[0]));last=group
 return rows[idx],[labels[i] for i in idx],exported,blocks

def run(out):
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 from matplotlib.backends.backend_pdf import PdfPages
 cfg=json.loads((out/'run.json').read_text());src=Path(cfg['source']);cut=cfg['cutoff'];low=cfg['low_cutoff']
 for p,h in cfg['hashes'].items():
  if sha(p)!=h:raise RuntimeError('Completed input changed: '+p)
 worker=load(src/'worker.py','reviewed_source');renderer=load(YBV/'worker.py','frozen_renderer');helper=load(YBV/'figure3_helper.py','frozen_helper')
 style=load(ROOT/'Ybv_all_populations_20261102/worker.py','frozen_style');focal_helper=load(APPEND,'existing_corrected_focal')
 scores=pd.read_csv(src/'tables/recalculated_regional_BA_all_SNPs.tsv.gz',sep='\t');previous=pd.read_csv(src/'tables/fixed_selected_SNPs_updated_BA.tsv',sep='\t')
 selected,all_scores=choose(scores,previous,cut)
 selected.to_csv(out/'tables/selected_high_low_BA_sites.tsv',sep='\t',index=False)
 all_scores.to_csv(out/'tables/all_BA_tail_classifications.tsv.gz',sep='\t',index=False)
 inventory=pd.read_csv(src/'tables/reconciled_exact_allele_inventory.tsv',sep='\t')
 if inventory.duplicated(KEYS).any():raise RuntimeError('Duplicate genotype inventory')
 with gzip.open(src/'tables/all_populations_ALT_dosages.tsv.gz','rt') as f:header=next(csv.reader(f,delimiter='\t'))[1:]
 if header!=inventory.pos.astype(str).tolist():raise RuntimeError('Saved genotype matrix column order differs from exact allele inventory')
 matrix=pd.read_csv(src/'tables/all_populations_ALT_dosages.tsv.gz',sep='\t',index_col=0)
 meta=pd.read_csv(src/'tables/all_populations_membership.tsv',sep='\t',dtype=str,keep_default_na=False)
 if matrix.index.duplicated().any() or meta.specimen_id.duplicated().any() or matrix.index.astype(str).tolist()!=meta.specimen_id.tolist():raise RuntimeError('Specimen order changed')
 lookup={tuple(r):j for j,r in enumerate(inventory[KEYS].itertuples(index=False,name=None))}
 indices=[lookup[tuple(r)] for r in selected[KEYS].itertuples(index=False,name=None)]
 values=matrix.to_numpy(float)[:,indices]
 rows,labels,exports,blocks=worker.sorted_rows(meta,values,selected,helper)
 rows,labels,exports,newblocks=reorder_population_rows(rows,labels,exports,cfg.get('row_order','geography-phenotype'))
 if newblocks is not None:blocks=newblocks
 lr,audit=focal_helper.focal(selected,'Ybv')
 # Preserve every previous corrected L/R cell, including explicit missingness.
 frozen=pd.read_csv(ROOT/'supp_genoplots_LR_cydno_historical_BDR_20261027/tables/corrected_LR_calls.tsv',sep='\t');frozen=frozen.loc[frozen.locus.eq('Ybv')]
 for _,g in frozen.groupby(KEYS):
  if len(g[['Left_ALT_dosage','Right_ALT_dosage']].drop_duplicates())!=1:raise RuntimeError('Conflicting previous focal calls')
 old={tuple(r[c] for c in KEYS):r for r in frozen.drop_duplicates(KEYS).to_dict('records')}
 for j,r in enumerate(selected.to_dict('records')):
  k=tuple(r[c] for c in KEYS)
  if r['previously_displayed']:
   if k not in old:raise RuntimeError('Previous focal allele absent')
   lr[:,j]=[old[k]['Left_ALT_dosage'],old[k]['Right_ALT_dosage']]
   audit[j].update(Left_ALT_dosage=float(lr[0,j]) if np.isfinite(lr[0,j]) else None,Right_ALT_dosage=float(lr[1,j]) if np.isfinite(lr[1,j]) else None,source='unchanged_previous_display_corrected_call')
 pd.DataFrame(audit).to_csv(out/'tables/corrected_LR_calls.tsv',sep='\t',index=False)
 exports.to_csv(out/'tables/population_frequencies_and_membership.tsv.gz',sep='\t',index=False)
 saved=np.vstack([rows,lr]),labels+['Left (corrected)','Right (corrected)'],exports,blocks
 renderer.grouped_rows=lambda *args:saved
 original=PdfPages.savefig;page=0
 def finish(self,figure=None,**kwargs):
  nonlocal page
  ss=selected.iloc[page*160:(page+1)*160];page+=1
  for region in ['West','East']:
   ax=next(a for a in figure.axes if a.get_ylabel()==region+' BA')
   if len(ax.collections)!=2:raise RuntimeError('Unexpected BA layer structure')
   ax.collections[1].remove()
   b=ss[region+'_fixed_direction_BA'];hi=b.ge(cut);lo=b.le(low);middle=b.notna()&~hi&~lo
   ax.scatter(ss.loc[hi,'pos'],b[hi],s=8,c='#1764AB',lw=0,zorder=4)
   ax.scatter(ss.loc[lo,'pos'],b[lo],s=8,c='#B6404F',lw=0,zorder=4)
   ax.axhline(low,c='#999999',ls='--',lw=.5)
   if region=='West':
    from matplotlib.lines import Line2D
    handles=[
     Line2D([],[],marker='o',ls='none',color='#1764AB',ms=4,label=f'Displayed SNP: BA >= {cut:g}'),
     Line2D([],[],marker='o',ls='none',color='#B6404F',ms=4,label=f'Displayed SNP: BA <= {low:g}'),
     Line2D([],[],marker='o',ls='none',color='#555555',alpha=.55,ms=3,label='All scoreable SNPs')]
    ax.legend(handles=handles,loc='lower left',bbox_to_anchor=(0,1.08),ncol=3,frameon=False,fontsize=6,handletextpad=.4,columnspacing=1,borderaxespad=0)
  style.remove_original_marker_squares(figure);style.add_phenotype_table(figure,'Ybv')
  for t in figure.texts:
   if 'West v+, West v−' in t.get_text():t.set_text(f'Rows: {cfg.get("row_order","geography-phenotype")}; no consensus rows. Previous SNPs retained plus BA >= {cut:g} or <= {low:g} in either region.\nBA points: blue=displayed SNP with regional BA >= {cut:g}; red=displayed SNP with regional BA <= {low:g}. Grey background=all scoreable SNPs, including retained intermediate scores. Same West-derived direction across regions, with labelled East fallback.\nHighlights refer to columns on this page, with chromosome/position/REF/ALT matched exactly. Missing BA is not zero. Colours do not encode REF/ALT or zygosity.\nExisting >=80% per-class callability gate and reviewed P cohorts retained. Genoplot remains 0=REF, 1=ALT; no allele flipping.\nNew L/R cells use existing qualified corrected calls only; unavailable calls remain grey. Descriptive association, not independent validation.')
  return original(self,figure,**kwargs)
 PdfPages.savefig=finish;pdf=out/'Ybv_all_populations_high_low_BA_genoplot.pdf'
 try:_,pages=renderer.draw(pdf,selected,all_scores,None,None,helper,'Yb(v): strong associations in both directions',cut)
 finally:PdfPages.savefig=original
 save(out/'complete.json',dict(status='COMPLETE',PDF=str(pdf),pages=pages,SNPs=len(selected),previous_SNPs=len(previous),added_SNPs=int((~selected.previously_displayed).sum()),low_tail_SNPs=int((selected.West_low_tail|selected.East_low_tail).sum()),high_cutoff=cut,low_cutoff=low,population_rows=len(rows),row_order=cfg.get('row_order','geography-phenotype'),BA_recalculated=False,new_variant_calls=False,new_mixture_fit=False,previous_corrected_LR_unchanged=True,LR_called_per_side=np.isfinite(lr).sum(axis=1).tolist()))
 print(f'PDF\t{pdf}\nSNPs\t{len(selected)}\nADDED\t{int((~selected.previously_displayed).sum())}\nSTATUS\tCOMPLETE')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=SOURCE);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--cutoff',type=float);p.add_argument('--order',choices=['geography-phenotype','phenotype-geography'],default='geography-phenotype');p.add_argument('--submit',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
 if a.submit:prepare(a)
 elif a.run:run(a.output)
 else:p.error('Use --submit')
if __name__=='__main__':main()
