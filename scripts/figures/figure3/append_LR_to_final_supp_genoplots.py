#!/usr/bin/env python3
"""Append corrected L/R to the completed final supplemental genoplots, on Pegasus.

Frozen exported comparative frequencies, marker sets, order and renderers are
reused verbatim. No BA calculation, selection, grouping or variant calling.
"""
import argparse, hashlib, importlib.util, json, shlex, shutil, subprocess, sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
SOURCE=ROOT/'analysis/Figure3_for_Arnaud'
MORPH=SOURCE/'final_morph_supplements'
YBV=SOURCE/'Ybv_geography_supplements'
HO=ROOT/'analysis/LR_window_heterozygosity_v1'
OPTIX=ROOT/'analysis/rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc/optix_consistent_v2/reciprocal_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz'
OUTPUT=ROOT/'analysis/supp_genoplots_with_LR_20261027'
KEYS=['chrom','pos','ref','alt']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def load(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def prepare(a):
 if a.output.exists():raise RuntimeError('New output directory required; refusing to replace '+str(a.output))
 inputs=[MORPH/'Dennis_complete.json',MORPH/'worker.py',MORPH/'figure3_helper.py',
         YBV/'complete.json',YBV/'run.json',YBV/'worker.py',YBV/'figure3_helper.py',
         YBV/'all_regional_BA_and_denominators.tsv.gz',OPTIX,
         HO/'worker.py',HO/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz']
 for suffix in ['BA093','compact']:
  inputs += [MORPH/f'Dennis_{suffix}_{tail}' for tail in ['selected_SNPs.tsv','frequencies.tsv.gz','morph_frequencies.pdf']]
 for stem in ['Ybv_geography_27_original','Ybv_geography_regional_union']:
  inputs += [YBV/(stem+tail) for tail in ['_sites.tsv','_frequencies_and_callability.tsv.gz','.pdf']]
 for locus in ['Dennis','Ybv']:
  inputs += [SOURCE/'tables'/f'{locus}_focus_markers_and_updated_BA.tsv',SOURCE/'data'/locus/'corrected_focal_ALT_dosages.tsv']
 for p in inputs:
  if not p.is_file():raise FileNotFoundError('Missing frozen source: '+str(p))
 for p in [MORPH/'Dennis_complete.json',YBV/'complete.json']:
  if json.loads(p.read_text()).get('status')!='COMPLETE':raise RuntimeError('Source incomplete: '+str(p))
 out=a.output;out.mkdir(parents=True);(out/'logs').mkdir();(out/'tables').mkdir()
 shutil.copy2(__file__,out/'worker.py')
 save(out/'run.json',dict(workflow='frozen_supplement_genoplots_append_LR',hashes={str(p):sha(p) for p in inputs},
   only_change='Append corrected Left/Right rows at bottom, using the same ALT scale',
   comparative_values_changed=False,BA_recalculated=False,SNPs_reselected=False,new_calls=False))
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=8G','--cpus-per-task=1',
  '--job-name=Supp_LR','--output='+str(out/'logs/plot.%j.out'),'--error='+str(out/'logs/plot.%j.err'),
  '--wrap='+shlex.join([sys.executable,str(out/'worker.py'),'--run','--output',str(out)])],text=True).strip().split(';')[0]
 save(out/'jobs.json',dict(plot=job));print('JOB\t'+job+'\nOUTPUT\t'+str(out))

def focal(sites,locus):
 import numpy as np,pandas as pd
 normalizer=load(HO/'worker.py','frozen_ho_normalizer')
 p=OPTIX if locus=='Dennis' else HO/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz'
 calls=normalizer.normalize(pd.read_csv(p,sep='\t'),'Optix' if locus=='Dennis' else 'Ivory_Yb')
 # Exact REF/ALT identity is mandatory. No orientation guessed from position.
 lookup={tuple(r[k] for k in KEYS):r for r in calls.to_dict('records')}
 if len(lookup)!=len(calls):raise RuntimeError('Duplicate corrected allele identity')
 focus=pd.read_csv(SOURCE/'tables'/f'{locus}_focus_markers_and_updated_BA.tsv',sep='\t')
 frozen=pd.read_csv(SOURCE/'data'/locus/'corrected_focal_ALT_dosages.tsv',sep='\t',index_col=0)
 focus_lookup={tuple(r[k] for k in KEYS):i for i,r in enumerate(focus.to_dict('records'))}
 if list(map(int,frozen.columns))!=focus.pos.astype(int).tolist():raise RuntimeError('Focal cache order differs')
 values=np.full((2,len(sites)),np.nan);audit=[];gtmap={'0/0':0.,'0/1':.5,'1/0':.5,'1/1':1.}
 for j,r in enumerate(sites.to_dict('records')):
  k=tuple(r[c] for c in KEYS);source='no_exact_qualified_corrected_call'
  if k in focus_lookup:
   values[:,j]=frozen.iloc[:,focus_lookup[k]].to_numpy(float);source='unchanged_main_focus_corrected_cache'
  elif k in lookup:
   c=lookup[k]
   good=c['pair_posterior']>=.95 and c['left_depth']>=10 and c['right_depth']>=10 and bool(c['other_ok'])
   if good:
    values[:,j]=[gtmap.get(str(c['mixture_'+side+'_GT']).replace('|','/'),np.nan) for side in ['left','right']]
    source='existing_full_region_corrected_pair_P95_DP10'
  audit.append(dict(locus=locus,**{c:r[c] for c in KEYS},Left_ALT_dosage=float(values[0,j]) if np.isfinite(values[0,j]) else None,
    Right_ALT_dosage=float(values[1,j]) if np.isfinite(values[1,j]) else None,source=source))
 return values,audit

def frozen_rows(exports,sites,group_columns,label):
 import numpy as np
 rows=[];labels=[]
 for _,g in exports.groupby(group_columns,sort=False,dropna=False):
  # Exports are already in the saved figure's row order; preserve it.
  if len(g)!=len(sites) or g[KEYS].to_dict('records')!=sites[KEYS].to_dict('records'):
   raise RuntimeError('Saved frequency allele/order differs from saved selected sites')
  rows.append(g.ALT_frequency.to_numpy(float));labels.append(label(g.iloc[0]))
 return np.asarray(rows),labels

def run(out):
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 cfg=json.loads((out/'run.json').read_text())
 for p,h in cfg['hashes'].items():
  if sha(p)!=h:raise RuntimeError('Frozen input changed: '+p)
 results=[];audits=[]
 morph=load(MORPH/'worker.py','frozen_morph_renderer');helper=load(MORPH/'figure3_helper.py','frozen_morph_style')
 versions=json.loads((MORPH/'Dennis_complete.json').read_text())['versions']
 for suffix in ['BA093','compact']:
  stem='Dennis_'+suffix;sites=pd.read_csv(MORPH/(stem+'_selected_SNPs.tsv'),sep='\t')
  exp=pd.read_csv(MORPH/(stem+'_frequencies.tsv.gz'),sep='\t',keep_default_na=False,na_values=[''])
  rows,labels=frozen_rows(exp,sites,['population','morph','geography','phenotype'],
    lambda r:helper.short(r.population)+f' | {r.morph} | n={int(r.n_total)}')
  lr,audit=focal(sites,'Dennis');audits.extend(audit)
  v=next(v for v in versions if v['version']==suffix)
  path=out/(stem+'_morph_frequencies_with_LR.pdf')
  morph.render(path,sites,np.vstack([rows,lr]),labels+['Left (corrected)','Right (corrected)'],
   'Optix / Dennis','Dennis',v['threshold'],len(sites),helper.SPANS['Dennis'],helper)
  shutil.copy2(MORPH/(stem+'_selected_SNPs.tsv'),out/'tables'/(stem+'_selected_SNPs.tsv'))
  shutil.copy2(MORPH/(stem+'_frequencies.tsv.gz'),out/'tables'/(stem+'_frequencies.tsv.gz'))
  results.append(dict(PDF=str(path),SNPs=len(sites),population_rows=len(labels),focal_called_per_side=np.isfinite(lr).sum(axis=1).tolist()))
 yb=load(YBV/'worker.py','frozen_ybv_renderer');helper=load(YBV/'figure3_helper.py','frozen_ybv_style')
 scores=pd.read_csv(YBV/'all_regional_BA_and_denominators.tsv.gz',sep='\t');cutoff=json.loads((YBV/'run.json').read_text())['cutoff']
 for stem,title in [('Ybv_geography_27_original','Yb(v): original 27 SNPs by geography and morph'),
   ('Ybv_geography_regional_union',f'Yb(v): original markers + regional BA ≥{cutoff:.2f} union')]:
  sites=pd.read_csv(YBV/(stem+'_sites.tsv'),sep='\t');exp=pd.read_csv(YBV/(stem+'_frequencies_and_callability.tsv.gz'),sep='\t',keep_default_na=False,na_values=[''])
  def label(r):
   base=f'All {r.geography} v'+('+' if int(r.phenotype)==1 else '−') if r.row_type=='region_phenotype' else '    '+helper.short(r.population)+' | '+r.morph
   return base+f' | n={int(r.n_total)}'
  rows,labels=frozen_rows(exp,sites,['row'],label)
  blocks=[]
  for i,g in exp.groupby('row',sort=False):
   r=g.iloc[0]
   if r.row_type=='region_phenotype':blocks.append((int(i),int(i)+1,r.geography,int(r.phenotype)))
  lr,audit=focal(sites,'Ybv');audits.extend(audit)
  saved=(np.vstack([rows,lr]),labels+['Left (corrected)','Right (corrected)'],exp,blocks)
  yb.grouped_rows=lambda *args,saved=saved:saved
  path=out/(stem+'_with_LR.pdf');_,pages=yb.draw(path,sites,scores,None,None,helper,title,cutoff)
  for tail in ['_sites.tsv','_frequencies_and_callability.tsv.gz']:shutil.copy2(YBV/(stem+tail),out/'tables'/(stem+tail))
  results.append(dict(PDF=str(path),SNPs=len(sites),pages=pages,population_rows=len(labels),focal_called_per_side=np.isfinite(lr).sum(axis=1).tolist()))
 pd.DataFrame(audits).to_csv(out/'tables/added_corrected_LR_calls.tsv',sep='\t',index=False)
 save(out/'complete.json',dict(status='COMPLETE',outputs=results,comparative_tables='Byte-identical copies of frozen original exports',
   only_change='Two corrected focal rows at bottom of each genoplot',BA_recalculated=False,SNPs_reselected=False,new_calls=False))
 for r in results:print('PDF\t'+r['PDF']+'\tSNPs='+str(r['SNPs'])+'\tL/R called='+str(r['focal_called_per_side']))
 print('STATUS\tCOMPLETE')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--submit',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
 if a.submit:prepare(a)
 elif a.run:run(a.output)
 else:p.error('Use --submit')
if __name__=='__main__':main()
