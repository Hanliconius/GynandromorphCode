#!/usr/bin/env python3
"""Compare BA estimands using exactly the retained lower-panel Ybv calls."""
import argparse,getpass,hashlib,importlib.util,json,shlex,shutil,subprocess,sys
from pathlib import Path
DEFAULT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/Figure3_for_Arnaud')
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def load(p):
 spec=importlib.util.spec_from_file_location('frozen_ba',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def prepare(a):
 source=a.source.resolve();out=source/'Ybv_BA_rule_comparison'
 paths=[source/'complete.json',source/'data/Ybv/membership.tsv',source/'data/Ybv/focus_ALT_dosages.tsv.gz',source/'tables/Ybv_focus_markers_and_updated_BA.tsv',source/'code/geography_helper.py']
 for p in paths:
  if not p.is_file():raise RuntimeError('Missing completed input: '+str(p))
 if json.loads(paths[0].read_text()).get('status')!='COMPLETE':raise RuntimeError('Source package is not complete')
 if out.exists():
  if not a.overwrite:raise RuntimeError('Use --overwrite for the existing diagnostic')
  if json.loads((out/'run.json').read_text()).get('workflow')!='Ybv_focus_BA_rules':raise RuntimeError('Output belongs to another workflow')
  jobs=json.loads((out/'jobs.json').read_text()) if (out/'jobs.json').exists() else {}
  active=set(subprocess.check_output(['squeue','-h','-u',getpass.getuser(),'-o','%A'],text=True).split())
  if active & set(jobs.values()):raise RuntimeError('Diagnostic job is still active')
 out.mkdir(exist_ok=True);(out/'logs').mkdir(exist_ok=True)
 (out/'complete.json').unlink(missing_ok=True)
 shutil.copy2(__file__,out/'worker.py');shutil.copy2(paths[-1],out/'geography_helper.py')
 save(out/'run.json',dict(workflow='Ybv_focus_BA_rules',source=str(source),input_sha256={str(p):digest(p) for p in paths},script_sha256=digest(Path(__file__)),
  rule='Pooled best allele orientation per rule, held fixed for West; >=80% callability per phenotype class; pooled ties unscored',
  calls='Exact cached 27-SNP lower-panel matrix; no new calls, selection, metadata edits or main-figure changes'))
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:10:00','--cpus-per-task=1','--mem=4G','--job-name=Ybv_BArules',
  '--output='+str(out/'logs/compare.%j.out'),'--error='+str(out/'logs/compare.%j.err'),
  '--wrap='+shlex.join([sys.executable,str(out/'worker.py'),'--run','--output',str(out)])],text=True).strip().split(';')[0]
 save(out/'jobs.json',dict(compare=job));print('JOB\t'+job+'\nOUTPUT\t'+str(out))

def run(out):
 import numpy as np,pandas as pd
 cfg=json.loads((out/'run.json').read_text());source=Path(cfg['source'])
 for path,expected in cfg['input_sha256'].items():
  if digest(Path(path))!=expected:raise RuntimeError('Input changed after submission: '+path)
 m=pd.read_csv(source/'data/Ybv/membership.tsv',sep='\t',dtype={'specimen_id':str})
 f=pd.read_csv(source/'tables/Ybv_focus_markers_and_updated_BA.tsv',sep='\t')
 x=pd.read_csv(source/'data/Ybv/focus_ALT_dosages.tsv.gz',sep='\t',index_col=0)
 if len(f)!=27 or m.specimen_id.duplicated().any() or f.pos.duplicated().any():raise RuntimeError('Expected 27 unique focus sites and independent specimens')
 if list(map(int,x.columns))!=f.pos.astype(int).tolist() or set(x.index)!=set(m.specimen_id):raise RuntimeError('Matrix/membership/marker identity mismatch')
 if m.astype(str).agg(' | '.join,axis=1).str.contains('numata|pardalinus',case=False).any():raise RuntimeError('Excluded taxon in comparison')
 if not m.phenotype.isin([0,1]).all() or not m.geography_used_Fig3b.isin(['East','West']).all():raise RuntimeError('Unresolved phenotype/geography')
 matrix=x.loc[m.specimen_id].to_numpy(float)
 if not np.isin(matrix[np.isfinite(matrix)],[0,.5,1]).all():raise RuntimeError('Unexpected dosage')
 y=m.phenotype.to_numpy(int);geo=m.geography_used_Fig3b.to_numpy();helper=load(out/'geography_helper.py')
 wide=f[['chrom','pos','ref','alt']].copy();records=[];summaries=[]
 for rule in ['carrier','homozygote']:
  raw,results,orientation,_=helper.calculate(matrix,y,geo,rule)
  comparisons=['pooled'] if rule=='carrier' else ['pooled','West_only']
  for comparison in comparisons:
   name=('pooled_' if comparison=='pooled' else 'western_')+rule
   score=results[comparison];r=raw[comparison]
   wide[name+'_BA']=score;wide[name+'_BA90']=np.isfinite(score)&(score>=.9);wide[name+'_orientation']=orientation
   for j,site in enumerate(f.itertuples()):
    records.append(dict(comparison=name,chrom=site.chrom,pos=int(site.pos),ref=site.ref,alt=site.alt,BA=None if not np.isfinite(score[j]) else float(score[j]),
     orientation=orientation[j],n_total_positive=int(r['n_total_1']),n_total_negative=int(r['n_total_0']),n_called_positive=int(r['n_called_1'][j]),n_called_negative=int(r['n_called_0'][j]),
     callable_fraction_positive=float(r['rate_1'][j]),callable_fraction_negative=float(r['rate_0'][j]),
     ALT_predicts_positive_BA=None if not np.isfinite(r['ALT'][j]) else float(r['ALT'][j]),REF_predicts_positive_BA=None if not np.isfinite(r['REF'][j]) else float(r['REF'][j])))
   summaries.append(dict(comparison=name,positive_n=int(r['n_total_1']),negative_n=int(r['n_total_0']),scoreable_SNPs=int(np.isfinite(score).sum()),BA_ge_0_90=int((np.isfinite(score)&(score>=.9)).sum())))
 wide['previous_western_homozygote_BA']=f.balanced_accuracy.to_numpy(float)
 a=wide.western_homozygote_BA.to_numpy();b=wide.previous_western_homozygote_BA.to_numpy()
 wide['western_score_matches_previous']=np.isclose(a,b,equal_nan=True,atol=1e-10,rtol=0)
 wide['pooled_homozygote_minus_western']=wide.pooled_homozygote_BA-wide.western_homozygote_BA
 wide['pooled_carrier_minus_homozygote']=wide.pooled_carrier_BA-wide.pooled_homozygote_BA
 wide.to_csv(out/'Ybv_27_SNP_BA_comparison.tsv',sep='\t',index=False)
 pd.DataFrame(records).to_csv(out/'Ybv_BA_denominators_and_orientations.tsv',sep='\t',index=False)
 pd.DataFrame(summaries).to_csv(out/'summary.tsv',sep='\t',index=False)
 m.to_csv(out/'membership.tsv',sep='\t',index=False)
 plot(out,wide,summaries)
 (out/'README.txt').write_text('All three analyses use the exact cached lower-panel 27-SNP genotypes and final reviewed specimen identities. Pooled means all eligible regions, not only the eight example rows. BA averages sensitivity and specificity, giving equal weight to the two phenotype classes. Each site requires >=80% callability in both classes. Carrier: at least one copy of the associated allele predicts v+. Homozygote: two copies predict v+. Missing genotypes are excluded, not treated as homozygous REF. For each rule, allele direction is chosen from pooled BA and held fixed within West; ties are unscored. Per-site class denominators and both allele orientations are exported. Prior western BA is included to detect differences between cached focus and broad-region calls. No main plots or scientific inputs are modified. High BA here is descriptive in-sample association, not an independent validation.\n')
 save(out/'complete.json',dict(status='COMPLETE',comparisons=summaries,SNPs=27,western_previous_score_mismatches=int((~wide.western_score_matches_previous).sum()),new_calls=False,main_figures_changed=False))
 print(pd.DataFrame(summaries).to_string(index=False));print('WESTERN_PREVIOUS_SCORE_MISMATCHES\t'+str(int((~wide.western_score_matches_previous).sum())))
 print('TABLE\t'+str(out/'Ybv_27_SNP_BA_comparison.tsv')+'\nPDF\t'+str(out/'Ybv_27_SNP_BA_comparison.pdf'))

def plot(out,wide,summary):
 import numpy as np,matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.colors import Normalize
 cols=['pooled_carrier','pooled_homozygote','western_homozygote'];x=wide[[c+'_BA' for c in cols]].to_numpy(float)
 fig,ax=plt.subplots(figsize=(8,10));fig.subplots_adjust(left=.25,right=.88,bottom=.17,top=.86)
 cmap=plt.get_cmap('viridis').copy();cmap.set_bad('#DDDDDD');im=ax.imshow(np.ma.masked_invalid(x),cmap=cmap,norm=Normalize(0,1),aspect='auto')
 ax.set_xticks(range(3),['Pooled\ncarrier','Pooled\nhomozygote','Western\nhomozygote']);ax.xaxis.tick_top()
 ax.set_yticks(range(len(wide)),[f'{r.pos:,} {r.ref}>{r.alt}' for r in wide.itertuples()],fontsize=8)
 for i in range(len(wide)):
  for j in range(3):
   v=x[i,j];ax.text(j,i,'NA' if not np.isfinite(v) else f'{v:.3f}'+(' *' if v>=.9 else ''),ha='center',va='center',fontsize=8,color='black' if not np.isfinite(v) or v>.6 else 'white')
 fig.colorbar(im,ax=ax,fraction=.035,pad=.035,label='Balanced accuracy')
 fig.suptitle('Yb(v): the same 27 focus SNPs under three comparisons',fontsize=12)
 fig.text(.06,.07,'* BA >=0.90. Grey = unscored (callability or pooled orientation tie).\n'+'; '.join(r['comparison']+f": v+ {r['positive_n']}, v- {r['negative_n']}" for r in summary)+'\nAllele orientation fixed from pooled score separately for each genotype rule.',fontsize=8,wrap=True)
 plt.rcParams['pdf.fonttype']=42;fig.savefig(out/'Ybv_27_SNP_BA_comparison.pdf');plt.close(fig)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=DEFAULT);p.add_argument('--submit',action='store_true');p.add_argument('--overwrite',action='store_true');p.add_argument('--run',action='store_true');p.add_argument('--output',type=Path);a=p.parse_args()
 if a.submit:prepare(a)
 elif a.run and a.output:run(a.output)
 else:p.error('Use --submit or --run --output PATH')
if __name__=='__main__':main()
