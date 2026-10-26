#!/usr/bin/env python3
"""Compact morph-frequency supplements from frozen reviewed Figure 3 calls."""
import argparse,hashlib,importlib.util,json,shlex,shutil,subprocess,sys,zipfile
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
SOURCE=ROOT/'analysis/Figure3_for_Arnaud'
THRESHOLDS=[.90,.91,.92,.93,.94,.95,.96,.97,.98,.99,1.]

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def load(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def submit(a):
 out=a.output;out.mkdir(parents=True,exist_ok=True);(out/'logs').mkdir(exist_ok=True)
 if (out/'complete.json').exists() and not a.overwrite:raise RuntimeError('Use --overwrite for completed output')
 for p in [a.source/'complete.json',a.source/'reviewed_specimens.tsv',a.source/'code/worker.py']:
  if not p.is_file():raise FileNotFoundError(p)
 for locus in ['Dennis','Ybv']:
  for n in ['complete.json','membership.tsv','all_SNP_BA.tsv.gz','all_ALT_dosages.tsv.gz']:
   p=a.source/'data'/locus/n
   if not p.is_file():raise FileNotFoundError(p)
 shutil.copy2(__file__,out/'worker.py');shutil.copy2(a.source/'code/worker.py',out/'figure3_helper.py')
 save(out/'run.json',dict(source=str(a.source),output=str(out),target_columns=80,
      metadata_sha256=sha(a.source/'reviewed_specimens.tsv'),worker_sha256=sha(out/'worker.py')))
 (out/'complete.json').unlink(missing_ok=True)
 worker=[sys.executable,str(out/'worker.py'),'--output',str(out),'--stage','locus']
 array=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=8G',
  '--cpus-per-task=1','--array=0-1','--job-name=Fig3_morph','--output='+str(out/'logs/morph.%A_%a.out'),
  '--error='+str(out/'logs/morph.%A_%a.err'),'--wrap='+shlex.join(worker)],text=True).strip().split(';')[0]
 report=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=2G',
  '--job-name=Fig3_pack','--dependency=afterok:'+array,'--output='+str(out/'logs/report.%j.out'),
  '--error='+str(out/'logs/report.%j.err'),'--wrap='+shlex.join([sys.executable,str(out/'worker.py'),'--output',str(out),'--stage','report'])],text=True).strip().split(';')[0]
 save(out/'jobs.json',dict(locus_array=array,report=report))
 print('MORPH_ARRAY\t'+array+'\nREPORT\t'+report+'\nOUTPUT\t'+str(out))

def matrix(path,meta,sites):
 import pandas as pd
 d=pd.read_csv(path,sep='\t',index_col=0);d.index=d.index.astype(str)
 if d.shape!=(len(meta),len(sites)):raise RuntimeError('Matrix/site/membership shape mismatch: '+str(path))
 if set(d.index)!=set(meta.specimen_id):raise RuntimeError('Specimen membership differs from cached matrix')
 return d.loc[meta.specimen_id].to_numpy(float)

def yb_display(helper,src,scores,selection,selected_matrix):
 import pandas as pd,numpy as np,re
 reviewed=helper.read_metadata(src/'reviewed_specimens.tsv')
 sm=pd.read_csv(helper.BROAD/'membership.tsv',sep='\t').fillna('')
 sites=pd.read_csv(helper.BROAD/'SNPs.tsv',sep='\t')
 raw=pd.read_csv(helper.BROAD/'individual_ALT_dosages.tsv.gz',sep='\t',index_col=0)
 records=[(r['sample'],set(re.split(r'[;, ]+',r.get('aliases','')))|{r['sample']}) for r in sm.to_dict('records')]
 available=set().union(*(ids for _,ids in records))
 # Retain BA-excluded geographical exceptions for display when broad calls exist.
 # PKs stay outside these comparative population pools; global exclusions still apply.
 text=reviewed.population.astype(str)+' '+reviewed.taxon.astype(str)
 pk=text.str.contains(r'piano|PK_dark|PK_pale|Butterfly\s+world',case=False,regex=True)
 meta=reviewed.loc[~reviewed.globally_excluded & ~pk & reviewed.Yb_ventral.isin(['0','1'])].copy()
 meta=meta.loc[[bool(helper.aliases(r)&available) for r in meta.to_dict('records')]].reset_index(drop=True)
 meta['phenotype']=meta.Yb_ventral.astype(int)
 meta['geography_used_Fig3b']=meta.geography_used_Fig3b.replace('','Unresolved')
 mapping=helper.map_sources(meta,records)
 lookup={helper.key(r.chrom,r.pos,r.ref,r.alt):i for i,r in enumerate(sites.itertuples())}
 focus=pd.read_csv(src/'tables/Ybv_focus_markers_and_updated_BA.tsv',sep='\t')
 cached=pd.read_csv(helper.YF/'individual_ALT_dosages.tsv.gz',sep='\t',index_col=0)
 focuslookup={helper.key(r.chrom,r.pos,r.ref,r.alt):i for i,r in enumerate(focus.itertuples())}
 if list(map(int,cached.columns))!=focus.pos.astype(int).tolist():raise RuntimeError('Yb focus cache order differs')
 if raw.shape[1]!=len(sites):raise RuntimeError('Broad Yb cache/site order differs')
 vectors=[];conflicts=[]
 for i,r in meta.iterrows():
  samples=mapping[i]
  broad,n=helper.union(raw.loc[samples].to_numpy(float))
  small_samples=[s for s in samples if s in cached.index]
  small,ns=helper.union(cached.loc[small_samples].to_numpy(float)) if small_samples else (None,0)
  row=[]
  for s in scores.itertuples():
   k=helper.key(s.chrom,s.pos,s.ref,s.alt)
   if k in focuslookup and small is not None:row.append(small[focuslookup[k]])
   elif k in lookup:row.append(broad[lookup[k]])
   else:raise RuntimeError('Display source lacks selected Yb allele: '+str(k))
  vectors.append(row);conflicts.append(dict(specimen_id=r.specimen_id,broad_conflicts=n,focus_conflicts=ns))
 values=np.asarray(vectors)
 # Ensure the displayed selection members use exactly the frozen main-BA calls.
 idx={s:i for i,s in enumerate(meta.specimen_id)}
 for i,s in enumerate(selection.specimen_id):
  if s not in idx:raise RuntimeError('BA specimen absent from broader display')
  current=values[idx[s]];frozen=selected_matrix[i];both=np.isfinite(current)&np.isfinite(frozen)
  if not np.allclose(current[both],frozen[both],atol=1e-12):raise RuntimeError('Display/frozen BA calls disagree: '+s)
  values[idx[s]]=frozen
 return meta,values,conflicts

def render(path,sites,rows,labels,title,locus,threshold,n_pass,span,helper):
 import numpy as np,matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.backends.backend_pdf import PdfPages
 from matplotlib.colors import LinearSegmentedColormap,Normalize
 from matplotlib.patches import ConnectionPatch,Rectangle
 cmap=LinearSegmentedColormap.from_list('previous_palette',helper.SUPP,N=257);cmap.set_bad('#C8C8C8')
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42})
 # Pagination is exhaustive, not a cap or rank-thinned subset.
 with PdfPages(path) as pdf:
  for first in range(0,len(sites),80):
   ss=sites.iloc[first:first+80];v=rows[:,first:first+80];nr=len(labels)
   fig=plt.figure(figsize=(17,max(6,3.3+.22*nr)))
   heat=fig.add_axes([.33,.17,.65,.62]);ann=fig.add_axes([.33,.87,.65,.055])
   heat.pcolormesh(np.arange(len(ss)+1)-.5,np.arange(nr+1)-.5,np.ma.masked_invalid(v),cmap=cmap,vmin=0,vmax=1,rasterized=False)
   heat.set(xlim=(-.5,len(ss)-.5),ylim=(nr-.5,-.5),xticks=[])
   heat.set_yticks(np.arange(nr),labels,fontsize=7);heat.tick_params(length=0,pad=6)
   lo,hi=span[1:];ann.set(xlim=(lo,hi),ylim=(0,1.5),yticks=[])
   ann.axhline(.15,color='#777777',lw=.6)
   intervals=[('B',1112651,1114314),('R1',1156746,1160220),('R2',1174547,1187833),('D1',1215210,1216124),('D2',1218557,1221494)] if locus=='Dennis' else [('cortex',1642343,1653300),('Yb(v)',1658079,1659072),('ivory*',1659008,1741475)]
   for j,(label,a,b) in enumerate(intervals):
    a,b=max(a,lo),min(b,hi)
    if a<b:
     y=.4+(j%2)*.4;ann.add_patch(Rectangle((a,y),b-a,.13,color='#777777'));ann.text((a+b)/2,y+.2,label,ha='center',fontsize=7)
   ann.set_xticks(np.linspace(lo,hi,5),[f'{x/1e6:.3f}' for x in np.linspace(lo,hi,5)],fontsize=7)
   ann.set_xlabel('Rosina chromosome '+('18' if locus=='Dennis' else '15')+' (Mb)',fontsize=7)
   for j,r in enumerate(ss.itertuples()):
    fig.add_artist(ConnectionPatch((r.pos,.15),(j,-.5),coordsA=ann.transData,coordsB=heat.transData,lw=.3,color='#BBBBBB',alpha=.5))
   for ax in [heat,ann]:
    for sp in ax.spines.values():sp.set_visible(False)
   fig.text(.03,.98,title,fontsize=12,weight='bold',va='top')
   fig.text(.03,.95,f'BA >= {threshold:.2f} | {n_pass} passing SNPs, all retained | columns {first+1}-{first+len(ss)} | no site cap',fontsize=8)
   cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=cmap),cax=fig.add_axes([.53,.105,.4,.014]),orientation='horizontal',ticks=[0,.5,1])
   cb.set_label('ALT allele frequency (mean individual genotype dosage / 2)',fontsize=8)
   fig.text(.03,.035,'One row per population/morph/phenotype; Yb additionally split by recorded geography. n = specimens; per-SNP n_called is exported.\nPale blue / gold / purple is a frequency scale: the middle colour does not imply that every specimen is heterozygous.\nBA is an in-sample descriptive classifier score; selection membership/rule are frozen from the main figure. Ivory annotation is partial.',fontsize=7)
   pdf.savefig(fig);plt.close(fig)

def locus(out):
 import os,numpy as np,pandas as pd
 cfg=json.loads((out/'run.json').read_text());src=Path(cfg['source'])
 if sha(src/'reviewed_specimens.tsv')!=cfg['metadata_sha256']:raise RuntimeError('Reviewed metadata changed')
 helper=load(out/'figure3_helper.py','final_package');name=['Dennis','Ybv'][int(os.environ['SLURM_ARRAY_TASK_ID'])]
 folder=src/'data'/name;done=json.loads((folder/'complete.json').read_text())
 if done['status']!='COMPLETE':raise RuntimeError('Main statistics incomplete')
 scores=pd.read_csv(folder/'all_SNP_BA.tsv.gz',sep='\t')
 keys=['chrom','pos','ref','alt']
 if scores.duplicated(keys).any():raise RuntimeError('Duplicate allele identity')
 selection=pd.read_csv(folder/'membership.tsv',sep='\t',dtype=str,keep_default_na=False)
 m=matrix(folder/'all_ALT_dosages.tsv.gz',selection,scores)
 if name=='Dennis':meta=selection;display=m;conflicts=[]
 else:meta,display,conflicts=yb_display(helper,src,scores,selection,m)
 if any(helper.blocked(r) for r in meta.to_dict('records')):raise RuntimeError('Excluded taxon in display')
 scores['_index']=np.arange(len(scores));scores=scores.sort_values(keys)
 countrows=[dict(locus=name,threshold=t,passing_SNPs=int(scores.balanced_accuracy.ge(t).sum())) for t in THRESHOLDS]
 pd.DataFrame(countrows).to_csv(out/(name+'_threshold_counts.tsv'),sep='\t',index=False)
 viable=[r for r in countrows if 0<r['passing_SNPs']<=cfg['target_columns']]
 nonempty=[r for r in countrows if r['passing_SNPs']>0]
 if not nonempty:raise RuntimeError('No SNPs pass BA >= 0.90')
 compact=viable[0]['threshold'] if viable else nonempty[-1]['threshold']
 membership=meta.copy();membership['used_in_main_BA']=membership.specimen_id.isin(selection.specimen_id)
 membership.to_csv(out/(name+'_display_membership.tsv'),sep='\t',index=False)
 pd.DataFrame(conflicts).to_csv(out/(name+'_source_conflicts.tsv'),sep='\t',index=False)
 results=[]
 for suffix,cut in [('BA093',.93),('compact',compact)]:
  sites=scores.loc[scores.balanced_accuracy.ge(cut)].copy();sites.to_csv(out/f'{name}_{suffix}_selected_SNPs.tsv',sep='\t',index=False)
  if sites.empty:
   results.append(dict(version=suffix,threshold=cut,SNPs=0,PDF=None));continue
  values=display[:,sites['_index'].to_numpy(int)];rows=[];labels=[];exports=[]
  morph_traits=['B','D','R'] if name=='Dennis' else ['Yb_dorsal','Yb_ventral']
  meta['morph']=[' '.join(c.replace('Yb_dorsal','d').replace('Yb_ventral','v')+str(r.get(c,'?')) for c in morph_traits) for r in meta.to_dict('records')]
  cols=['population','morph','phenotype'] if name=='Dennis' else ['geography_used_Fig3b','population','morph','phenotype']
  # Matrix rows are aligned to the unmodified meta row order.
  meta=meta.reset_index(drop=True)
  for group,g in meta.groupby(cols,sort=True,dropna=False):
   x=values[g.index];n=np.isfinite(x).sum(axis=0)
   freq=np.divide(np.nansum(x,axis=0),n,out=np.full(len(sites),np.nan),where=n>0)
   rows.append(freq);geo='' if name=='Dennis' else str(group[0])+' | '
   pop=group[0] if name=='Dennis' else group[1];state=int(float(group[-1]));trait='D' if name=='Dennis' else 'v'
   labels.append(geo+helper.short(pop)+f' | {group[-2]} | n={len(g)}')
   for j,r in enumerate(sites.itertuples()):
    exports.append(dict(population=pop,morph=group[-2],geography='' if name=='Dennis' else group[0],phenotype=state,
     chrom=r.chrom,pos=int(r.pos),ref=r.ref,alt=r.alt,BA=float(r.balanced_accuracy),n_total=len(g),n_called=int(n[j]),
     ALT_frequency=None if n[j]==0 else float(freq[j]),
     REF_REF_fraction=None if n[j]==0 else float(np.sum(x[:,j]==0)/n[j]),
     REF_ALT_fraction=None if n[j]==0 else float(np.sum(x[:,j]==.5)/n[j]),
     ALT_ALT_fraction=None if n[j]==0 else float(np.sum(x[:,j]==1)/n[j])))
  pdf=out/f'{name}_{suffix}_morph_frequencies.pdf'
  render(pdf,sites,np.asarray(rows),labels,'Optix / Dennis' if name=='Dennis' else 'Yb(v): geography and morph',name,cut,len(sites),helper.SPANS[name],helper)
  pd.DataFrame(exports).to_csv(out/f'{name}_{suffix}_frequencies.tsv.gz',sep='\t',index=False)
  results.append(dict(version=suffix,threshold=cut,SNPs=len(sites),rows=len(rows),PDF=str(pdf)))
 save(out/(name+'_complete.json'),dict(status='COMPLETE',locus=name,BA_source=str(folder),BA_cohort=done,
  display_individuals=len(meta),BA_individuals=len(selection),versions=results,
  cap=False,selection='Every site passing declared BA threshold; pagination exhaustive at 80 columns/page',
  all_threshold_counts=countrows,compact_achieves_target=bool(viable),
  note='Main BA selection membership and score rule unchanged; broader Yb rows are display-only when outside that cohort'))

def report(out):
 summary={n:json.loads((out/(n+'_complete.json')).read_text()) for n in ['Dennis','Ybv']}
 save(out/'complete.json',dict(status='COMPLETE',loci=summary,main_figures_changed=False,downloads=0,alignments=0,calling=0))
 with zipfile.ZipFile(out/'Figure3_morph_supplements.zip','w',zipfile.ZIP_DEFLATED) as z:
  for p in sorted(out.iterdir()):
   if p.is_file() and p.suffix not in ['.zip','.py'] and p.name!='jobs.json':z.write(p,p.name)
 for n,d in summary.items():
  for r in d['versions']:print(n,r['version'],'BA='+str(r['threshold']),'SNPs='+str(r['SNPs']),r['PDF'],sep='\t')
 print('ZIP\t'+str(out/'Figure3_morph_supplements.zip'))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=SOURCE)
 p.add_argument('--output',type=Path,default=SOURCE/'final_morph_supplements');p.add_argument('--submit',action='store_true')
 p.add_argument('--overwrite',action='store_true');p.add_argument('--stage',choices=['locus','report']);a=p.parse_args()
 if a.submit:submit(a)
 elif a.stage=='locus':locus(a.output)
 elif a.stage=='report':report(a.output)
 else:p.error('Use --submit')
if __name__=='__main__':main()
