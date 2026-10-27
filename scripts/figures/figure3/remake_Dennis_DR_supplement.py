#!/usr/bin/env python3
"""Remake only the saved BA093 Optix supplement: D/R tracks, concise labels."""
import argparse,hashlib,importlib.util,json,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis')
SOURCE=ROOT/'supp_genoplots_LR_cydno_historical_BDR_20261027'
OUTPUT=ROOT/'Dennis_BA093_LR_DR_clean_20261027'
KEYS=['chrom','pos','ref','alt']
def load(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def inputs(src):
 return [src/'complete.json',src/'worker.py',src/'append_helper.py',src/'run.json',
  src/'tables/Dennis_BA093_selected_SNPs.tsv',src/'tables/Dennis_BA093_frequencies.tsv.gz',
  src/'tables/corrected_LR_calls.tsv',src/'tables/Optix_September_B_D_R_BA.tsv.gz']
def prepare(a):
 if a.output.exists():raise RuntimeError('Output exists; preserving previous figures: '+str(a.output))
 paths=inputs(a.source)
 for p in paths:
  if not p.is_file():raise FileNotFoundError(p)
 if json.loads(paths[0].read_text()).get('status')!='COMPLETE':raise RuntimeError('Source package incomplete')
 old=load(a.source/'worker.py','source_config')
 paths += [old.MORPH/'worker.py',old.MORPH/'figure3_helper.py']
 for p in paths:
  if not p.is_file():raise FileNotFoundError(p)
 a.output.mkdir(parents=True);(a.output/'logs').mkdir()
 shutil.copy2(__file__,a.output/'worker.py')
 save(a.output/'run.json',dict(source=str(a.source),hashes={str(p):digest(p) for p in paths}))
 cmd=[sys.executable,str(a.output/'worker.py'),'--run','--output',str(a.output)]
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:15:00','--mem=8G','--cpus-per-task=1','--job-name=Dennis_DR',
  '--output='+str(a.output/'logs/plot.%j.out'),'--error='+str(a.output/'logs/plot.%j.err'),'--wrap='+shlex.join(cmd)],text=True).strip().split(';')[0]
 save(a.output/'jobs.json',dict(plot=job));print('JOB\t'+job+'\nOUTPUT\t'+str(a.output))

def add_tracks(fig,scores):
 import numpy as np
 from matplotlib.patches import Rectangle
 w,h=fig.get_size_inches();newh=h+2.65;fig.set_size_inches(w,newh)
 for ax in list(fig.axes):
  p=ax.get_position();ax.set_position([p.x0,p.y0*h/newh,p.width,p.height*h/newh])
 for t in list(fig.texts):
  x,y=t.get_position()
  if y>.97:t.set_text('Optix');t.set_position((x,(newh-.22)/newh))
  elif y>.94:t.set_text('176 SNPs · BA ≥ 0.93');t.set_position((x,(newh-.48)/newh))
  else:t.set_position((x,y*h/newh))
 palette={'D':'#E69F00','R':'#CC79A7'}
 for trait,bottom in [('D',h+.80),('R',h-.05)]:
  s=scores.loc[scores.trait.eq(trait)];ax=fig.add_axes([.33,bottom/newh,.65,.62/newh])
  y=np.clip(s.balanced_accuracy.to_numpy(float),.5,1.)
  # One opaque image per dense background cloud, rather than thousands of
  # transparent editable point objects. Highlighted sites remain vector.
  ax.scatter(s.pos,y,s=2.2,c='#E2E6EA',lw=0,rasterized=True,clip_on=False)
  selected=s.historical_highlight.to_numpy(bool)
  ax.scatter(s.loc[selected,'pos'],y[selected],s=10,c=palette[trait],lw=0,clip_on=False)
  ax.axhline(.5,c='#AEB6BF',lw=.5);ax.axhline(.8,c='#666666',ls=(0,(2,2)),lw=.6)
  ax.set(xlim=(1060000,1240000),ylim=(.495,1.015),yticks=[.5,1.],xticks=[])
  ax.tick_params(axis='y',length=0,labelsize=6,pad=3)
  ax.set_ylabel(trait,rotation=0,labelpad=13,fontsize=10,weight='bold')
  for sp in ax.spines.values():sp.set_visible(False)
 ann=fig.add_axes([.33,(h+1.65)/newh,.65,.40/newh]);ann.set(xlim=(1060000,1240000),ylim=(0,1.5),yticks=[],xticks=[])
 ann.add_patch(Rectangle((1088602,.98),9244,.18,facecolor='#3A7CA5',lw=0));ann.text(1093224,1.22,'optix',ha='center',fontsize=8,color='#3A7CA5',style='italic')
 for name,a,b,trait,y in [('R1',1156746,1160220,'R',.4),('R2',1174547,1187833,'R',.4),('D1',1215210,1216124,'D',.4),('D2',1218557,1221494,'D',.1)]:
  ann.add_patch(Rectangle((a,y),b-a,.16,facecolor=palette[trait],lw=0));ann.text((a+b)/2,y+.19,name,ha='center',fontsize=7,color=palette[trait])
 for sp in ann.spines.values():sp.set_visible(False)

def simplify_objects(fig):
 from matplotlib.patches import ConnectionPatch
 from matplotlib.collections import LineCollection
 segments=[]
 for artist in list(fig.artists):
  if isinstance(artist,ConnectionPatch):
   # Resolve both source axes into figure coordinates, preserving endpoints.
   a=fig.transFigure.inverted().transform(artist.coords1.transform(artist.xy1))
   b=fig.transFigure.inverted().transform(artist.coords2.transform(artist.xy2))
   segments.append([a,b]);artist.remove()
 if segments:fig.add_artist(LineCollection(segments,transform=fig.transFigure,colors='#DDDDDD',linewidths=.3,clip_on=False))
 # All displayed data are bounded explicitly; there is no reason to export
 # an axes clipping mask for every cell, connector or annotation artist.
 for artist in fig.findobj():
  if hasattr(artist,'set_clip_on'):artist.set_clip_on(False)
 fig.set_dpi(600)
 return len(segments)

def run(out):
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 from matplotlib.backends.backend_pdf import PdfPages
 cfg=json.loads((out/'run.json').read_text());src=Path(cfg['source'])
 for p,h in cfg['hashes'].items():
  if digest(Path(p))!=h:raise RuntimeError('Input changed: '+p)
 old=load(src/'worker.py','old_style');base=load(src/'append_helper.py','frozen_rows')
 morph=load(old.MORPH/'worker.py','morph_renderer');helper=load(old.MORPH/'figure3_helper.py','style_helper')
 sites=pd.read_csv(src/'tables/Dennis_BA093_selected_SNPs.tsv',sep='\t')
 if len(sites)!=176:raise RuntimeError('Not the supplied 176-SNP figure')
 exp=pd.read_csv(src/'tables/Dennis_BA093_frequencies.tsv.gz',sep='\t',keep_default_na=False,na_values=[''])
 rows,labels=base.frozen_rows(exp,sites,['population','morph','geography','phenotype'],lambda r:helper.short(r.population)+f' | {r.morph} | n={int(r.n_total)}')
 focal=pd.read_csv(src/'tables/corrected_LR_calls.tsv',sep='\t');focal=focal.loc[focal.locus.eq('Dennis')]
 # Two source versions duplicate some sites. Retain only identical duplicates.
 for _,g in focal.groupby(KEYS,dropna=False):
  if len(g[['Left_ALT_dosage','Right_ALT_dosage']].drop_duplicates())!=1:raise RuntimeError('Inconsistent saved focal calls')
 focal=focal.drop_duplicates(KEYS)
 joined=sites[KEYS].merge(focal,on=KEYS,how='left',validate='one_to_one',indicator=True)
 if not joined['_merge'].eq('both').all():raise RuntimeError('Saved focal site missing')
 lr=joined[['Left_ALT_dosage','Right_ALT_dosage']].to_numpy(float).T
 scores=pd.read_csv(src/'tables/Optix_September_B_D_R_BA.tsv.gz',sep='\t')
 if not {'D','R'}<=set(scores.trait):raise RuntimeError('Missing saved D/R BA')
 if scores.historical_highlight.dtype!=bool:
  scores['historical_highlight']=scores.historical_highlight.astype(str).str.lower().map({'true':True,'false':False})
  if scores.historical_highlight.isna().any():raise RuntimeError('Invalid highlight flag')
 path=out/'Dennis_BA093_LR_September_DR_BA.pdf';counts=[];original=PdfPages.savefig
 def finish(self,figure=None,**kwargs):
  add_tracks(figure,scores);old.add_phenotype_table(figure,'Dennis')
  counts.append(simplify_objects(figure))
  return original(self,figure,dpi=600,**kwargs)
 PdfPages.savefig=finish
 try:morph.render(path,sites,np.vstack([rows,lr]),labels+['Left (corrected)','Right (corrected)'],'Optix','Dennis',.93,len(sites),helper.SPANS['Dennis'],helper)
 finally:PdfPages.savefig=original
 (out/'README.txt').write_text('Only the supplied BA093 Optix supplement is remade. Saved frequency rows, selected SNPs, corrected focal calls and September D/R scores/highlights are reused without recalculation. B BA track removed; B phenotype column and original lower genomic annotation retained. Long provenance notes remain here rather than above the figure. Dashed BA line=0.8; coloured points=historical focal-compatible BA>0.99, not independent validation. Gray BA point clouds are 600-dpi raster layers. Genotype cells, connectors, labels, phenotype boxes and highlighted SNPs remain vector. Connectors form one collection; axes clipping masks are disabled. Original files unchanged.\n')
 save(out/'complete.json',dict(status='COMPLETE',PDF=str(path),SNPs=len(sites),BA_tracks=['D','R'],BA_recalculated=False,genotypes_changed=False,connector_segments=counts,raster_layers='gray BA clouds only; 600dpi',source=str(src)))
 print('PDF\t'+str(path)+'\nSTATUS\tCOMPLETE')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=SOURCE);p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--submit',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
 if a.submit:prepare(a)
 elif a.run:run(a.output)
 else:p.error('Use --submit')
if __name__=='__main__':main()
