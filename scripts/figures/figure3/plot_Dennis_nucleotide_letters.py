#!/usr/bin/env python3
"""Display actual nucleotide alleles from the recorded Dennis panel."""
import argparse, hashlib, json, shlex, shutil, subprocess, sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis')
RUN=ROOT/'rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc'

SOURCE=RUN/'Optix_Arnaud_review_nonPK_v2/merged/Dennis_BA90_intervals_frequencies.tsv'

def extract(path):
 import re
 import pandas as pd
 wide=pd.read_csv(path,sep='\t',index_col=0)
 labels=list(wide.index.astype(str))
 if len(labels)!=4 or not labels[0].startswith('Dennis+') or labels[1:3]!=['TNLeft','TNRight'] or not labels[3].startswith('Dennis−'):
  raise RuntimeError('Unexpected source row labels: '+repr(labels))
 for label,n in [(labels[0],142),(labels[3],281)]:
  match=re.search(r'n\s*=\s*(\d+)',label)
  if match is None or int(match.group(1))!=n:raise RuntimeError('Source cohort differs from supplied figure: '+label)
 keys=[]
 for col in wide.columns:
  chrom,pos,ref,alt=col.rsplit(':',3);keys.append((chrom,int(pos),ref,alt))
 if len(keys)!=24 or len(set(keys))!=24:raise RuntimeError('Expected exactly 24 distinct original SNPs')
 # Preserve the source column order and every frequency/dosage unchanged.
 m=wide.apply(pd.to_numeric,errors='raise').copy()
 if ((m<0)|(m>1)).any().any():raise RuntimeError('Dosage outside [0,1]')
 m.index=['present','TN_Left','TN_Right','absent'];m.columns=[k[1] for k in keys]
 if m.columns.duplicated().any():raise RuntimeError('Multiple selected variants at one position')
 records=[]
 for label,row in zip(labels,m.index):
  for key,x in zip(keys,m.loc[row]):
   records.append(dict(chrom=key[0],pos=key[1],ref=key[2],alt=key[3],row=row,original_label=label,frequency=x))
 return pd.DataFrame(records),m

def run(a):
 import numpy as np
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.patches import Ellipse,Rectangle,ConnectionPatch
 from matplotlib.colors import LinearSegmentedColormap
 source=a.values
 if not source.is_file():raise FileNotFoundError('Recorded original plot matrix missing: '+str(source))
 t,m=extract(source)
 import pandas as pd
 sites_path=source.with_name('Dennis_BA90_intervals_sites.tsv')
 sites=pd.read_csv(sites_path,sep='\t')
 keys=['chrom','pos','ref','alt']
 marker=t[keys].drop_duplicates()
 if sites.duplicated(keys).any():raise RuntimeError('Duplicate site orientation')
 marker=marker.merge(sites[keys+['pattern_is_alt']],on=keys,how='left',validate='one_to_one')
 if marker.pattern_is_alt.isna().any():raise RuntimeError('Missing recorded allele orientation')
 def as_bool(x):
  if str(x).lower() in ('true','1'):return True
  if str(x).lower() in ('false','0'):return False
  raise RuntimeError('Unrecognised allele orientation: '+str(x))
 orientations={r.pos:as_bool(r.pattern_is_alt) for r in marker.itertuples()}
 alleles={r.pos:(r.alt,r.ref) if orientations[r.pos] else (r.ref,r.alt) for r in marker.itertuples()}
 letters=[]
 a.output.mkdir(parents=True,exist_ok=False)
 plt.rcParams.update({'pdf.fonttype':42,'font.size':8,'font.family':'DejaVu Sans'})
 fig=plt.figure(figsize=(12,4.2));heat=fig.add_axes([.38,.26,.57,.38]);locus=fig.add_axes([.38,.77,.57,.09]);wing=fig.add_axes([.02,.26,.12,.38]);wing.set_xlim(-1,1);wing.set_ylim(3.5,-.5);wing.axis('off')
 # Small generic forewing cartoons, explicitly phenotype schematics.
 for y,present in [(0,True),(3,False)]:
  wing.add_patch(Ellipse((0,y),1.35,.65,angle=-12,facecolor='#252525',edgecolor='none'))
  if present:wing.add_patch(Ellipse((-.40,y+.07),.38,.26,angle=-12,facecolor='#D44936',edgecolor='none'))
 wing.text(0,-.75,'Forewing base',ha='center',fontsize=7)
 cmap=LinearSegmentedColormap.from_list('Dennis',[(0,'#D7EDF5'),(.5,'#ECCB63'),(1,'#662A83')]);cmap.set_bad('#CCCCCC')
 heat.set_xlim(-.5,23.5);heat.set_ylim(3.5,-.5)
 for j,pos in enumerate(m.columns):
  associated,other=alleles[pos]
  for i,row in enumerate(m.index):
   value=m.loc[row,pos]
   if not np.isfinite(value):
    heat.text(j,i,'–',ha='center',va='center',color='#999999');continue
   if i in (1,2):
    if not any(np.isclose(value,d,atol=1e-8) for d in (0,.5,1)):raise RuntimeError('Focal dosage is not a discrete diploid call')
    bases=(other,other) if np.isclose(value,0) else (associated,associated) if np.isclose(value,1) else (other,associated)
    # Slash denotes an unphased genotype, not a haplotype assignment.
    heat.text(j,i,'/'.join(bases),ha='center',va='center',fontsize=9,fontfamily='monospace',color='#333333')
    display='/'.join(bases)
   else:
    heat.text(j,i-.17,associated,ha='center',va='center',fontsize=11,fontfamily='monospace',weight='bold',color='#662A83')
    heat.text(j,i+.05,f'{value:.0%}',ha='center',va='center',fontsize=6.5,color='#662A83')
    heat.text(j,i+.31,f'{other} {(1-value):.0%}',ha='center',va='center',fontsize=6.5,color='#28708A')
    display=f'{associated} {value:.6f}; {other} {1-value:.6f}'
   letters.append(dict(pos=int(pos),row=row,associated_base=associated,other_base=other,associated_frequency=float(value),display=display))
 heat.set_yticks(range(4),['Dennis present (142)','Left wing','Right wing','Dennis absent (281)']);heat.tick_params(length=0,pad=9);heat.set_xticks([])
 for y in [.5,1.5,2.5]:heat.axhline(y,color='white',lw=2)
 split=sum(m.columns<1218557)
 if 0<split<24:heat.axvline(split-.5,color='white',lw=2)
 for name,lo,hi in [('D1',1215210,1216124),('D2',1218557,1221494)]:
  locus.add_patch(Rectangle((lo,.1),hi-lo,.18,facecolor='#555555',edgecolor='none'));locus.text((lo+hi)/2,.45,name,ha='center',fontsize=8)
 locus.set_xlim(1214500,1222500);locus.set_ylim(0,1);locus.set_yticks([]);locus.set_xticks([1215000,1220000],["1.215","1.220"]);locus.tick_params(length=2,labelsize=7);locus.set_xlabel('Chromosome 18 (Mb)',fontsize=7,labelpad=2)
 for j,pos in enumerate(m.columns):fig.add_artist(ConnectionPatch(xyA=(pos,.1),coordsA=locus.transData,xyB=(j,-.5),coordsB=heat.transData,color='#B3B8BC',alpha=.4,lw=.35,clip_on=False,zorder=0))
 for ax in [heat,locus]:
  for s in ax.spines.values():s.set_visible(False)
 fig.text(.38,.13,'Purple: Dennis-associated nucleotide     Blue: other nucleotide',fontsize=8)
 fig.text(.02,.94,'Dennis',weight='bold',fontsize=13);fig.text(.95,.94,'24 associated SNPs',ha='right',fontsize=9)
 fig.text(.38,.025,'Controls: both nucleotide frequencies (rounded percentages). Left/Right: corrected unphased genotypes.\nBA ≥ 0.90; all 24 original SNPs retained. Homozygotes repeat the base; dash = unavailable. Wing icons schematic.',fontsize=6.5)
 dest=a.output/'Dennis_nucleotide_letters.pdf';fig.savefig(dest);plt.close(fig)
 pd.DataFrame(letters).to_csv(a.output/'nucleotide_display_audit.tsv',sep='\t',index=False)
 t.to_csv(a.output/'unchanged_plot_values.tsv',sep='\t',index=False)
 (a.output/'complete.json').write_text(json.dumps({'status':'COMPLETE','input':str(source),'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'orientation_table':str(sites_path),'orientation_sha256':hashlib.sha256(sites_path.read_bytes()).hexdigest(),'SNPs':24,'cohort_n':[142,281],'values_changed':False,'orientation':'original Dennis-associated allele, NOT necessarily VCF ALT','phase_inferred':False,'PDF':str(dest)},indent=2)+'\n')
 print('SOURCE\t'+str(source));print('PDF\t'+str(dest))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--values',type=Path,default=SOURCE);p.add_argument('--output',type=Path,default=ROOT/'Dennis_nucleotide_letters_v1');p.add_argument('--submit',action='store_true');a=p.parse_args()
 if a.submit:
  stage=ROOT/'Dennis_nucleotide_letters_code';stage.mkdir(parents=True,exist_ok=True);f=stage/('plot_'+hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:12]+'.py');shutil.copy2(__file__,f)
  args=[sys.executable,str(f),'--output',str(a.output)]
  if a.values:args+=['--values',str(a.values)]
  print('JOB\t'+subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=6G','--cpus-per-task=1','--job-name=Dennis_letters','--output='+str(stage/'plot.%j.out'),'--error='+str(stage/'plot.%j.err'),'--wrap='+shlex.join(args)],text=True).strip())
 else:run(a)
if __name__=='__main__':main()
