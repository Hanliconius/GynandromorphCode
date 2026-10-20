#!/usr/bin/env python3
"""Frozen-cohort Dennis diagnostics and unphased publication alternatives."""
import argparse,hashlib,json,re,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis')
RUN=ROOT/'rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc'
SOURCE=RUN/'Optix_Arnaud_review_nonPK_v2/merged'
KEY=['chrom','pos','ref','alt']

def boolean(x):
 if str(x).lower() in ('true','1','1.0'):return True
 if str(x).lower() in ('false','0','0.0'):return False
 raise ValueError('Invalid Boolean '+str(x))
def clean(x):return ' '.join(str(x).split())
def qualified(c):
 import numpy as np
 gt=c.get('GT') or ()
 if len(gt)!=2 or any(v not in (0,1) for v in gt) or (c.get('DP') or 0)<2 or c.get('GQ') is None or c.get('GQ')<20:return np.nan
 return sum(gt)/2

def run(a):
 import numpy as np,pandas as pd,pysam,matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.patches import Rectangle,Polygon,ConnectionPatch
 from matplotlib.colors import LinearSegmentedColormap
 from matplotlib.backends.backend_pdf import PdfPages
 from plot_Dennis_nucleotide_letters import extract
 source=a.source;values=source/'Dennis_BA90_intervals_frequencies.tsv';sites_path=source/'Dennis_BA90_intervals_sites.tsv';meta_path=source/'phenotype_manifest.tsv';membership_path=source/'Dennis_frequency_group_membership.json'
 original,focal=extract(values)
 sites=pd.read_csv(sites_path,sep='\t');sites=original[KEY].drop_duplicates().merge(sites,on=KEY,validate='one_to_one')
 if len(sites)!=24:raise RuntimeError('Marker identity changed')
 sites['pattern_is_alt']=sites.pattern_is_alt.map(boolean)
 meta=pd.read_csv(meta_path,sep='\t').fillna('');meta['taxon']=meta.taxon.map(clean);meta['group']=meta['group'].map(clean)
 if meta['sample'].duplicated().any():raise RuntimeError('Duplicate individuals')
 for col in ['is_PK','Dennis_selection_eligible']:meta[col]=meta[col].map(boolean)
 meta['D']=pd.to_numeric(meta.D,errors='coerce')
 members=json.loads(membership_path.read_text());positive=next(g for g in members if g['phenotype']==1);negative=next(g for g in members if g['phenotype']==0)
 selected=set(positive['samples'])|set(negative['samples'])
 if selected!=set(meta.loc[meta.Dennis_selection_eligible,'sample']):raise RuntimeError('Selection membership mismatch')
 if (positive['n'],negative['n'])!=(142,281):raise RuntimeError('Cohorts changed')
 vcf=RUN/'Optix_legacy_merge_v1/Optix_legacy_plus_recalled.vcf.gz'
 keys=[tuple(r[k] for k in KEY) for r in sites.to_dict('records')];indices={k:j for j,k in enumerate(keys)}
 names=meta['sample'].tolist();matrix=np.full((len(meta),24),np.nan);seen=set()
 with pysam.VariantFile(str(vcf)) as vf:
  if set(names)!=set(vf.header.samples):raise RuntimeError('Metadata and comparative VCF sample sets differ')
  for rec in vf:
   if len(rec.alts or ())!=1:continue
   key=(rec.contig,rec.pos,rec.ref,rec.alts[0])
   if key not in indices:continue
   if key in seen:raise RuntimeError('Duplicate selected VCF variant')
   seen.add(key);j=indices[key];d=np.array([qualified(rec.samples[n]) for n in names]);matrix[:,j]=d if sites.iloc[j].pattern_is_alt else 1-d
 if seen!=set(keys):raise RuntimeError('Selected variants absent from VCF')
 lookup={n:i for i,n in enumerate(names)}
 def pool(ids,min_rate=.8):
  x=matrix[[lookup[n] for n in ids]];n=np.isfinite(x).sum(axis=0)
  return np.divide(np.nansum(x,axis=0),n,out=np.full(24,np.nan),where=(n>0)&(n>=len(ids)*min_rate)),n
 for group,row in [(positive,'present'),(negative,'absent')]:
  v,n=pool(group['samples'])
  if not np.allclose(v,focal.loc[row].to_numpy(),equal_nan=True,atol=1e-8):raise RuntimeError('Freshly read control values differ from frozen figure; stop before plotting')
 # Stable groups use metadata, never genotype concordance.
 meta['race']=meta.apply(lambda r:r['group'] if r.is_PK or not r.taxon else r.taxon,axis=1)
 meta['D_label']=meta.D.map(lambda x:'?' if pd.isna(x) else str(int(x)))
 meta['label_single_line']=meta.apply(lambda r:clean(r['race'])+' | '+r['sample']+' | D'+r.D_label+(' | PK' if r.is_PK else '')+(' | BA cohort' if r['sample'] in selected else ''),axis=1)
 order=meta.sort_values(['is_PK','race','D_label','sample']).index.tolist()
 race_rows=[];race_export=[]
 for (race,state),mm in meta.groupby(['race','D_label'],sort=True):
  ids=mm['sample'].tolist();v,n=pool(ids)
  race_rows.append((f'{race} | D{state} | n={len(ids)}',v,'frequency'))
  for j,key in enumerate(keys):race_export.append(dict(race=race,D=state,n_total=len(ids),n_called=int(n[j]),frequency=v[j],**dict(zip(KEY,key))))
 a.output.mkdir(parents=True,exist_ok=False)
 plt.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':8})
 colors=['#D55E00','#54278F'];cmap=LinearSegmentedColormap.from_list('associated',[colors[0],'#F4ECDD',colors[1]]);cmap.set_bad('#D0D0D0')
 pages=[]
 def draw(rows,title,pdf,style='heat',context=False,ba=False):
  nr=len(rows);h=max(4.2,2.4+nr*.23);label_width=max(4.2,max(len(r[0]) for r in rows)*.045);width=max(13,label_width+7);fig=plt.figure(figsize=(width,h));left=label_width/width
  heat=fig.add_axes([left,.21, .98-left,.53]);loc=fig.add_axes([left,.82,.98-left,.09])
  positions=sites.pos.to_numpy();lo,hi=(1060000,1240000) if context else (1214500,1223500)
  if context and a.optix_gene:lo=min(lo,a.optix_gene[0]-3000);hi=max(hi,a.optix_gene[1]+3000)
  loc.set_xlim(lo,hi);loc.set_ylim(0,1.5);loc.axhline(.12,color='#999999',lw=.5)
  windows=[('B',1111500,1117000),('R1',1156000,1161000),('R2',1174000,1189000),('D1',1215000,1217000),('D2',1218000,1223000)]
  for label,x,y in windows:
   if y<lo or x>hi:continue
   loc.add_patch(Rectangle((x,.2),y-x,.2,color='#555555' if label.startswith('D') else '#CCCCCC',lw=0));loc.text((x+y)/2,.55 if label!='D2' else .95,label,ha='center',fontsize=7)
  if a.optix_gene and context:
   x,y=a.optix_gene;loc.add_patch(Rectangle((x,1.05),y-x,.15,color='#333333'));loc.text((x+y)/2,1.25,'optix',ha='center',fontsize=7)
  loc.set_yticks([]);ticks=np.linspace(lo,hi,4);loc.set_xticks(ticks,[f'{x/1e6:.3f}' for x in ticks],fontsize=7);loc.set_xlabel('Chromosome 18 (Mb)',fontsize=7,labelpad=1)
  for j,pos in enumerate(positions):fig.add_artist(ConnectionPatch(xyA=(pos,.12),coordsA=loc.transData,xyB=(j,-.5),coordsB=heat.transData,color='#BFC3C7',lw=.3,alpha=.4,clip_on=False))
  heat.set_xlim(-.5,23.5);heat.set_ylim(nr-.5,-.5);heat.set_xticks([]);heat.set_yticks(range(nr),[r[0] for r in rows],fontsize=7);heat.tick_params(length=0,pad=8)
  if style=='heat':heat.imshow(np.ma.masked_invalid(np.array([r[1] for r in rows])),vmin=0,vmax=1,cmap=cmap,aspect='auto',interpolation='nearest')
  else:
   for i,(label,vals,kind) in enumerate(rows):
    for j,x in enumerate(vals):
     r=sites.iloc[j];associated,other=(r.alt,r.ref) if r.pattern_is_alt else (r.ref,r.alt)
     if not np.isfinite(x):heat.text(j,i,'–',ha='center',va='center',color='#999999');continue
     if kind=='genotype':
      if x not in (0,.5,1):raise RuntimeError('Non-discrete focal genotype')
      if x==.5:
       heat.add_patch(Polygon([(j-.47,i-.43),(j+.47,i-.43),(j-.47,i+.43)],facecolor='#E4DDED',edgecolor='none'))
       heat.add_patch(Polygon([(j+.47,i-.43),(j+.47,i+.43),(j-.47,i+.43)],facecolor='#FAE4D4',edgecolor='none'))
       heat.text(j-.17,i-.15,associated,color=colors[1],ha='center',va='center',fontsize=9,weight='bold');heat.text(j+.17,i+.17,other,color=colors[0],ha='center',va='center',fontsize=9,weight='bold')
      else:
       c=1 if x==1 else 0;heat.add_patch(Rectangle((j-.47,i-.43),.94,.86,facecolor='#E4DDED' if c else '#FAE4D4',edgecolor='none'));heat.text(j,i,associated if c else other,ha='center',va='center',color=colors[c],weight='bold',fontsize=10)
     elif kind=='consensus':
      base=associated if x>.5 else other if x<.5 else associated+'/'+other;heat.text(j,i,base,ha='center',va='center',color=colors[int(x>.5)],fontsize=10,weight='bold')
     else:
      heat.add_patch(Rectangle((j-.47,i-.43),.94,.86,facecolor=cmap(x),edgecolor='white',lw=.3));heat.text(j,i,associated if x>.5 else other,ha='center',va='center',fontsize=8,color='white' if x>.8 else '#333333')
  split=sum(positions<1218000);heat.axvline(split-.5,color='#AAAAAA',lw=.6)
  for ax in (heat,loc):
   for spine in ax.spines.values():spine.set_visible(False)
  fig.text(.025,.96,title,weight='bold',fontsize=12)
  if ba:
   for j,x in enumerate(sites.balanced_accuracy):heat.text(j,1.025,f'{x:.2f}',transform=heat.get_xaxis_transform(),ha='center',fontsize=6)
   heat.text(-.025,1.025,'BA',transform=heat.transAxes,ha='right',fontsize=6)
  note='24 frozen BA ≥ 0.90 SNPs; PK excluded from selection. Missing calls are grey/dashes; DP ≥ 2 and reported GQ ≥ 20.\n'
  note+='Purple = Dennis-present-associated allele; orange = Dennis-absent-associated allele. Frequencies average qualified individual dosages.\n'
  note+='Consensus = per-site majority base across selection controls, not a haplotype. Split cells = unphased heterozygotes; single base = homozygote. Corrected L/R.'
  fig.text(.025,.055,note,fontsize=7)
  if style=='heat':
   cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap=cmap),cax=fig.add_axes([.73,.11,.24,.015]),orientation='horizontal',ticks=[0,.5,1]);cb.ax.tick_params(labelsize=7);cb.set_label('Associated allele frequency / dosage',fontsize=7)
  pdf.savefig(fig);plt.close(fig);pages.append(dict(title=title,rows=nr,context=context,style=style))
 def diagnostic(filename,items,title,ba=False):
  with PdfPages(a.output/filename) as pdf:
   for start in range(0,len(items),42):
    draw(items[start:start+42],title+f' | rows {start+1}–{min(start+42,len(items))}',pdf,ba=ba)
 focal_rows=[('Left (corrected)',focal.loc['TN_Left'].to_numpy(),'genotype'),('Right (corrected)',focal.loc['TN_Right'].to_numpy(),'genotype')]
 allrows=[(meta.loc[i,'label_single_line'],matrix[i],'genotype') for i in order]
 diagnostic('01_Dennis_all_individuals.pdf',focal_rows+allrows,'Dennis | every included comparative individual')
 diagnostic('02_Dennis_race_frequencies.pdf',focal_rows+race_rows,'Dennis | race and phenotype frequencies')
 selrows=[(meta.loc[i,'label_single_line'],matrix[i],'genotype') for i in order if meta.loc[i,'sample'] in selected]
 diagnostic('03_Dennis_BA_selection_individuals.pdf',focal_rows+selrows,'Dennis | individuals used for BA selection',True)
 # Choose four morphs per phenotype by sample count, without inspecting SNP agreement.
 chosen={};morphrows={}
 for state in [1,0]:
  mm=meta[meta['sample'].isin(selected)&meta.D.eq(state)];counts=mm.groupby('race').size().sort_values(ascending=False,kind='stable');races=counts.index[:4].tolist();chosen[state]=races;morphrows[state]=[]
  for race in races:
   ids=mm.loc[mm.race.eq(race),'sample'].tolist();v,n=pool(ids);morphrows[state].append((f'{race} | D{state} | n={len(ids)}',v,'frequency'))
 bread=[('Dennis present consensus (n=142)',focal.loc['present'].to_numpy(),'consensus'),('Dennis absent consensus (n=281)',focal.loc['absent'].to_numpy(),'consensus')]
 with PdfPages(a.output/'04_Dennis_publication_alternatives.pdf') as pdf:
  for context in [False,True]:
   draw([bread[0]]+focal_rows+[bread[1]],'Dennis | consensus / corrected genotypes',pdf,'letters',context)
   draw(morphrows[1]+[bread[0]]+focal_rows+[bread[1]]+morphrows[0],'Dennis | morph frequencies explain the consensus',pdf,'letters',context)
 meta['included_in_BA']=meta['sample'].isin(selected);meta.to_csv(a.output/'individual_membership.tsv',sep='\t',index=False)
 pd.DataFrame(matrix,index=names,columns=[':'.join(map(str,k)) for k in keys]).to_csv(a.output/'all_individual_associated_dosages.tsv.gz',sep='\t')
 pd.DataFrame(race_export).to_csv(a.output/'race_frequencies_and_callability.tsv.gz',sep='\t',index=False)
 sites.to_csv(a.output/'frozen_markers_and_BA.tsv',sep='\t',index=False)
 pd.DataFrame([dict(phenotype=state,race=race,selection_rule='rank by eligible sample count, no genotype matching') for state,races in chosen.items() for race in races]).to_csv(a.output/'publication_morph_membership.tsv',sep='\t',index=False)
 provenance={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [values,sites_path,meta_path,membership_path]}
 (a.output/'complete.json').write_text(json.dumps(dict(status='COMPLETE',sources=provenance,vcf=str(vcf),comparative_individuals=len(meta),BA_individuals=len(selected),SNPs=24,corrected_focal_matches_frozen_matrix=True,control_frequencies_verified=True,phase_inferred=False,chosen_morphs=chosen,morph_selection='four largest eligible race groups per D phenotype, no genotype matching',optix_gene_coordinates=a.optix_gene,optix_gene_annotation_status='explicit coordinates supplied' if a.optix_gene else 'not supplied; wide locus view contains B R D windows only',pages=pages,PDFs=[str(p) for p in sorted(a.output.glob('*.pdf'))]),indent=2)+'\n')
 for p in sorted(a.output.glob('*.pdf')):print('PDF\t'+str(p))
 print(f'COMPLETE\t{len(meta)} individuals; {len(selected)} BA controls; 24 unchanged SNPs')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=SOURCE);p.add_argument('--output',type=Path,default=ROOT/'Dennis_figure_set_v1');p.add_argument('--optix-gene',type=int,nargs=2,metavar=('START','END'));p.add_argument('--submit',action='store_true');a=p.parse_args()
 if a.submit:
  helper=Path(__file__).with_name('plot_Dennis_nucleotide_letters.py');digest=hashlib.sha256(Path(__file__).read_bytes()+helper.read_bytes()).hexdigest()[:12];stage=ROOT/'Dennis_figure_set_code'/digest;stage.mkdir(parents=True,exist_ok=True)
  for f in [Path(__file__),helper]:shutil.copy2(f,stage/f.name)
  args=[sys.executable,str(stage/Path(__file__).name),'--source',str(a.source),'--output',str(a.output)]
  if a.optix_gene:args+=['--optix-gene',*map(str,a.optix_gene)]
  job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=8G','--cpus-per-task=1','--job-name=Dennis_set','--output='+str(stage/'plot.%j.out'),'--error='+str(stage/'plot.%j.err'),'--wrap='+shlex.join(args)],text=True).strip();print('JOB\t'+job);print('LOG_DIR\t'+str(stage));print('OUTPUT\t'+str(a.output))
 else:run(a)
if __name__=='__main__':main()
