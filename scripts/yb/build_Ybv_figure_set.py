#!/usr/bin/env python3
"""Yb(v) figure set from frozen v7 membership and original 27 marker definitions."""
import argparse,hashlib,json,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
BASES={'A':'#26943A','C':'#2463C5','G':'#D89500','T':'#D52C35'}
def clean(s):return ' '.join(str(s).split())
def run(a):
 import numpy as np,pandas as pd
 frozen=a.source
 audit=json.loads((frozen/'audit.json').read_text())
 if audit['status']!='COMPLETE' or audit['allele_scale']!='0=REF,1=ALT':raise RuntimeError('Require complete v7 ALT panel')
 members=pd.read_csv(frozen/'membership_geography.tsv',sep='\t').fillna('')
 if members['sample'].duplicated().any():raise RuntimeError('Duplicate representative individual')
 members['v']=pd.to_numeric(members.v,errors='coerce')
 source=Path(audit['source'])
 old=pd.read_csv(ROOT/'analysis/yb_completed_reads_recovery_v1/frozen/HmRos_Ybv_unchanged_focal_calls.tsv.gz',sep='\t').sort_values('pos1')
 if len(old)!=27:raise RuntimeError('Original short marker count differs')
 if a.with_gene_annotations:
  cached=ROOT/'analysis/Ybv_figure_set_v1/individual_ALT_dosages.tsv.gz'
  completed=json.loads((cached.parent/'complete.json').read_text())
  if completed['status']!='COMPLETE':raise RuntimeError('Original figure set incomplete')
  stored=pd.read_csv(cached,sep='\t',index_col=0)
  positions=old.pos1.astype(int).tolist()
  if set(stored.index)!=set(members['sample']) or list(map(int,stored.columns))!=positions:raise RuntimeError('Cached matrix identity differs from frozen v7')
  matrix=stored.loc[members['sample']].to_numpy(dtype=float);marker_audit=[]
  gene_features=json.loads(Path(__file__).with_name('Ybv_gene_features.json').read_text())
 else:
  import pysam,importlib.util
  spec=importlib.util.spec_from_file_location('yb_quality',ROOT/'analysis/yb_completed_reads_recovery_v1/engine.py');engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
  positions=old.pos1.astype(int).tolist();matrix=np.full((len(members),len(old)),np.nan)
  wanted={int(r.pos1):(j,r) for j,r in enumerate(old.itertuples())};seen=set();marker_audit=[];raw_audit=[];record_vectors={};incompatible=np.zeros(matrix.shape,dtype=bool)
  with pysam.VariantFile(str(source/'full_window.vcf.gz')) as vf:
   for rec in vf.fetch('NC_140403.1',min(positions)-1,max(positions)):
    if rec.pos not in wanted:continue
    j,r=wanted[rec.pos]
    if len(rec.ref)!=1 or any(len(x)!=1 or x not in 'ACGT' for x in (rec.alts or ())):
     marker_audit.append(dict(pos1=rec.pos,REF=rec.ref,old_ALT=r.alt,joint_ALT=','.join(rec.alts or ()),status='separate_indel_or_symbolic_record_not_used_for_SNP',n_called=0));continue
    recvector=np.full(len(members),np.nan)
    seen.add(rec.pos)
    if rec.ref!=r.ref:raise RuntimeError('Reference mismatch at original marker '+str(rec.pos))
    alleles=rec.alleles
    status='exact_old_REF_ALT' if tuple(rec.alts or ())==(r.alt,) else 'reference_only' if not rec.alts else 'different_or_multiple_ALT'
    for i,member in enumerate(members.to_dict('records')):
     vals=[]
     for alias in member['aliases'].split(';'):
      c=rec.samples[alias];value,quality=engine.quality(c,allow_pl=True)
      gt=c.get('GT') or ();bases=[alleles[k] for k in gt if k is not None and k<len(alleles)]
      compatible=len(gt)==2 and len(bases)==2 and all(x in (r.ref,r.alt) for x in bases)
      if value is not None and not compatible:incompatible[i,j]=True
      alt_value=sum(x==r.alt for x in bases)/2 if compatible and value is not None else np.nan
      if np.isfinite(alt_value):vals.append(alt_value)
      raw_audit.append(dict(sample=alias,representative=member['sample'],pos1=rec.pos,old_REF=r.ref,old_ALT=r.alt,joint_alleles=','.join(alleles),GT=str(gt),DP=c.get('DP'),GQ=c.get('GQ'),quality=quality,target_ALT_dosage=alt_value))
     if vals and len(set(vals))==1:recvector[i]=vals[0]
     elif vals:incompatible[i,j]=True
    marker_audit.append(dict(pos1=rec.pos,REF=r.ref,old_ALT=r.alt,joint_ALT=','.join(rec.alts or ()),status=status,n_called=int(np.isfinite(recvector).sum())))
    record_vectors.setdefault(j,[]).append(recvector)
  for j,vectors in record_vectors.items():
   x=np.asarray(vectors);n=np.isfinite(x).sum(axis=0)
   low=np.min(np.where(np.isfinite(x),x,np.inf),axis=0);high=np.max(np.where(np.isfinite(x),x,-np.inf),axis=0)
   matrix[:,j]=np.where((n>0)&(low==high)&~incompatible[:,j],low,np.nan)
   if len(vectors)>1:marker_audit.append(dict(pos1=positions[j],status='duplicate_SNP_records_concordant_union_conflicts_missing',records=len(vectors),n_conflicts=int(((n>0)&((low!=high)|incompatible[:,j])).sum()),n_called=int(np.isfinite(matrix[:,j]).sum())))
  for pos in set(positions)-seen:marker_audit.append(dict(pos1=pos,status='absent_from_complete_VCF_kept_missing',n_called=0))
 def pool(ids):
  x=matrix[ids];n=np.isfinite(x).sum(axis=0)
  return np.divide(np.nansum(x,axis=0),n,out=np.full(27,np.nan),where=n>0),n
 frozen_freq=pd.read_csv(frozen/'region_race_ALT_frequencies.tsv',sep='\t')
 summaries={}
 for state in [1,0]:
  ids=members.index[members.v.eq(state)].tolist()
  f,n=pool(ids);summaries[state]=(ids,f,n)
  label='All retained populations | v'+('+' if state else '-')
  before=frozen_freq[frozen_freq.row.eq(label)].set_index('pos1').reindex(positions)
  if len(before)!=27 or not np.allclose(before.ALT_frequency.to_numpy(),f,equal_nan=True,atol=1e-8):raise RuntimeError('Control frequencies differ from completed v7 plot: '+label)
 focal=[]
 for label,col in [('Left (corrected)','left_rosina_dosage'),('Right (corrected)','right_rosina_dosage')]:
  x=old[col].to_numpy(dtype=float,copy=True);flip=old.rosina_base.eq(old.ref).to_numpy();x[flip]=1-x[flip]
  if any(np.isfinite(v) and not any(np.isclose(v,z) for z in [0,.5,1]) for v in x):raise RuntimeError('Non-discrete corrected focal genotype')
  focal.append((label,x,'genotype'))
 members['single_line_label']=members.apply(lambda r:clean(r['geography'])+' | '+clean(r['race'])+' | '+r['sample']+' | v'+('?' if pd.isna(r.v) else str(int(r.v))),axis=1)
 order=members.sort_values(['geography','race','v','sample']).index.tolist()
 race_rows=[];export=[]
 for (geo,race,state),g in members.groupby(['geography','race','v'],sort=True,dropna=False):
  f,n=pool(g.index.tolist());race_rows.append((clean(geo)+' | '+clean(race)+f' | v{int(state)} | n={len(g)}',f,'frequency'))
  for j,pos in enumerate(positions):export.append(dict(geography=geo,race=race,v=state,pos1=pos,n_total=len(g),n_called=int(n[j]),ALT_frequency=f[j]))
 morphs={};chosen=[]
 for state in [1,0]:
  mm=members[members.v.eq(state)]
  ranks=mm.groupby(['geography','race'],sort=True).size().sort_values(ascending=False,kind='stable')
  morphs[state]=[]
  for geo,race in list(ranks.index)[:4]:
   ids=mm.index[mm.geography.eq(geo)&mm.race.eq(race)].tolist();f,n=pool(ids)
   morphs[state].append((clean(race)+f' | v{state} | n={len(ids)}',f,'frequency'))
   chosen.extend(dict(sample=members.loc[i,'sample'],race=race,geography=geo,v=state) for i in ids)
 a.output.mkdir(parents=True,exist_ok=False)
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.patches import Rectangle,Polygon,ConnectionPatch
 from matplotlib.colors import LinearSegmentedColormap,to_rgba
 from matplotlib.backends.backend_pdf import PdfPages
 plt.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':8})
 cmap=LinearSegmentedColormap.from_list('ALT',['#D55E00','#F4ECDD','#54278F']);cmap.set_bad('#D0D0D0')
 pages=[]
 def draw(rows,title,pdf,letters=False,nt=False,wide=False):
  nr=len(rows);lw=max(4.3,max(len(r[0]) for r in rows)*.044);width=max(13,lw+7);left=lw/width
  heights=np.array([2.7 if a.bold_focal and letters and row[2]=='genotype' else 1. for row in rows]);centres=np.cumsum(heights)-heights/2-.5;total_height=float(heights.sum())
  fig=plt.figure(figsize=(width,max(4.5,2.5+total_height*.25)));ax=fig.add_axes([left,.23,.98-left,.49]);locus=fig.add_axes([left,.83,.98-left,.09])
  lo,hi=((1630000,1692000) if a.with_gene_annotations else (1630000,1661500)) if wide else (min(positions)-180,max(positions)+180)
  locus.set_xlim(lo,hi);locus.set_ylim(0,2.6 if a.with_gene_annotations and wide else 1.5);locus.axhline(.1,color='#AAAAAA',lw=.5)
  for label,start,end in [('Yb(d)',1633041,1635646),('Yb(v)',1658079,1659072)]:
   if end<lo or start>hi:continue
   locus.add_patch(Rectangle((start,.2),end-start,.2,color='#555555'));locus.text((start+end)/2,.65,label,ha='center',fontsize=8)
  if a.with_gene_annotations and wide:
   for name,y,color in [('cortex',1.25,'#454545'),('ivory',2.0,'#267B82')]:
    feature=gene_features[name];exons=feature['exons'];start=min(x[0] for x in exons);end=max(x[1] for x in exons)
    locus.plot([max(lo,start),min(hi,end)],[y,y],color=color,lw=.7)
    for ex_start,ex_end in exons:
     if ex_end<lo or ex_start>hi:continue
     locus.add_patch(Rectangle((max(lo,ex_start),y-.09),min(hi,ex_end)-max(lo,ex_start)+1,.18,facecolor=color,lw=0))
    if name=='ivory':
     locus.annotate('',xy=(hi-500,y),xytext=(min(end,hi-3500),y),arrowprops=dict(arrowstyle='->',color=color,lw=.8))
     locus.text(start+3000,y+.22,'ivory → (continues)',fontsize=7,color=color,ha='left')
    else:
     locus.annotate('',xy=(start,y),xytext=(end,y),arrowprops=dict(arrowstyle='->',color=color,lw=.7))
     locus.text((start+end)/2,y+.22,'cortex',fontsize=7,color=color,ha='center')
  ticks=np.linspace(lo,hi,4);locus.set_xticks(ticks,[f'{p/1e6:.4f}' for p in ticks],fontsize=7);locus.set_yticks([]);locus.set_xlabel('Chromosome 15 (Rosina, Mb)',fontsize=7,labelpad=1)
  for j,pos in enumerate(positions):fig.add_artist(ConnectionPatch(xyA=(pos,.1),coordsA=locus.transData,xyB=(j,-.5),coordsB=ax.transData,color='#BFC3C7',lw=.3,alpha=.4,clip_on=False))
  ax.set_xlim(-.5,26.5);ax.set_ylim(total_height-.5,-.5);ax.set_xticks([]);ax.set_yticks(centres,[r[0] for r in rows],fontsize=7);ax.tick_params(length=0,pad=8)
  if not letters:ax.imshow(np.ma.masked_invalid(np.asarray([r[1] for r in rows])),cmap=cmap,vmin=0,vmax=1,aspect='auto',interpolation='nearest')
  else:
   for i,(label,vals,kind) in enumerate(rows):
    cy=centres[i];rh=heights[i];bold=a.bold_focal and kind=='genotype'
    for j,x in enumerate(vals):
     ref,alt=old.iloc[j]['ref'],old.iloc[j]['alt']
     colour=lambda base:BASES[base] if nt else '#54278F' if base==alt else '#D55E00'
     if not np.isfinite(x):ax.text(j,cy,'–',ha='center',va='center',color='#AAAAAA');continue
     if kind=='genotype' and np.isclose(x,.5):
      ax.add_patch(Polygon([(j-.47,cy-.43*rh),(j+.47,cy-.43*rh),(j-.47,cy+.43*rh)],facecolor=to_rgba(colour(ref),.14),edgecolor='none'))
      ax.add_patch(Polygon([(j+.47,cy-.43*rh),(j+.47,cy+.43*rh),(j-.47,cy+.43*rh)],facecolor=to_rgba(colour(alt),.14),edgecolor='none'))
      ax.text(j-.17,cy-.16*rh,ref,color=colour(ref),ha='center',va='center',fontsize=14 if bold else 9,weight='bold');ax.text(j+.17,cy+.17*rh,alt,color=colour(alt),ha='center',va='center',fontsize=14 if bold else 9,weight='bold')
     else:
      base=alt if x>.5 else ref if x<.5 else ref+'/'+alt
      if kind=='genotype':
       ax.add_patch(Rectangle((j-.47,cy-.43*rh),.94,.86*rh,facecolor=to_rgba(colour(base),.14),edgecolor='none'))
      elif kind=='frequency':ax.add_patch(Rectangle((j-.47,cy-.43*rh),.94,.86*rh,facecolor=to_rgba(cmap(x),.18),edgecolor='white',lw=.4))
      ax.text(j,cy,base,ha='center',va='center',color=colour(base) if base in BASES else '#555555',fontsize=15 if bold else 10 if kind!='frequency' else 8,weight='bold')
  if a.bold_focal and letters:
   for i,row in enumerate(rows):
    if row[2]!='genotype':continue
    for j in range(27):ax.add_patch(Rectangle((j-.47,centres[i]-.43*heights[i]),.94,.86*heights[i],facecolor='none',edgecolor='#777777',lw=.65))
    tick=ax.get_yticklabels()[i];tick.set_fontweight('bold');tick.set_fontsize(10)
  for z in [ax,locus]:
   for spine in z.spines.values():spine.set_visible(False)
  fig.text(.025,.96,title,weight='bold',fontsize=12)
  if not letters:
   cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap=cmap),cax=fig.add_axes([.73,.12,.24,.015]),orientation='horizontal',ticks=[0,.5,1]);cb.ax.tick_params(labelsize=7);cb.set_label('ALT frequency / corrected focal dosage',fontsize=7)
  legend='A green; C blue; G orange; T red.' if nt else 'Purple = ALT; orange = REF. Colours do not encode phenotype.'
  fig.text(.025,.045,legend+'\nOriginal 27 SNPs; v7 individuals/exclusions unchanged. Consensus = per-site majority base, not a haplotype.\nCorrected Left/Right: one genotype per row; split cell = unphased heterozygote, one base = homozygote. Missing = grey/dash.\nMorph rows: majority letter and pale ALT-frequency background; four largest geography/race groups per phenotype, not chosen for agreement.',fontsize=7)
  pdf.savefig(fig);plt.close(fig);pages.append(dict(title=title,rows=nr,ATGC_colours=nt,wide=wide))
 def diagnostic(filename,rows,title):
  with PdfPages(a.output/filename) as pdf:
   for start in range(0,len(rows),42):draw(rows[start:start+42],title+f' | rows {start+1}–{min(start+42,len(rows))}',pdf)
 if not a.with_gene_annotations:diagnostic('01_Ybv_all_retained_individuals.pdf',focal+[(members.loc[i,'single_line_label'],matrix[i],'genotype') for i in order],'Yb(v) | all retained individuals')
 if not a.with_gene_annotations:diagnostic('02_Ybv_race_frequencies.pdf',focal+race_rows,'Yb(v) | geography and race frequencies')
 if not a.with_gene_annotations:diagnostic('03_Ybv_phenotype_comparison_individuals.pdf',focal+[(members.loc[i,'single_line_label'],matrix[i],'genotype') for i in order if members.loc[i,'v'] in (0,1)],'Yb(v) | phenotype-comparison individuals; no new BA selection')
 bread=[(f'Yb(v)+ consensus (n={len(summaries[1][0])})',summaries[1][1],'consensus'),(f'Yb(v)− consensus (n={len(summaries[0][0])})',summaries[0][1],'consensus')]
 if a.with_gene_annotations:
  with PdfPages(a.output/('Ybv_publication_allele_colours_cortex_ivory.pdf' if a.allele_colours else 'Ybv_publication_ATGC_cortex_ivory.pdf')) as pdf:
   draw(morphs[1]+[bread[0]]+focal+[bread[1]]+morphs[0],'Yb(v) | cortex / ivory context',pdf,True,not a.allele_colours,True)
 else:
  with PdfPages(a.output/'04_Ybv_publication_alternatives.pdf') as pdf:
   for wide in [False,True]:
    for nt in [False,True]:
     draw([bread[0]]+focal+[bread[1]],'Yb(v) | consensus and corrected genotypes',pdf,True,nt,wide)
     draw(morphs[1]+[bread[0]]+focal+[bread[1]]+morphs[0],'Yb(v) | morph frequencies and corrected genotypes',pdf,True,nt,wide)
  with PdfPages(a.output/'05_Ybv_publication_last_panel_ATGC.pdf') as pdf:draw(morphs[1]+[bread[0]]+focal+[bread[1]]+morphs[0],'Yb(v) | morph frequencies and corrected genotypes',pdf,True,True,True)
 members.to_csv(a.output/'individual_membership.tsv',sep='\t',index=False)
 pd.DataFrame(matrix,index=members['sample'],columns=positions).to_csv(a.output/'individual_ALT_dosages.tsv.gz',sep='\t')
 pd.DataFrame(export).to_csv(a.output/'race_frequencies_and_callability.tsv.gz',sep='\t',index=False)
 pd.DataFrame(chosen).to_csv(a.output/'publication_morph_membership.tsv',sep='\t',index=False)
 pd.DataFrame(marker_audit).to_csv(a.output/'original_marker_mapping_audit.tsv',sep='\t',index=False)
 provenance={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [frozen/'audit.json',frozen/'membership_geography.tsv',frozen/'region_race_ALT_frequencies.tsv',ROOT/'analysis/yb_completed_reads_recovery_v1/frozen/HmRos_Ybv_unchanged_focal_calls.tsv.gz']}
 (a.output/'complete.json').write_text(json.dumps(dict(status='COMPLETE',source=str(frozen),source_hashes=provenance,individuals=len(members),SNPs=27,scale='0=REF,1=ALT',no_new_calls=True,no_phenotype_or_marker_selection_changes=True,control_frequencies_match_v7=True,focal='original corrected calls converted to ALT',phase_inferred=False,bold_focal=a.bold_focal,allele_colours=a.allele_colours,gene_annotation=gene_features if a.with_gene_annotations else None,BA_selection='Not recomputed; diagnostic is phenotype comparator cohort',pages=pages,PDFs=[str(p) for p in sorted(a.output.glob('*.pdf'))]),indent=2)+'\n')
 for p in sorted(a.output.glob('*.pdf')):print('PDF\t'+str(p))
 print('COMPLETE\t'+str(len(members))+' individuals; 27 unchanged markers')
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=ROOT/'analysis/ybv_short_ALT_geography_v7');p.add_argument('--output',type=Path,default=ROOT/'analysis/Ybv_figure_set_v1');p.add_argument('--allele-colours',action='store_true');p.add_argument('--bold-focal',action='store_true');p.add_argument('--with-gene-annotations',action='store_true');p.add_argument('--submit',action='store_true');a=p.parse_args()
 if a.bold_focal:a.with_gene_annotations=True
 if not a.submit:return run(a)
 tag=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:12];stage=ROOT/'analysis/Ybv_figure_set_code'/tag;stage.mkdir(parents=True,exist_ok=True);script=stage/Path(__file__).name;shutil.copy2(__file__,script)
 cmd=[sys.executable,str(script),'--source',str(a.source),'--output',str(a.output)]
 if a.bold_focal:cmd+=['--bold-focal']
 if a.allele_colours:cmd+=['--allele-colours']
 if a.with_gene_annotations:
  annotation=Path(__file__).with_name('Ybv_gene_features.json');shutil.copy2(annotation,stage/annotation.name);cmd+=['--with-gene-annotations']
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=8G','--cpus-per-task=1','--job-name=Ybv_set','--output='+str(stage/'plot.%j.out'),'--error='+str(stage/'plot.%j.err'),'--wrap='+shlex.join(cmd)],text=True).strip()
 print('JOB\t'+job);print('LOG_DIR\t'+str(stage));print('OUTPUT\t'+str(a.output))
if __name__=='__main__':main()
