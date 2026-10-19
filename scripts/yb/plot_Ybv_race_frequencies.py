#!/usr/bin/env python3
"""Race-resolved Yb(v) frequencies using completed calls and unchanged frozen markers."""
import argparse,csv,hashlib,importlib.util,json,shlex,shutil,subprocess,sys
from collections import defaultdict
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
DEFAULT=ROOT/'analysis/yb_completed_reads_recovery_v1'
def read(p):
 with Path(p).open() as h:return list(csv.DictReader(h,delimiter='\t'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def combine(vectors):
 import numpy as np
 mat=np.asarray(vectors,dtype=float);n=np.isfinite(mat).sum(axis=0)
 lo=np.min(np.where(np.isfinite(mat),mat,np.inf),axis=0);hi=np.max(np.where(np.isfinite(mat),mat,-np.inf),axis=0)
 return np.where((n>0)&(lo==hi),lo,np.nan)
def draw(matrix,labels,styles,positions,out,pdf,callability=False):
 import numpy as np,matplotlib.pyplot as plt
 from matplotlib.colors import LinearSegmentedColormap
 from matplotlib.patches import ConnectionPatch,Rectangle
 n=len(labels);height=max(5,2.3+n*.185);fig=plt.figure(figsize=(9.4,height))
 left=3.9;right=.3;bottom=1.05;heat_height=n*.185
 ax=fig.add_axes([left/9.4,bottom/height,(9.4-left-right)/9.4,heat_height/height])
 cmap=plt.get_cmap('Greys') if callability else LinearSegmentedColormap.from_list('Rosina',['#D7EDF5','#ECCB63','#662A83'])
 cmap=cmap.copy();cmap.set_bad('#C8C8C8')
 im=ax.pcolormesh(np.arange(len(positions)+1)-.5,np.arange(n+1)-.5,np.ma.masked_invalid(matrix),vmin=0,vmax=1,cmap=cmap,rasterized=False)
 ax.set_ylim(n-.5,-.5);ax.set_xlim(-.5,len(positions)-.5);ax.set_xticks([]);ax.set_yticks(range(n),labels,fontsize=6.4);ax.tick_params(length=0,pad=6)
 for tick,style in zip(ax.get_yticklabels(),styles):
  if style=='focal':tick.set_fontweight('bold')
  elif style=='positive':tick.set_color('#662A83')
  elif style=='negative':tick.set_color('#24748E')
 for j in range(1,n):
  if styles[j]!=styles[j-1]:ax.axhline(j-.5,color='white',lw=1.8)
 locus=fig.add_axes([left/9.4,(bottom+heat_height+.16)/height,(9.4-left-right)/9.4,.43/height])
 lo,hi=1658079,1659072;locus.set_xlim(min(lo-400,min(positions)-30),max(hi+400,max(positions)+30));locus.set_ylim(-.12,.8)
 locus.axhline(0,color='#8a8a8a',lw=.5);locus.add_patch(Rectangle((lo,.08),hi-lo,.26,facecolor='#777',edgecolor='none'))
 locus.set_yticks([]);locus.ticklabel_format(axis='x',style='plain',useOffset=False);locus.tick_params(axis='x',labelsize=6,length=2)
 locus.text((lo+hi)/2,.47,'Yb(v) interval',ha='center',fontsize=7)
 for j,pos in enumerate(positions):fig.add_artist(ConnectionPatch(xyA=(pos,0),coordsA=locus.transData,xyB=(j,-.5),coordsB=ax.transData,color='#9ca2a8',alpha=.24,lw=.3,clip_on=False))
 for aa in (ax,locus):
  for sp in aa.spines.values():sp.set_visible(False)
 fig.text(.025,.97,'Yb(v) | frequencies by race' if not callability else 'Yb(v) | callable individuals by race',fontsize=11,weight='bold',va='top')
 fig.text(.97,.97,'Rosina | NC_140403.1',fontsize=7,ha='right',va='top')
 cb=fig.colorbar(im,cax=fig.add_axes([left/9.4,.61/height,.31,.10/height]),orientation='horizontal')
 cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=6,length=2);cb.outline.set_visible(False)
 cb.set_label('Fraction callable at each SNP' if callability else 'Rosina allele frequency / focal dosage',fontsize=6.5,labelpad=1)
 note='27 unchanged SNP columns; DP >= 2, reported GQ or recorded PL-derived Q >= 20; PK excluded.'
 fig.text(.025,.09/height,note,fontsize=5.8)
 fig.text(.025,.24/height,'n = unique individuals with >=1 callable marker; per-SNP denominators vary. d/v labels describe each wing surface.',fontsize=5.8)
 pdf.savefig(fig);plt.close(fig)
def build(a):
 import numpy as np,pandas as pd
 source=a.source;out=a.output;cfg=json.loads((source/'run.json').read_text());complete=json.loads((source/'complete.json').read_text())
 if complete.get('status')!='COMPLETE':raise RuntimeError('Require complete recovery; do not silently plot failed tasks')
 engine=source/'engine.py'
 if sha(engine)!=cfg['engine_sha']:raise RuntimeError('Frozen extraction engine changed')
 spec=importlib.util.spec_from_file_location('yb_extractor',engine);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 focal=source/'frozen/HmRos_Ybv_unchanged_focal_calls.tsv.gz'
 if sha(focal)!=cfg['focal_hashes'][focal.name]:raise RuntimeError('Focal marker/call hash mismatch')
 sites=pd.read_csv(focal,sep='\t').sort_values('pos1');positions=sites.pos1.astype(int).tolist()
 selected=pd.read_csv(source/'selected_sources.tsv',sep='\t');selected=selected[selected.reference.eq('HmRos')&selected.interval.eq('v')].copy()
 if selected['sample'].duplicated().any():raise RuntimeError('Duplicate sample in final selection')
 members=pd.read_csv(source/'control_membership.tsv',sep='\t');members=members[members.reference.eq('HmRos')&members.interval.eq('v')]
 included=set(members.loc[members.included.isin([True,1,'True','true']),'sample'])
 ph={r['sample']:r for r in read(source/'resolved_phenotypes.tsv')}
 # Reconstruct the same individual unions used in the completed recovery report.
 allpaths={}
 for r in read(Path(cfg['parent'])/'selected_sources.tsv'):
  if r['reference']=='HmRos' and r['interval']=='v':allpaths[r['sample']]=(r['path'],r['source'])
 for t in cfg['new_tasks']:
  p=source/'calls'/f'{t["call_task"]}.vcf.gz'
  if p.is_file():allpaths[t['sample']]=(str(p),'forced_BAM_call')
 merge={}
 mergefile=source/'individual_merges.tsv'
 if mergefile.is_file():
  for r in read(mergefile):
   if r.get('reference')=='HmRos' and r.get('interval')=='v':merge[r['representative']]=r['runs'].split(';')
 cache={}
 def extract(path,sample,kind):
  key=(path,sample,kind)
  if key not in cache:
   x,_=m.extract(path,[sample],sites,legacy=kind=='legacy_VCF',allow_pl=True);cache[key]=x[sample]
  return cache[key]
 values={};identity={};membership=[]
 for r in selected.to_dict('records'):
  s=r['sample'];aliases=merge.get(s,[s]);vectors=[]
  for alias in aliases:
   if alias==s:path,kind=r['path'],r['source']
   else:
    if alias not in allpaths:raise RuntimeError('Missing merged individual source: '+alias)
    path,kind=allpaths[alias]
   vectors.append(extract(path,alias,kind))
  vector=combine(vectors);has=bool(np.isfinite(vector).any())
  if has!=(s in included):raise RuntimeError('Callability disagrees with completed figure: '+s)
  if not has:continue
  if ph.get(s,{}).get('exclusion'):raise RuntimeError('Excluded individual entered plot: '+s)
  ind=ph.get(s,{}).get('individual_id',s)
  if ind in identity:raise RuntimeError('Biological individual counted twice: '+ind)
  identity[ind]=s;values[s]=vector
  race=' '.join(str(r['taxon']).split());d=int(r['d']) if pd.notna(r['d']) else None;v=int(r['v'])
  membership.append(dict(sample=s,individual_id=ind,race=race,d=d,v=v,aliases=';'.join(aliases),source=r['source'],path=r['path'],callable_markers=int(np.isfinite(vector).sum())))
 groups=defaultdict(list)
 for r in membership:groups[(r['race'],r['d'],r['v'])].append(r['sample'])
 def aggregate(ids):
  mat=np.array([values[s] for s in ids]);n=np.isfinite(mat).sum(axis=0)
  mean=np.divide(np.nansum(mat,axis=0),n,out=np.full(len(sites),np.nan),where=n>0)
  return mean,n
 labels=[];freq=[];callable_fraction=[];styles=[];counts=[];rows=[]
 def add(label,ids,style,race='',d=None,v=None):
  mean,n=aggregate(ids);labels.append(label);freq.append(mean);callable_fraction.append(n/len(ids));styles.append(style);counts.append(n)
  for j,pos in enumerate(positions):rows.append(dict(row=label,race=race,d=d,v=v,pos1=pos,n_individuals=len(ids),n_called=int(n[j]),frequency=mean[j]))
 for v in (1,0):
  ss=[r['sample'] for r in membership if r['v']==v];add(f'All races | v{"+" if v else "-"} (n={len(ss)})',ss,'summary',v=v)
  expected=next(p for p in complete['pages'] if p['reference']=='HmRos' and p['interval']=='v')['positive_n' if v else 'negative_n']
  if len(ss)!=expected:raise RuntimeError('Cohort n mismatch')
 # Pooled summaries exactly match the completed figure.
 old=pd.read_csv(source/'HmRos_Ybv_plotted_values.tsv',sep='\t',index_col=0)
 np.testing.assert_allclose(freq[0],old.iloc[0].to_numpy(float),equal_nan=True,rtol=0,atol=1e-12)
 np.testing.assert_allclose(freq[1],old.iloc[3].to_numpy(float),equal_nan=True,rtol=0,atol=1e-12)
 for name,column in [('TNLeft','left_rosina_dosage'),('TNRight','right_rosina_dosage')]:
  vector=sites[column].to_numpy(float);labels.append(name);freq.append(vector);callable_fraction.append(np.where(np.isfinite(vector),1.,np.nan));styles.append('focal');counts.append(np.isfinite(vector).astype(int))
 # A comparison based on the previous cohort, using the current calls, helps
 # distinguish changed population weights from changed genotype processing.
 baseline_ids=set()
 base=a.baseline/'control_membership.tsv'
 if base.is_file():
  bm=pd.read_csv(base,sep='\t');baseline_ids=set(bm.loc[bm.reference.eq('HmRos')&bm.interval.eq('v')&bm.included.isin([True,1,'True','true']),'sample'])
  for v in (1,0):
   ss=[r['sample'] for r in membership if r['v']==v and set(r['aliases'].split(';'))&baseline_ids]
   if ss:add(f'Previous cohort, current calls | v{"+" if v else "-"} (n={len(ss)})',ss,'previous',v=v)
 for r in membership:r['in_previous_cohort']=bool(set(r['aliases'].split(';'))&baseline_ids)
 for (race,d,v),ss in sorted(groups.items(),key=lambda x:(-x[0][2],x[0][0],-1 if x[0][1] is None else x[0][1])):
  ds='?' if d is None else '+' if d else '-'
  label=f'{race.replace("Heliconius ","H. ").replace("ssp. ","")} | d{ds}/v{"+" if v else "-"} (n={len(ss)})'
  add(label,ss,'positive' if v else 'negative',race,d,v)
 out.mkdir(parents=True,exist_ok=False)
 pd.DataFrame(membership).to_csv(out/'race_membership.tsv',sep='\t',index=False)
 pd.DataFrame(rows).to_csv(out/'race_site_frequencies.tsv',sep='\t',index=False,na_rep='NA')
 pd.DataFrame(np.array(freq),index=labels,columns=positions).to_csv(out/'plotted_frequencies.tsv',sep='\t',na_rep='NA')
 pd.DataFrame(np.array(counts),index=labels,columns=positions).to_csv(out/'plotted_callable_counts.tsv',sep='\t')
 import matplotlib;matplotlib.use('Agg');matplotlib.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Sans'})
 from matplotlib.backends.backend_pdf import PdfPages
 with PdfPages(out/'Ybv_frequency_by_race_Rosina.pdf') as pdf:
  draw(np.asarray(freq),labels,styles,positions,out,pdf)
  draw(np.asarray(callable_fraction),labels,styles,positions,out,pdf,True)
 shutil.copy2(__file__,out/'plot_Ybv_race_frequencies.py')
 audit=dict(status='COMPLETE',source=str(source),source_complete_sha256=sha(source/'complete.json'),markers=positions,focal_sha256=sha(focal),focal_calls_unchanged=True,marker_selection_unchanged=True,genotype_filters_unchanged=True,pooled_frequencies_match_completed_figure=True,unique_individuals=len(membership),races=len({r['race'] for r in membership}),race_surface_rows=len(groups),previous_cohort_rows_use_current_calls=True,pages=['race frequencies','per-site callable fractions'],PK_in_controls=False,PDF=str(out/'Ybv_frequency_by_race_Rosina.pdf'))
 (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
 print('INDIVIDUALS\t'+str(len(membership))+'\nRACE_ROWS\t'+str(len(groups))+'\nPDF\t'+audit['PDF'])
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=DEFAULT);p.add_argument('--baseline',type=Path,default=ROOT/'analysis/yb_matched_cohort_v1');p.add_argument('--output',type=Path,required=True);p.add_argument('--submit',action='store_true');a=p.parse_args()
 if a.submit:
  stage=ROOT/'analysis/yb_race_frequency_code';stage.mkdir(exist_ok=True)
  frozen=stage/('plot_'+sha(__file__)[:12]+'.py');shutil.copy2(__file__,frozen)
  args=[sys.executable,str(frozen),'--source',str(a.source),'--baseline',str(a.baseline),'--output',str(a.output)]
  cmd=['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=8G','--time=00:30:00','--job-name=Ybv_race','--output='+str(stage/'plot.%j.out'),'--error='+str(stage/'plot.%j.err'),'--wrap='+shlex.join(args)]
  print('JOB\t'+subprocess.check_output(cmd,text=True).strip());print('EXPECTED_PDF\t'+str(a.output/'Ybv_frequency_by_race_Rosina.pdf'))
 else:build(a)
if __name__=='__main__':main()
