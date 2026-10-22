#!/usr/bin/env python3
"""Geographic sensitivity of Yb(v) BA using the completed 62-kb call matrix.

No downloads, alignment, calling, mixture fitting, or publication SNP selection.
Regions and phenotypes come from the frozen cohort used by the completed scan.
"""
import argparse,getpass,hashlib,importlib.util,json,re,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
PARENT=ROOT/'analysis/Optix_Ybv_integrated_figure'
COHORT=ROOT/'analysis/ybv_short_ALT_geography_v7'
YF=ROOT/'analysis/Ybv_figure_set_v1'
HO=ROOT/'analysis/LR_window_heterozygosity_v1'
OLD=ROOT/'analysis/yb_completed_reads_recovery_v1/frozen/HmRos_Ybv_unchanged_focal_calls.tsv.gz'
CHROM='NC_140403.1';LO=1630000;HI=1692000
BLOCKED={'ERR1303880','ERR1303882','ERR1303884','ERR1303885','ERR260303','ERR260304'}
COMPARISONS=['pooled','West_only','East_only','equal_regions','West_positive_all_negative']
LABELS={'pooled':'Pooled: all retained v+ versus all retained v-',
 'West_only':'Within West: western v+ versus western v-',
 'East_only':'Within East: eastern v+ versus eastern v-',
 'equal_regions':'Equal regions: mean of West and East within-region BA',
 'West_positive_all_negative':'Western v+ versus all retained v-'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def load(path,name):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def ensure_allowed(members):
 text=members.astype(str).agg(' | '.join,axis=1)
 bad=text.str.contains('numata|pardalinus|piano|bwpk|pk_dark|pk_pale|gynandromorph|tn2l|tn2r|tn3l|tn3r|tn_left|tn_right',case=False,regex=True)
 ids=members['sample'].astype(str).isin(BLOCKED)
 if 'aliases' in members:
  ids=ids|members.aliases.map(lambda s:bool(set(re.split(r'[;, ]+',str(s)))&BLOCKED))
 if (bad|ids).any():raise RuntimeError('Excluded specimen in comparative cohort: '+', '.join(members.loc[bad|ids,'sample']))

def raw_scores(matrix,y,mask,rule,min_rate=.8):
 """Both allele orientations, with the existing per-class 80% callability gate."""
 import numpy as np
 matrix=np.asarray(matrix,float);y=np.asarray(y,int);mask=np.asarray(mask,bool)
 called=np.isfinite(matrix);nsite=matrix.shape[1];ans={}
 predictions=[matrix>0,matrix<1] if rule=='carrier' else [matrix==1,matrix==0] if rule=='homozygote' else None
 if predictions is None:raise ValueError(rule)
 for state in [0,1]:
  group=mask&(y==state);n=int(group.sum());ncalled=called[group].sum(axis=0)
  ans['n_total_'+str(state)]=n;ans['n_called_'+str(state)]=ncalled
  ans['rate_'+str(state)]=ncalled/n if n else np.zeros(nsite)
 good=(ans['rate_0']>=min_rate)&(ans['rate_1']>=min_rate)
 for orientation,pred in zip(['ALT','REF'],predictions):
  sensitivity=np.divide((pred[mask&(y==1)]&called[mask&(y==1)]).sum(axis=0),ans['n_called_1'],out=np.full(nsite,np.nan),where=ans['n_called_1']>0)
  specificity=np.divide(((~pred[mask&(y==0)])&called[mask&(y==0)]).sum(axis=0),ans['n_called_0'],out=np.full(nsite,np.nan),where=ans['n_called_0']>0)
  ans[orientation]=np.where(good,(sensitivity+specificity)/2,np.nan)
 return ans

def calculate(matrix,y,geography,rule):
 import numpy as np
 y=np.asarray(y,int);geography=np.asarray(geography);all_samples=np.ones(len(y),bool)
 masks={'pooled':all_samples,'West_only':geography=='West','East_only':geography=='East',
        'West_positive_all_negative':((geography=='West')&(y==1))|(y==0)}
 raw={k:raw_scores(matrix,y,mask,rule) for k,mask in masks.items()}
 pooled=raw['pooled'];tied=np.isclose(pooled['ALT'],pooled['REF'],atol=1e-12,rtol=0)
 valid=np.isfinite(pooled['ALT'])&np.isfinite(pooled['REF'])&~tied
 orientation=np.where(pooled['ALT']>pooled['REF'],'ALT','REF')
 result={k:np.where(valid,np.where(orientation=='ALT',r['ALT'],r['REF']),np.nan) for k,r in raw.items()}
 result['equal_regions']=(result['West_only']+result['East_only'])/2
 equal_alt=(raw['West_only']['ALT']+raw['East_only']['ALT'])/2
 equal_ref=(raw['West_only']['REF']+raw['East_only']['REF'])/2
 # Export a separately labelled re-optimised sensitivity; never silently flip a regional score.
 equal_tie=np.isclose(equal_alt,equal_ref,atol=1e-12,rtol=0)
 best_equal=np.where(equal_tie,np.nan,np.maximum(equal_alt,equal_ref))
 return raw,result,np.where(valid,orientation,'unscored'),best_equal

def prepare(a):
 source=a.source;out=a.output
 paths=[source/'complete.json',source/'run.json',source/'membership.tsv',source/'SNPs.tsv',
 source/'individual_ALT_dosages.tsv.gz',source/'Ybv_all_SNP_BA.tsv.gz',
 a.parent/'worker.py',a.parent/'Ybv_gene_features.json',OLD,YF/'individual_ALT_dosages.tsv.gz',
 HO/'worker.py',HO/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz']
 for p in paths:
  if not p.is_file():raise RuntimeError('Missing documented input '+str(p))
 complete=json.loads((source/'complete.json').read_text())
 if complete.get('status')!='COMPLETE' or complete.get('span')!=[LO,HI]:raise RuntimeError('Completed full-span Yb calls required')
 if out.exists():
  if not a.overwrite:raise RuntimeError('Use --overwrite to update this diagnostic in place')
  cfg=json.loads((out/'run.json').read_text())
  if cfg.get('workflow')!='Ybv_geography_BA':raise RuntimeError('Output belongs to another workflow')
  job=str(json.loads((out/'jobs.json').read_text()).get('plot','')) if (out/'jobs.json').exists() else ''
  active=set(subprocess.check_output(['squeue','-h','-u',getpass.getuser(),'-o','%A'],text=True).split())
  if job and job in active:raise RuntimeError('Geography diagnostic still active')
 out.mkdir(parents=True,exist_ok=True);(out/'logs').mkdir(exist_ok=True)
 for p in [out/'complete.json',*out.glob('Ybv_BA_geography*.pdf')]:p.unlink(missing_ok=True)
 shutil.copy2(__file__,out/'worker.py');shutil.copy2(a.parent/'worker.py',out/'panel_helper.py')
 shutil.copy2(a.parent/'Ybv_gene_features.json',out/'Ybv_gene_features.json');shutil.copy2(HO/'worker.py',out/'Ho_helper.py')
 save(out/'run.json',dict(workflow='Ybv_geography_BA',source=str(source),parent=str(a.parent),
  input_hashes={str(p):sha(p) for p in paths},script_sha256=sha(__file__),spans=[LO,HI],
  cohort='Frozen completed BA cohort; sheet edits reserved for final rebuild',
  orientation='One pooled orientation per SNP/rule, held fixed across comparisons; pooled ties remain unscored',
  geographic_estimand='Equal-region mean requires both regions scoreable; does not adjust for race/relatedness',
  excluded_taxa=['numata','pardalinus'],no_new_calls=True,publication_markers_changed=False))
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--cpus-per-task=1','--mem=8G',
  '--job-name=Ybv_BAgeo','--output='+str(out/'logs/plot.%j.out'),'--error='+str(out/'logs/plot.%j.err'),
  '--wrap='+shlex.join([sys.executable,str(out/'worker.py'),'--run','--output',str(out)])],text=True).strip().split(';')[0]
 save(out/'jobs.json',dict(plot=job));print('JOB\t'+job+'\nOUTPUT\t'+str(out)+'\nPDF\t'+str(out/'Ybv_BA_geography_Ho_genotypes.pdf'))

def nucleotide_panel(helper,members):
 import numpy as np,pandas as pd
 old=pd.read_csv(OLD,sep='\t').sort_values('pos1');pos=old.pos1.astype(int).tolist()
 dos=pd.read_csv(YF/'individual_ALT_dosages.tsv.gz',sep='\t',index_col=0)
 if list(map(int,dos.columns))!=pos or not set(members['sample'])<=set(dos.index):raise RuntimeError('Publication marker identity mismatch')
 if len(pos)!=27:raise RuntimeError('Expected 27 unchanged publication markers')
 rows={};bread={}
 for state in [1,0]:
  group=members.loc[members.v.eq(state)];rows[state]=[]
  bread[state]=(f'Original pooled v{"+" if state else "-"} consensus (n={len(group)})',dos.loc[group['sample']].mean().to_numpy(float),'consensus')
  examples,_=helper.yb_example_groups(group)
  for geo,race in examples:
   ids=group.loc[group.geography.eq(geo)&group.race.eq(race),'sample']
   rows[state].append((f'{race} | v{state} | n={len(ids)}',dos.loc[ids].mean().to_numpy(float),'frequency'))
 focal=[]
 for side in ['left','right']:
  v=old[side+'_rosina_dosage'].to_numpy(float,copy=True);flip=old.rosina_base.eq(old.ref).to_numpy();v[flip]=1-v[flip]
  focal.append((side.title()+' (corrected)',v,'genotype'))
 return dict(sites=np.asarray(pos),rows=rows[1]+[bread[1]]+focal+[bread[0]]+rows[0],a=old.alt.tolist(),b=old.ref.tolist(),palette_type='ALT')

def figure(out,scores,members,counts,summary):
 import numpy as np,pandas as pd,matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.backends.backend_pdf import PdfPages
 from matplotlib.patches import Rectangle,ConnectionPatch
 plt.rcParams.update({'pdf.fonttype':42,'font.size':8,'font.family':'DejaVu Sans'})
 helper=load(out/'panel_helper.py','frozen_panel');ho=load(out/'Ho_helper.py','frozen_ho')
 letters=nucleotide_panel(helper,members);native=[g for g in helper.genes() if g['chrom']==CHROM]
 ivory=json.loads((out/'Ybv_gene_features.json').read_text())
 pairs=ho.normalize(pd.read_csv(HO/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz',sep='\t'),'Ivory_Yb')
 windows={q:ho.windows(pairs,LO,HI,'bp',10000,q) for q in [.95,.99]}
 for q,w in windows.items():w.to_csv(out/f'Ybv_Ho_10kb_step2500_P{q}.tsv.gz',sep='\t',index=False)
 ymax=max(float(w[['left_Ho','right_Ho']].max().max()) for w in windows.values())
 ylim=min(1.,max(.3,float(np.ceil(ymax*1.1*10)/10)))
 pdf=out/'Ybv_BA_geography_Ho_genotypes.pdf'
 colours=['#333333','#1764AB','#C7503B','#71439A','#3B8064']
 with PdfPages(pdf) as pages:
  for rule in ['carrier','homozygote']:
   fig=plt.figure(figsize=(14,17));left=.35;width=.62
   fig.text(.025,.975,'Ivory / Yb(v): geographic BA sensitivity - '+rule+' rule',fontsize=14,weight='bold')
   fig.text(left,.952,'Fixed allele orientation from pooled BA; identical calls and genomic span',fontsize=8)
   for i,(comparison,col) in enumerate(zip(COMPARISONS,colours)):
    ax=fig.add_axes([left,.832-i*.101,width,.075]);s=scores.loc[scores.rule.eq(rule)&scores.comparison.eq(comparison)]
    ax.scatter(s.pos,s.balanced_accuracy,s=5,c=col,alpha=.65,linewidths=0)
    ax.axhline(.5,color='#aaaaaa',ls=':',lw=.5);ax.axhline(.9,color='#B15422',ls='--',lw=.6)
    ax.set_xlim(LO,HI);ax.set_ylim(0,1.02);ax.set_ylabel('BA');ax.tick_params(axis='x',labelbottom=False);ax.spines[['top','right']].set_visible(False)
    cc=counts[comparison];ntext='two within-region scores, equally weighted' if comparison=='equal_regions' else f'v+ n={cc[1]}, v- n={cc[0]}'
    ax.set_title(LABELS[comparison]+' | '+ntext,fontsize=8,loc='left',pad=3)
   ax=fig.add_axes([left,.314,width,.095])
   for q,ls in [(.95,'-'),(.99,'--')]:
    w=windows[q]
    for side,col in [('left','#1764AB'),('right','#D04B35')]:ax.plot(w.mid,w[side+'_Ho'],c=col,ls=ls,lw=1.1,label=side.title()+f' P>={q}')
   ax.set_xlim(LO,HI);ax.set_ylim(0,ylim);ax.set_ylabel('Heterozygous fraction');ax.legend(fontsize=6,ncol=2,loc='upper right');ax.tick_params(axis='x',labelbottom=False);ax.spines[['top','right']].set_visible(False)
   ann=fig.add_axes([left,.245,width,.06]);ann.set_xlim(LO,HI);ann.set_ylim(-.2,4.5);ann.axhline(0,c='#999999',lw=.5)
   for label,a,b in helper.INTERVALS['Yb(v)']:
    ann.add_patch(Rectangle((a,.15),b-a,.32,facecolor='#777777'));ann.text((a+b)/2,.6,label,ha='center',fontsize=7)
   c=ivory['cortex'];gs=[g for g in native if not (g['start']==c['start'] and g['end']==c['end'])]
   gs+=[dict(chrom=CHROM,label='cortex',**{k:c[k] for k in ['start','end','strand','exons']})]
   iv=ivory['ivory'];gs+=[dict(chrom=CHROM,start=min(x[0] for x in iv['exons']),end=max(x[1] for x in iv['exons']),strand='+',label='ivory (continues)',exons=iv['exons'])]
   lanes=[LO-1]*3
   for g in sorted(gs,key=lambda z:z['start']):
    lane=next((i for i,end in enumerate(lanes) if g['start']>end),int(np.argmin(lanes)))
    lanes[lane]=g['end']+(HI-LO)*.04;yy=1.5+lane*.95;a=max(LO,g['start']);b=min(HI,g['end'])
    if b<a:continue
    col='#267B82' if g['label'].startswith('ivory') else '#444444'
    ann.plot([a,b],[yy,yy],c=col,lw=.65)
    for a1,b1 in g['exons']:
     if b1<LO or a1>HI:continue
     ann.add_patch(Rectangle((max(LO,a1),yy-.11),min(HI,b1)-max(LO,a1)+1,.22,facecolor=col,lw=0))
    x1,x2=(a,b) if g['strand']=='+' else (b,a)
    ann.annotate('',xy=(x2,yy),xytext=(x1,yy),arrowprops=dict(arrowstyle='->',color=col,lw=.6))
    ann.text((a+b)/2,yy+.2,g['label'],fontsize=6.5,ha='center',color=col)
   ticks=np.linspace(LO,HI,5);ann.set_xticks(ticks,[f'{v/1e6:.3f}' for v in ticks]);ann.set_yticks([]);ann.set_xlabel('Rosina chromosome 15 (Mb)',fontsize=7)
   for sp in ann.spines.values():sp.set_visible(False)
   heat=fig.add_axes([left,.047,width,.17]);helper.draw_letters(heat,letters)
   for j,pos in enumerate(letters['sites']):fig.add_artist(ConnectionPatch(xyA=(int(pos),0),coordsA=ann.transData,xyB=(j,-.5),coordsB=heat.transData,color='#BFC3C7',lw=.35,alpha=.5,clip_on=False))
   fig.text(.025,.017,'BA: >=80% called in each phenotype class within each tested cohort. Missing/unscoreable values are not zero.\nHo and the original 27-SNP pooled nucleotide panel are held constant for comparison; no marker reselection or focal inference.\nEqual-region BA requires both regions; geography bins retain the documented population assignments, not a literal Andes-side boundary.',fontsize=6.5)
   pages.savefig(fig);plt.close(fig)
 # Compact matched-SNP changes so shifts can be assessed without eyeballing peaks.
 delta_pdf=out/'Ybv_BA_geography_shifts.pdf'
 with PdfPages(delta_pdf) as pages:
  for rule in ['carrier','homozygote']:
   pivot=scores.loc[scores.rule.eq(rule)].pivot(index='pos',columns='comparison',values='balanced_accuracy')
   fig,axes=plt.subplots(4,2,figsize=(11,12),gridspec_kw={'width_ratios':[2.5,1]})
   fig.subplots_adjust(left=.11,right=.97,top=.91,bottom=.10,hspace=.65,wspace=.3)
   for i,comparison in enumerate(COMPARISONS[1:]):
    x=pivot['pooled'];y=pivot[comparison];good=x.notna()&y.notna();delta=y[good]-x[good]
    axes[i,0].scatter(pivot.index[good]/1e6,delta,s=5,c=colours[i+1],alpha=.65,linewidths=0);axes[i,0].axhline(0,c='#999999',lw=.5)
    axes[i,0].set_xlim(LO/1e6,HI/1e6);axes[i,0].set_ylim(-1,1);axes[i,0].set_ylabel('BA minus pooled BA');axes[i,0].set_title(LABELS[comparison],loc='left',fontsize=8)
    axes[i,1].scatter(x[good],y[good],s=5,c=colours[i+1],alpha=.4,linewidths=0);axes[i,1].plot([0,1],[0,1],c='#999999',lw=.5)
    axes[i,1].set_xlim(0,1);axes[i,1].set_ylim(0,1);axes[i,1].set_title(f'Matched scoreable SNPs: {int(good.sum())}',fontsize=8)
    for ax in axes[i]:ax.spines[['top','right']].set_visible(False)
   axes[-1,0].set_xlabel('Chromosome 15 (Mb)');axes[-1,1].set_xlabel('Pooled BA')
   fig.suptitle('Yb(v) BA shift with geographic comparison - '+rule+' rule',fontsize=13)
   fig.text(.11,.035,'Same calls and fixed pooled allele orientation; regional BA can fall below 0.5 when the association reverses.\nEqual-region scores require scoreable West and East classes. Race composition and relatedness are not adjusted.',fontsize=8)
   pages.savefig(fig);plt.close(fig)
 return [str(pdf),str(delta_pdf)]

def run(out):
 import numpy as np,pandas as pd
 cfg=json.loads((out/'run.json').read_text());source=Path(cfg['source'])
 for p,h in cfg['input_hashes'].items():
  if sha(p)!=h:raise RuntimeError('Input changed since submission: '+p)
 complete=json.loads((source/'complete.json').read_text())
 members=pd.read_csv(source/'membership.tsv',sep='\t').fillna('');members['v']=pd.to_numeric(members.v,errors='raise')
 ensure_allowed(members)
 if not members.v.isin([0,1]).all() or not members.geography.isin(['West','East']).all():raise RuntimeError('Resolve missing phenotype/geography before comparison')
 sites=pd.read_csv(source/'SNPs.tsv',sep='\t');dos=pd.read_csv(source/'individual_ALT_dosages.tsv.gz',sep='\t',index_col=0)
 if members['sample'].duplicated().any() or dos.index.duplicated().any() or sites.pos.duplicated().any():raise RuntimeError('Duplicate identity/position')
 if set(dos.index)!=set(members['sample']) or list(map(int,dos.columns))!=sites.pos.astype(int).tolist():raise RuntimeError('Call matrix identity mismatch')
 if len(members)!=complete['individuals'] or len(sites)!=complete['SNPs']:raise RuntimeError('Completed call-set dimensions differ')
 matrix=dos.loc[members['sample']].to_numpy(float);y=members.v.to_numpy(int);geo=members.geography.to_numpy()
 if np.any(np.isfinite(matrix)&~np.isin(matrix,[0,.5,1])):raise RuntimeError('Invalid diploid dosage')
 counts={k:{state:int(sum(mask&(y==state))) for state in [0,1]} for k,mask in {'pooled':np.ones(len(y),bool),'West_only':geo=='West','East_only':geo=='East','West_positive_all_negative':((geo=='West')&(y==1))|(y==0)}.items()}
 counts['equal_regions']={0:None,1:None}
 members.to_csv(out/'membership.tsv',sep='\t',index=False)
 members.groupby(['geography','race','v']).size().reset_index(name='individuals').to_csv(out/'cohort_counts_by_race.tsv',sep='\t',index=False)
 pd.DataFrame([dict(comparison=k,negative_n=v[0],positive_n=v[1]) for k,v in counts.items()]).to_csv(out/'cohort_counts.tsv',sep='\t',index=False)
 prior=pd.read_csv(source/'Ybv_all_SNP_BA.tsv.gz',sep='\t');tables=[];shifts=[];summary=[]
 for rule in ['carrier','homozygote']:
  raw,values,orientation,best_equal=calculate(matrix,y,geo,rule)
  expected=prior.loc[prior.rule.eq(rule)].set_index('pos').reindex(sites.pos).balanced_accuracy.to_numpy(float)
  if not np.allclose(values['pooled'],expected,equal_nan=True,atol=1e-12,rtol=0):raise RuntimeError('Recomputed pooled BA differs from completed scan')
  for comp in COMPARISONS:
   tab=sites.copy();tab['rule']=rule;tab['comparison']=comp;tab['balanced_accuracy']=values[comp];tab['pooled_orientation']=orientation
   tab['delta_vs_pooled']=values[comp]-values['pooled'];tab['equal_region_best_orientation_BA']=best_equal if comp=='equal_regions' else np.nan
   if comp in raw:
    for name in ['ALT','REF','n_called_0','n_called_1','rate_0','rate_1']:tab[name]=raw[comp][name]
   tab['status']=np.where(np.isfinite(values[comp]),'scored',np.where(orientation=='unscored','pooled_tie_or_insufficient_callability','regional_classes_or_callability_unscoreable'))
   tables.append(tab)
   both=np.isfinite(values[comp])&np.isfinite(values['pooled']);delta=values[comp][both]-values['pooled'][both]
   summary.append(dict(rule=rule,comparison=comp,scoreable=int(np.isfinite(values[comp]).sum()),matched_SNPs=int(both.sum()),
    median_shift=float(np.median(delta)) if len(delta) else None,median_absolute_shift=float(np.median(np.abs(delta))) if len(delta) else None,
    maximum_absolute_shift=float(np.max(np.abs(delta))) if len(delta) else None,BA_ge_090=int(np.sum(values[comp]>=.9)),
    gained_ge_090=int(np.sum(both&(values[comp]>=.9)&(values['pooled']<.9))),lost_ge_090=int(np.sum(both&(values[comp]<.9)&(values['pooled']>=.9)))))
   if comp!='pooled':shifts.append(tab.loc[both].assign(absolute_shift=np.abs(tab.loc[both,'delta_vs_pooled'])).nlargest(25,'absolute_shift'))
 scores=pd.concat(tables,ignore_index=True);scores.to_csv(out/'Ybv_geography_BA_all_SNPs.tsv.gz',sep='\t',index=False)
 pd.concat(shifts,ignore_index=True).to_csv(out/'largest_BA_shifts.tsv',sep='\t',index=False)
 pd.DataFrame(summary).to_csv(out/'BA_shift_summary.tsv',sep='\t',index=False)
 # Also export every original marker, including missing/non-biallelic markers, without reselection.
 marker=pd.read_csv(OLD,sep='\t')[['pos1','ref','alt']].rename(columns={'pos1':'pos'})
 combinations=pd.DataFrame([dict(rule=rule,comparison=comp) for rule in ['carrier','homozygote'] for comp in COMPARISONS])
 marker=marker.merge(combinations,how='cross').merge(scores,on=['pos','ref','alt','rule','comparison'],how='left',validate='one_to_one')
 marker['status']=marker.status.fillna('original_marker_not_in_new_biallelic_SNP_scan')
 marker.to_csv(out/'original_27_marker_BA.tsv',sep='\t',index=False)
 pdfs=figure(out,scores,members,counts,summary)
 save(out/'complete.json',dict(status='COMPLETE',PDFs=pdfs,SNPs=len(sites),individuals=len(members),cohorts=counts,
  pooled_reproduces_completed_scan=True,summary=summary,no_new_calls=True,focal_calls_changed=False,
  publication_markers_reselected=False,genomic_span=[LO,HI],cohort_sheet_edits_applied=False,
  geography_adjustment='Separate West/East BA and equal-region average; pooled allele orientation fixed',
  geography_caveat='Population bins retained from membership.tsv; no race/kinship adjustment'))
 print('STATUS\tCOMPLETE')
 for p in pdfs:print('PDF\t'+p)
 print('SUMMARY\t'+str(out/'BA_shift_summary.tsv'))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=PARENT/'Ybv_BA_broad');p.add_argument('--parent',type=Path,default=PARENT)
 p.add_argument('--output',type=Path,default=PARENT/'Ybv_BA_geography');p.add_argument('--submit',action='store_true');p.add_argument('--run',action='store_true');p.add_argument('--overwrite',action='store_true');a=p.parse_args()
 if a.submit:prepare(a)
 elif a.run:run(a.output)
 else:p.error('Use --submit')
if __name__=='__main__':main()
