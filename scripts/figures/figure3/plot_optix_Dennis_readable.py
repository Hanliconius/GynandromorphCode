#!/usr/bin/env python3
"""Plot frozen Dennis audit calls with explicit row alignment and REF/ALT scale."""
import argparse, hashlib, json, shlex, shutil, subprocess, sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis')

def plot(a):
 import numpy as np
 import pandas as pd
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.colors import LinearSegmentedColormap
 from matplotlib.backends.backend_pdf import PdfPages
 plt.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Sans','font.size':8})
 a.output.mkdir(parents=True,exist_ok=False)
 meta=pd.read_csv(a.audit/'sample_identities.tsv',sep='\t').fillna('')
 sites=pd.read_csv(a.audit/'Dennis_site_summary.tsv',sep='\t')
 sites=sites[sites.source.eq('legacy')].sort_values('pos').copy()
 if sites.duplicated(['pos','ref','alt']).any():raise RuntimeError('Duplicate audited variant')
 sites['strong']=(sites.D1_n_called.ge(30)&sites.D0_n_called.ge(30)&(sites.D1_ALT_AF-sites.D0_ALT_AF).ge(.8))
 calls=pd.read_csv(a.audit/'Dennis_individual_raw_calls.tsv.gz',sep='\t')
 calls=calls[calls.source.eq('legacy')]
 races=pd.read_csv(a.audit/'race_site_ALT_frequencies.tsv.gz',sep='\t')
 races=races[races.source.eq('legacy')]
 focus=meta[meta.group.isin(['meriana_unquoted','meriana_quoted'])].sort_values(['group','sample'],ascending=[False,True])
 # Named group order is explicit; no dependence on alphabetic label sorting.
 focus=pd.concat([focus[focus.group.eq('meriana_unquoted')],focus[focus.group.eq('meriana_quoted')]])
 rows=[]
 for r in focus.itertuples():
  country=' '.join(str(r.Country).split())
  label=('meriana' if r.group=='meriana_unquoted' else "‘meriana’")+'  |  '+r.sample+'  |  '+country+'  |  D+'
  rows.append((label,'individual',r.sample))
 for group,label in [('rosina','rosina'),('melpomene_melpomene','melpomene'),('besckei','besckei')]:
  ids=set(meta.loc[meta.group.eq(group),'sample'])&set(calls['sample'])
  rows.append((f'{label}  |  n={len(ids)}','group',group))
 for taxon in sorted(races.taxon.unique()):
  if 'timareta' in taxon.lower():
   n=int(races.loc[races.taxon.eq(taxon),'n_total'].max())
   rows.append((taxon.replace('Heliconius ','').replace('ssp. ','')+f'  |  n={n}','race',taxon))
 lookup=calls.pivot(index='sample',columns='pos',values='qualified_alt_dosage')/2
 cmap=LinearSegmentedColormap.from_list('ref_alt',['#D7EDF5','#ECCB63','#662A83']);cmap.set_bad('#D0D0D0')
 pages=[]
 pdfpath=a.output/'Optix_Dennis_readable_REF_ALT.pdf'
 with PdfPages(pdfpath) as pdf:
  for short in [True,False]:
   ss=sites[sites.strong] if short else sites
   if ss.empty:raise RuntimeError('No qualifying markers')
   positions=ss.pos.tolist();data=[];counts=[]
   for label,kind,key in rows:
    if kind=='individual':
     vals=lookup.reindex(index=[key],columns=positions);v=vals.iloc[0].to_numpy();nc=vals.notna().sum().to_numpy()
    elif kind=='group':
     ids=meta.loc[meta.group.eq(key),'sample'];vals=lookup.reindex(index=ids,columns=positions);v=vals.mean().to_numpy();nc=vals.notna().sum().to_numpy()
    else:
     rr=races[races.taxon.eq(key)].set_index('pos').reindex(positions);v=rr.ALT_AF.to_numpy();nc=rr.n_called.fillna(0).to_numpy()
    data.append(v);counts.extend(dict(view='associated' if short else 'all',row=label,pos=int(p),n_called=int(n),ALT_frequency=float(x) if np.isfinite(x) else None) for p,n,x in zip(positions,nc,v))
   fig=plt.figure(figsize=(12,max(6,2.8+.29*len(rows))))
   gs=fig.add_gridspec(2,1,height_ratios=[1,len(rows)],hspace=.08,left=.40,right=.98,top=.86,bottom=.20)
   score=fig.add_subplot(gs[0]);ax=fig.add_subplot(gs[1],sharex=score)
   score.bar(np.arange(len(ss)),ss.D1_ALT_AF-ss.D0_ALT_AF,width=1,color='#555555')
   score.set_ylim(0,1);score.set_yticks([0,1]);score.set_ylabel('Δ ALT AF\nD+ − D−',rotation=0,ha='right',va='center',fontsize=7)
   score.tick_params(axis='x',bottom=False,labelbottom=False);score.spines[['top','right']].set_visible(False)
   im=ax.imshow(np.ma.masked_invalid(np.asarray(data,dtype=float)),vmin=0,vmax=1,cmap=cmap,aspect='auto',interpolation='nearest')
   ax.set_yticks(range(len(rows)),[r[0] for r in rows],fontsize=8);ax.tick_params(axis='y',length=0,pad=9)
   for y in np.arange(.5,len(rows),1):ax.axhline(y,color='white',lw=.6)
   for y in [2.5,6.5]:ax.axhline(y,color='#333333',lw=1)
   split=sum(ss.interval.eq('Dennis1'))
   if 0<split<len(ss):
    for axis in [ax,score]:axis.axvline(split-.5,color='black',lw=1.5)
   if short:
    ax.set_xticks(range(len(ss)),[f'{r.pos:,}\n{r.ref}→{r.alt}' for r in ss.itertuples()],rotation=90,fontsize=6)
   else:
    ticks=np.unique(np.linspace(0,len(ss)-1,10).astype(int));ax.set_xticks(ticks,[f'{positions[i]:,}' for i in ticks],fontsize=7)
   ax.set_xlabel('Chromosome 18 position (Rosina); equally spaced SNP columns. Black divider: Dennis1 / Dennis2.',fontsize=8)
   fig.suptitle('Dennis-associated variants' if short else 'Dennis intervals: all audited biallelic SNPs',x=.40,ha='left',fontsize=13,fontweight='bold')
   fig.text(.40,.91,f'{len(ss)} SNPs | original legacy calls | individual meriana rows; comparative rows are frequencies',fontsize=8)
   cb=fig.colorbar(im,ax=[score,ax],orientation='horizontal',fraction=.035,pad=.14,shrink=.65,ticks=[0,.5,1]);cb.set_label('ALT allele frequency: 0 = REF / REF     0.5 = REF / ALT     1 = ALT / ALT (individual rows)')
   fig.text(.03,.045,'Grey = no qualified call. DP ≥ 2; GQ ≥ 20, or explicitly audited PL-derived likelihood Q ≥ 20.\nCompact page: ALT AF(D+) − ALT AF(D−) ≥ 0.8; ≥30 called individuals in each phenotype group.\nComparative frequencies average individual dosages; n labels show available individuals, site-specific n is exported.\nQuoted and unquoted meriana are separate metadata identities. REF/ALT does not imply ancestral/derived or causal state.',fontsize=7)
   pdf.savefig(fig);plt.close(fig)
   pd.DataFrame(counts).to_csv(a.output/('associated_site_frequencies.tsv' if short else 'all_site_frequencies.tsv'),sep='\t',index=False)
   pages.append({'view':'associated' if short else 'all','SNPs':len(ss),'rows':len(rows)})
 sites.to_csv(a.output/'marker_selection.tsv',sep='\t',index=False)
 (a.output/'complete.json').write_text(json.dumps({'status':'COMPLETE','source':str(a.audit),'call_source':'legacy','scale':'0=REF, 1=ALT','pages':pages,'PDF':str(pdfpath)},indent=2)+'\n')
 print('PDF\t'+str(pdfpath))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--audit',type=Path,default=ROOT/'optix_meriana_Dennis_audit_v2');p.add_argument('--output',type=Path,default=ROOT/'optix_Dennis_readable_v1');p.add_argument('--submit',action='store_true');a=p.parse_args()
 if a.submit:
  stage=ROOT/'optix_meriana_audit_code';stage.mkdir(parents=True,exist_ok=True)
  frozen=stage/('readable_'+hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:12]+'.py');shutil.copy2(__file__,frozen)
  cmd=shlex.join([sys.executable,str(frozen),'--audit',str(a.audit),'--output',str(a.output)])
  job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=6G','--cpus-per-task=1','--job-name=Optix_Dread','--output='+str(stage/'readable.%j.out'),'--error='+str(stage/'readable.%j.err'),'--wrap='+cmd],text=True).strip();print('JOB\t'+job)
 else:plot(a)
if __name__=='__main__':main()
