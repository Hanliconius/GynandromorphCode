#!/usr/bin/env python3
"""Optix pooled positive/negative controls, separate PK groups and corrected focal dosage.
Uses the uniform-call renderer/scorer; new user annotations are recorded in
this output only. Existing data and the independent-genotype PDF stay intact.
"""
import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plot_optix_uniform_genoplots as base

ORIGINAL_MANIFEST=base.phenotype_manifest


def annotated_manifest(cfg,root,pk_path,overrides=None):
    meta,sources=ORIGINAL_MANIFEST(cfg,root,pk_path,overrides)
    meta['Yb_dorsal']=np.nan
    meta['Yb_ventral']=np.nan
    meta['Yb_annotation_source']=''
    race=meta['group'].str.lower()
    selected=race.isin(['cythera','vulcanus']) & ~meta.is_PK
    meta.loc[selected,['B','D','R']]=[1,0,0]
    for mod in base.MODULES:meta.loc[selected,mod+'_selection_eligible']=True
    meta.loc[selected,'phenotype_source']='Joe Hanly, explicit race phenotype confirmation, 2026-10-16'
    vul=race.eq('vulcanus') & ~meta.is_PK
    meta.loc[vul,['Yb_dorsal','Yb_ventral']]=[0,1]
    meta.loc[vul,'Yb_annotation_source']='Joe Hanly, 2026-10-16: dorsal yellow band absent; ventral yellow band present'
    if meta[['B','D','R']].isna().any().any():
        missing=meta.loc[meta[['B','D','R']].isna().any(axis=1),'sample'].tolist()
        raise RuntimeError('Unresolved B/D/R phenotypes: '+','.join(missing))
    return meta,sources+['User-confirmed cythera/vulcanus B1 D0 R0; vulcanus Yb dorsal0 ventral1, 2026-10-16']


def frequency_groups(meta,names,module):
    """One positive and one negative non-PK cohort; unchanged PK groups."""
    if meta['sample'].duplicated().any() or len(names)!=len(set(names)):
        raise RuntimeError('Duplicated individual in pooling input')
    lookup={name:i for i,name in enumerate(names)}
    eligible=meta.loc[~meta.is_PK & meta[module+'_selection_eligible']]
    trait=base.TRAITS[module]
    controls={}
    for value in (1,0):
        rows=eligible.loc[eligible[trait].eq(value)]
        if rows.empty:raise RuntimeError('Empty pooled control class')
        controls[value]=dict(group=module+('+' if value else '−'),is_PK=False,
            samples=rows['sample'].tolist(),indices=[lookup[n] for n in rows['sample']],
            n=len(rows),label=f'{module}{"+" if value else "−"} controls | n={len(rows)}',
            trait=module,phenotype=value)
    pk=[]
    for key,rows in meta.loc[meta.is_PK].groupby(['group','B','D','R'],sort=True,dropna=False):
        pop,b,d,r=key
        if any(pd.isna(x) for x in (b,d,r)):raise RuntimeError('Unknown PK phenotype')
        pk.append(dict(group=pop,B=int(b),D=int(d),R=int(r),is_PK=True,
            samples=rows['sample'].tolist(),indices=[lookup[n] for n in rows['sample']],n=len(rows),
            label=f'{pop.replace("PK_","PK ")} | B{int(b)} D{int(d)} R{int(r)} | n={len(rows)}'))
    groups=[controls[1],*pk,controls[0]]
    members=[name for g in groups for name in g['samples']]
    if len(members)!=len(set(members)):raise RuntimeError('An individual entered more than one row')
    return groups


def group_frequency(values,min_rate=.8):
    """Equal-individual genotype dosage mean; missing GT is never REF."""
    a=np.asarray(values,dtype=float)
    n=int(np.isfinite(a).sum())
    if not len(a) or n/len(a)<min_rate:return np.nan,n
    return float(np.nanmean(a)),n


def frequency_matrix(t,dosages,recall,groups):
    matrix=np.full((2+len(groups),len(t)),np.nan)
    called=np.zeros((len(groups),len(t)),dtype=int)
    for j,r in enumerate(t.to_dict('records')):
        d=dosages[int(r['site_index'])].copy()
        if not r['pattern_is_alt']:d=1-d
        focal=recall[tuple(r[k] for k in base.KEY)]
        for i,side in enumerate(('left','right')):
            dose=sum(int(x) for x in focal['mixture_'+side+'_GT'].split('/'))/2
            matrix[i,j]=dose if r['pattern_is_alt'] else 1-dose
        for i,g in enumerate(groups):
            matrix[2+i,j],called[i,j]=group_frequency(d[g['indices']])
    return matrix,called


def draw_frequency(t,dosages,recall,meta,names,cohort,module,cutoff,interval_only,out,pdf,counts):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Rectangle,ConnectionPatch
    groups=frequency_groups(meta,names,module)
    matrix,called=frequency_matrix(t,dosages,recall,groups)
    # Raw matrix order is L,R,+,PK...,-. Display +,L,R,PK...,-.
    order=[2,0,1,*range(3,len(groups)+2)]
    display_matrix=matrix[order]
    labels=[groups[0]['label'],'TNLeft corrected dosage','TNRight corrected dosage']+[g['label'] for g in groups[1:]]
    nrows=len(labels)
    fig=plt.figure(figsize=(12,max(7,3+nrows*.31)))
    fig.suptitle(f'optix — {module} allele frequencies',x=.34,ha='left',fontsize=14)
    fig.text(.34,.93,f'BA ≥ {cutoff:.2f} | '+('labelled intervals' if interval_only else 'full window')+
             f' | {counts["focal_callable"]} callable SNPs; {len(t)} displayed',fontsize=9)
    locus=fig.add_axes([.34,.82,.59,.065]);heat=fig.add_axes([.34,.16,.59,.55])
    if len(t):
        cmap=LinearSegmentedColormap.from_list('frequency',['#D7EDF5','#ECCB63','#662A83'])
        cmap.set_bad('#C8C8C8')
        im=heat.pcolormesh(np.arange(len(t)+1)-.5,np.arange(nrows+1)-.5,display_matrix,
                          cmap=cmap,vmin=0,vmax=1,shading='flat',rasterized=False)
        heat.set_xlim(-.5,len(t)-.5);heat.set_ylim(nrows-.5,-.5)
        heat.set_yticks(range(nrows),labels,fontsize=8)
        for boundary in (.5,2.5,nrows-1.5):heat.axhline(boundary,color='white',lw=2)
        bar=fig.colorbar(im,cax=fig.add_axes([.948,.25,.012,.30]))
        bar.set_ticks([0,.5,1]);bar.set_label('Trait-associated allele frequency / focal dosage',fontsize=8)
        for j,r in enumerate(t.itertuples()):
            fig.add_artist(ConnectionPatch(xyA=(r.pos/1e6,-.12),coordsA=locus.transData,
                xyB=(j,-.5),coordsB=heat.transData,color='#929aa3',lw=.35,alpha=.16,clip_on=False))
    else:
        heat.text(.5,.5,'No SNPs pass these filters',transform=heat.transAxes,ha='center');heat.set_yticks([])
    heat.set_xticks([]);heat.tick_params(axis='y',length=0)
    locus.set_xlim(1.06,1.24);locus.set_ylim(-.2,1.2);locus.axhline(0,color='#666',lw=.8)
    for mod in base.MODULES:
        for i,(lo,hi) in enumerate(base.INTERVALS[mod]):
            colour='#555' if mod==module else '#c8c8c8'
            locus.add_patch(Rectangle((lo/1e6,.1),(hi-lo)/1e6,.25,color=colour,lw=0))
            locus.text((lo+hi)/2e6,.63 if not (mod=='Dennis' and i==1) else 1.05,
                {'Band':'B','Ray':'R','Dennis':'D'}[mod]+('' if mod=='Band' else str(i+1)),ha='center',fontsize=8,color=colour)
    locus.set_yticks([]);locus.set_xlabel('Rosina chromosome 18 (NC_140406.1), Mb; analysis windows',fontsize=8)
    for ax in (heat,locus):
        for spine in ax.spines.values():spine.set_visible(False)
    n0=int(np.sum(cohort[1]==0));n1=int(np.sum(cohort[1]==1))
    fig.text(.34,.055,f'Trait scoring: non-PK controls only (absent n={n0}; present n={n1}); DP≥2, GQ≥20.\n'
        'Group frequency = mean called diploid allele dosage; ≥80% of group callable per SNP, otherwise grey.\n'
        'TNLeft/Right = corrected genotype dosage (pair posterior ≥0.95), not mixed-read frequency.\n'
        'Non-PK controls pooled once per individual; PK colour and B/D/R groups remain separate.',fontsize=8)
    stem=f'{module}_BA{int(cutoff*100)}_{"intervals" if interval_only else "full"}'
    pd.DataFrame(display_matrix,index=labels,columns=[f'{r.chrom}:{r.pos}:{r.ref}:{r.alt}' for r in t.itertuples()]).to_csv(
        out/(stem+'_frequencies.tsv'),sep='\t',na_rep='NA')
    rows=[]
    for i,g in enumerate(groups):
        for j,r in enumerate(t.itertuples()):
            rows.append(dict(group=g['label'],chrom=r.chrom,pos=r.pos,ref=r.ref,alt=r.alt,
                n_called=int(called[i,j]),n_total=g['n'],called_chromosomes=2*int(called[i,j]),
                frequency=None if np.isnan(matrix[i+2,j]) else float(matrix[i+2,j])))
    pd.DataFrame(rows,columns=['group','chrom','pos','ref','alt','n_called','n_total','called_chromosomes','frequency']).to_csv(
        out/(stem+'_group_callability.tsv.gz'),sep='\t',index=False,na_rep='NA')
    # Store full membership once; no inferred phenotype or read-depth weighting.
    (out/(module+'_frequency_group_membership.json')).write_text(json.dumps(groups,indent=2)+'\n')
    pdf.savefig(fig);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--pk-catalog',type=Path,required=True)
    p.add_argument('--root',type=Path,default=base.ROOT)
    p.add_argument('--max-sites',type=int,default=80)
    p.add_argument('--phenotype-overrides',type=Path)
    a=p.parse_args()
    base.phenotype_manifest=annotated_manifest
    base.draw=draw_frequency
    messages=io.StringIO()
    with contextlib.redirect_stdout(messages):base.main()
    (a.output_dir/'scoring.log').write_text(messages.getvalue())
    old=a.output_dir/'Optix_uniform_Band_Dennis_Ray_genoplots.pdf'
    dest=a.output_dir/'Optix_Band_Dennis_Ray_pooled_controls_PK.pdf'
    old.rename(dest)
    auditpath=a.output_dir/'audit.json';audit=json.loads(auditpath.read_text())
    audit.update(pdf=str(dest),pdf_sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),
        display='pooled trait-positive/negative non-PK allele frequency; separate PK groups; corrected focal dosage',
        frequency='sum associated allele copies / (2 * number of called diploid individuals)',
        min_group_call_rate=.8,read_depth_weighted=False,
        grouping=['non-PK eligible positive controls','non-PK eligible negative controls','PK dark/pale by B/D/R'],
        duplicate_individual_contributions=False,
        user_confirmed_BDR={'cythera':[1,0,0],'vulcanus':[1,0,0]},
        vulcanus_Yb={'dorsal':0,'ventral':1,'used_for_Optix_scoring':False})
    auditpath.write_text(json.dumps(audit,indent=2)+'\n')
    print('PDF\t'+str(dest));print('AUDIT\t'+str(auditpath))

if __name__=='__main__':main()
