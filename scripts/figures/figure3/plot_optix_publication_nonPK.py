"""Compact vector Optix panels; PK excluded from scoring and control summaries."""
import json
import numpy as np
import pandas as pd
import plot_optix_uniform_genoplots as base
import plot_optix_frequency_genoplots as freq

def frequency_groups(meta,names,module):
    if meta['sample'].duplicated().any():raise RuntimeError('Duplicate sample')
    lookup={n:i for i,n in enumerate(names)}
    t=meta.loc[~meta.is_PK & meta[module+'_selection_eligible']]
    groups=[]
    for value in (1,0):
        r=t.loc[t[base.TRAITS[module]].eq(value)]
        if r.empty:raise RuntimeError('Empty phenotype class')
        groups.append(dict(label=f'{module}{"+" if value else "−"} (n = {len(r)})',
            samples=r['sample'].tolist(),indices=[lookup[n] for n in r['sample']],n=len(r),n_PK=0,phenotype=value))
    return groups

def render(t,matrix,labels,module,cutoff,scope,pdf):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Rectangle,ConnectionPatch
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':7,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.5})
    n=len(labels);height=max(2.45,1.35+n*.16)
    fig=plt.figure(figsize=(7.2,height))
    fig.text(.025,.94,{'Band':'A','Dennis':'B','Ray':'C'}[module],weight='bold',fontsize=10)
    fig.text(.23,.94,module,fontsize=9,weight='bold')
    fig.text(.96,.94,f'BA ≥ {cutoff:.2f}'+(' · labelled intervals' if scope else ''),ha='right',fontsize=6.5)
    locus=fig.add_axes([.23,.73,.70,.09]);heat=fig.add_axes([.23,.25,.70,.39])
    cmap=LinearSegmentedColormap.from_list('associated',['#D7EDF5','#ECCB63','#662A83']);cmap.set_bad('#C8C8C8')
    if len(t):
        im=heat.pcolormesh(np.arange(len(t)+1)-.5,np.arange(n+1)-.5,matrix,cmap=cmap,vmin=0,vmax=1,shading='flat')
        heat.set_xlim(-.5,len(t)-.5)
        for j,r in enumerate(t.itertuples()):
            fig.add_artist(ConnectionPatch(xyA=(r.pos/1e6,0),coordsA=locus.transData,xyB=(j,-.5),coordsB=heat.transData,color='#9ca2a8',lw=.25,alpha=.22,clip_on=False))
        cb=fig.colorbar(im,cax=fig.add_axes([.69,.105,.23,.025]),orientation='horizontal')
        cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=6,length=2,pad=1);cb.outline.set_visible(False)
        cb.set_label('Associated allele frequency / focal dosage',fontsize=6,labelpad=1)
    else:heat.text(.5,.5,'No SNPs pass filters',transform=heat.transAxes,ha='center')
    heat.set_ylim(n-.5,-.5);heat.set_yticks(range(n),labels,fontsize=7);heat.set_xticks([]);heat.tick_params(length=0)
    locus.set_xlim(1.06,1.24);locus.set_ylim(-.1,1.15);locus.axhline(0,color='#888',lw=.45)
    for mod in base.MODULES:
        for i,(lo,hi) in enumerate(base.INTERVALS[mod]):
            color='#444' if mod==module else '#bfc2c5'
            locus.add_patch(Rectangle((lo/1e6,.05),(hi-lo)/1e6,.24,color=color,lw=0))
            locus.text((lo+hi)/2e6,.51+(i*.48 if mod=='Dennis' else 0),{'Band':'B','Dennis':'D','Ray':'R'}[mod]+('' if mod=='Band' else str(i+1)),ha='center',fontsize=6,color=color)
    locus.set_xticks([1.06,1.12,1.18,1.24]);locus.set_xticklabels(['1.06','1.12','1.18','1.24'],fontsize=6)
    locus.tick_params(axis='x',length=2,pad=1);locus.set_yticks([])
    fig.text(.23,.025,f'Chromosome 18 (Mb)   ·   {len(t)} SNPs shown   ·   grey: insufficient callability',fontsize=6)
    for ax in (heat,locus):
        for s in ax.spines.values():s.set_visible(False)
    pdf.savefig(fig);plt.close(fig)

def draw_frequency(t,dosages,recall,meta,names,cohort,module,cutoff,interval_only,out,pdf,counts):
    groups=frequency_groups(meta,names,module)
    raw,called=freq.frequency_matrix(t,dosages,recall,groups)
    matrix=raw[[2,0,1,3]]
    labels=[groups[0]['label'],'TNLeft','TNRight',groups[1]['label']]
    render(t,matrix,labels,module,cutoff,interval_only,pdf)
    stem=f'{module}_BA{int(cutoff*100)}_{"intervals" if interval_only else "full"}'
    cols=[f'{r.chrom}:{r.pos}:{r.ref}:{r.alt}' for r in t.itertuples()]
    pd.DataFrame(matrix,index=labels,columns=cols).to_csv(out/(stem+'_frequencies.tsv'),sep='\t',na_rep='NA')
    rows=[]
    for i,g in enumerate(groups):
        for j,r in enumerate(t.itertuples()):
            rows.append(dict(group=g['label'],pos=r.pos,n_called=int(called[i,j]),n_total=g['n'],frequency=raw[i+2,j]))
    pd.DataFrame(rows,columns=['group','pos','n_called','n_total','frequency']).to_csv(out/(stem+'_group_callability.tsv.gz'),sep='\t',index=False)
    (out/(module+'_frequency_group_membership.json')).write_text(json.dumps(groups,indent=2)+'\n')

def draw_individual(t,dosages,recall,meta,names,cohort,module,cutoff,interval_only,out,pdf,counts):
    lookup={n:i for i,n in enumerate(names)}
    groups=[dict(indices=[lookup[r.sample]],n=1,label=base.individual_label(r)) for r in meta.itertuples()]
    raw,_=freq.frequency_matrix(t,dosages,recall,groups)
    labels=['TNLeft','TNRight']+[g['label'] for g in groups]
    # Long diagnostic labels require space; use short specimen labels, full metadata saved separately.
    labels=['TNLeft','TNRight']+meta['label'].tolist()
    render(t,raw,labels,module,cutoff,interval_only,pdf)
