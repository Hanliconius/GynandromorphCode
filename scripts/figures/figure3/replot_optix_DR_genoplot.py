#!/usr/bin/env python3
"""Dennis/Ray tiled genoplots from cached independently selected values."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle,ConnectionPatch
from matplotlib.backends.backend_pdf import PdfPages
from replot_optix_pruned import select,INTERVALS

KEY=['chrom','pos','ref','alt']
COLORS={'Dennis':'#D79500','Ray':'#CC70A7'}
plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none','font.size':9})

def draw(values,cutoff,scope,out,combined,pk_samples=None,focal_callable=False):
    blocks=[];offset=0;maps=[];arrays=[]
    pk_samples=pk_samples or []
    rows=['present']+pk_samples+['TN_Left','TN_Right','absent']
    nrows=len(rows)
    for module in ['Ray','Dennis']:
        t=values.loc[values.module.eq(module)].copy()
        sites=t[KEY].drop_duplicates().sort_values(KEY).reset_index(drop=True)
        if sites.empty:continue
        sites['column']=np.arange(offset,offset+len(sites));sites['module']=module
        t=t.merge(sites.drop(columns='module'),on=KEY,validate='many_to_one')
        if t.duplicated(['row','column']).any():raise RuntimeError('Duplicate genoplot cell')
        arrays.append(t.pivot(index='row',columns='column',values='frequency').reindex(index=rows,columns=sites.column).to_numpy())
        blocks.append((module,offset,len(sites)));maps.append(sites);offset+=len(sites)
    fig=plt.figure(figsize=(11,5.5+len(pk_samples)*.18))
    fig.suptitle('optix: Dennis and Ray genoplot',x=.19,ha='left',fontsize=14)
    fig.text(.19,.90,f'BA >= {cutoff:.2f} | '+('Trait-specific labelled intervals' if scope else 'Full region')+(' | PK individual calls' if pk_samples else '')+(' | focal no-calls removed' if focal_callable else ''),fontsize=10)
    if not arrays:
        fig.text(.5,.5,'No SNPs pass these comparative filters',ha='center')
    else:
        mapping=pd.concat(maps,ignore_index=True);matrix=np.concatenate(arrays,axis=1)
        locus=fig.add_axes([.19,.69,.72,.12]);heat=fig.add_axes([.19,.30,.72,.26])
        locus.set_xlim(1.06,1.24);locus.set_ylim(-.2,1.4)
        locus.hlines(0,1.06,1.24,color='#777',lw=1)
        for module,labels in [('Ray',['R1','R2']),('Dennis',['D1','D2'])]:
            for label,(start,end) in zip(labels,INTERVALS[module]):
                locus.add_patch(Rectangle((start/1e6,.1),(end-start)/1e6,.25,color=COLORS[module],lw=0))
                locus.text((start+end)/2e6,.63 if label!='D2' else 1.05,label,color=COLORS[module],ha='center',fontsize=9)
        locus.set_yticks([]);locus.set_xlabel('Hmr18 locus (Mb)',fontsize=8)
        cmap=LinearSegmentedColormap.from_list('dosage',[(0,'#D7EDF5'),(.5,'#ECCB63'),(1,'#662A83')]);cmap.set_bad('#C8C8C8')
        image=heat.imshow(matrix,cmap=cmap,vmin=0,vmax=1,aspect='auto',interpolation='nearest',extent=(-.5,offset-.5,nrows-.5,-.5))
        heat.set_yticks(range(nrows),['Pattern present']+pk_samples+['Left adjusted','Right adjusted','Pattern absent'])
        heat.set_xticks([]);heat.tick_params(axis='y',length=0)
        for module,start,count in blocks:
            if start:heat.axvline(start-.5,color='white',lw=3)
            heat.text(start+(count-1)/2,1.06,f'{module} (n={count})',transform=heat.get_xaxis_transform(),ha='center',va='bottom',color=COLORS[module],weight='bold')
        for row in mapping.itertuples():
            fig.add_artist(ConnectionPatch(xyA=(row.pos/1e6,-.18),coordsA=locus.transData,xyB=(row.column,-.5),coordsB=heat.transData,color='#929aa3',lw=.3,alpha=.15,zorder=0,clip_on=False))
        for ax in [locus,heat]:
            for spine in ax.spines.values():spine.set_visible(False)
        bar=fig.colorbar(image,cax=fig.add_axes([.935,.30,.013,.26]));bar.set_ticks([0,.5,1]);bar.set_label('Associated allele frequency / dosage',fontsize=8)
        mapping.to_csv(out/f'BA{int(cutoff*100)}_{"intervals" if scope else "full"}_columns.tsv',sep='\t',index=False)
        # Paired glyphs are schematic local genotype expectations, not phased haplotypes.
        cartoon=fig.add_axes([.19,.12,.72,.12]);cartoon.axis('off')
        cartoon.text(0,.95,'Expected at Dennis / Ray-associated markers',weight='bold',fontsize=9)
        for x,side,cols,label in [(0,'Left',['#662A83','#662A83'],'DR/DR'),(.53,'Right',['#D7EDF5','#662A83'],'B/DR')]:
            cartoon.text(x,.38,side,va='center')
            for j,c in enumerate(cols):cartoon.plot([x+.075+j*.05,x+.11+j*.05],[.38,.38],color=c,lw=9,solid_capstyle='round')
            cartoon.text(x+.19,.38,label,va='center')
        cartoon.set_xlim(0,1);cartoon.set_ylim(0,1)
    fig.text(.19,.035,'Comparative rows: measured frequencies. Focal calls: reciprocal-mixture posterior >=0.95.\n'+('Unavailable focal sites removed. ' if focal_callable else 'Grey = unavailable. ')+('PK: DP>=2 GQ>=20; grey = unavailable PK call. ' if pk_samples else '')+'Conflicting calls retained; phase is not inferred.',fontsize=8)
    stem=f'optix_DR_genoplot_BA{int(cutoff*100)}_{"labelled_intervals" if scope else "full_region"}'+('_with_PK' if pk_samples else '')
    fig.savefig(out/(stem+'.pdf'));combined.savefig(fig);plt.close(fig)

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,required=True);a=p.parse_args()
    v=pd.read_csv(a.input_dir/'plot_values.tsv.gz',sep='\t');v=v.loc[v.module.isin(['Dennis','Ray'])]
    out=a.input_dir/'DR_genoplot';out.mkdir(exist_ok=True)
    records=[]
    with PdfPages(out/'optix_Dennis_Ray_genoplot_comparison.pdf') as combined:
        for cutoff in [.8,.9,.95]:
            for scope in [False,True]:
                subset=select(v,cutoff,scope)
                draw(subset,cutoff,scope,out,combined)
                records.append({'BA_min':cutoff,'labelled_intervals_only':scope,'counts':{m:len(subset.loc[subset.module.eq(m),KEY].drop_duplicates()) for m in ['Dennis','Ray']}})
    (out/'audit.json').write_text(json.dumps({'input':str(a.input_dir/'plot_values.tsv.gz'),'focal_state_selection':False,'cluster_summary':False,'comparisons':records},indent=2)+'\n')
    print('PDF\t'+str(out/'optix_Dennis_Ray_genoplot_comparison.pdf'))
if __name__=='__main__':main()
