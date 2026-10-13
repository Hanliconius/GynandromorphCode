#!/usr/bin/env python3
"""Full-range comparative-QC optix evidence; no BA or focal-state selection."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch,Rectangle
import Fig3_optix_BA099_both_callable as base
import Fig3_optix_full_range_frequency as full

COLORS={'Band':'#087DBB','Dennis':'#D79500','Ray':'#CC70A7'}
plt.rcParams.update({'pdf.fonttype':42,'svg.fonttype':'none','font.size':9})

def draw(values,out,compact=False,combined=None,subtitle=None):
    values=values.copy()
    mapping=values[base.KEY].drop_duplicates().sort_values(base.KEY).reset_index(drop=True)
    mapping["plot_x"]=np.arange(len(mapping)) if compact else mapping.pos/1e6
    values=values.merge(mapping,on=base.KEY,validate="many_to_one")
    fig=plt.figure(figsize=(11,7.5))
    grid=fig.add_gridspec(6,1,height_ratios=[.35,1,.35,1,.35,1],
        left=.14,right=.96,bottom=.50 if compact else .40,top=.90,hspace=.25)
    for i,module in enumerate(full.MODULES):
        t=values.loc[values.module.eq(module)].copy()
        assoc=fig.add_subplot(grid[2*i]); ax=fig.add_subplot(grid[2*i+1],sharex=assoc)
        sites=t[base.KEY+['balanced_accuracy','plot_x']].drop_duplicates()
        assoc.scatter(sites.plot_x,sites.balanced_accuracy,s=3,color=COLORS[module],alpha=.6)
        assoc.set_ylim(.48,1.02);assoc.set_yticks([.5,1],['.5','1'])
        assoc.text(-.025,.5,{'Band':'B','Dennis':'D','Ray':'R'}[module],
            transform=assoc.transAxes,ha='right',va='center',color=COLORS[module],weight='bold')
        assoc.tick_params(axis='x',bottom=False,labelbottom=False)
        for side,marker in [('TN_Left','o'),('TN_Right','v')]:
            s=t.loc[t.row.eq(side)].copy()
            # Each trait is independently polarized; labels express local marker dosage.
            s['dose']=1-s.frequency if module=='Band' else s.frequency
            good=s.dose.notna()
            ax.scatter(s.loc[good,'plot_x'],s.loc[good,'dose'],marker=marker,
                s=12 if marker=='o' else 19,facecolors=COLORS[module] if marker=='o' else 'white',
                edgecolors=COLORS[module],linewidths=.7,alpha=.75,zorder=3)
            ax.scatter(s.loc[~good,'plot_x'],np.full((~good).sum(),-.18),marker='|',s=10,color='#bcbcbc')
        for y in [0,.5,1]:ax.axhline(y,color='#e7e7e7',lw=.6,zorder=0)
        ax.set_ylim(-.25,1.15); ax.set_yticks([0,.5,1],['B-like / B-like','heterozygous','DR-like / DR-like'])
        ax.tick_params(axis='y',length=0,labelsize=8)
        ax.set_xlim((-.5,len(mapping)-.5) if compact else (full.START/1e6,full.END/1e6))
        ax.tick_params(axis='x',labelbottom=i==2,bottom=i==2)
        if i==2:
            ax.set_xlabel('Consecutive SNP columns (shared across B / D / R)' if compact else 'Position on Hmr18 (Mb)')
            if compact:ax.set_xticks([])
        for a in [assoc,ax]:
            for sp in a.spines.values():sp.set_visible(False)
    if compact:
        locus=fig.add_axes([.14,.34,.82,.065])
        locus.set_xlim(full.START/1e6,full.END/1e6);locus.set_ylim(-.2,1.2)
        locus.axhline(0,color='#8b8b8b',lw=1)
        windows=[('B',1111500,1117000,COLORS['Band']),('R1',1156000,1161000,COLORS['Ray']),('R2',1174000,1189000,COLORS['Ray']),('D1',1215000,1217000,COLORS['Dennis']),('D2',1218000,1223000,COLORS['Dennis'])]
        for label,start,end,color in windows:
            locus.add_patch(Rectangle((start/1e6,.1),(end-start)/1e6,.25,color=color,lw=0))
            locus.text((start+end)/2e6,.55 if label!='D2' else .95,label,color=color,ha='center',fontsize=8)
        locus.set_yticks([]);locus.set_xlabel('Hmr18 locus cartoon (Mb); annotated analysis windows',fontsize=8)
        for spine in locus.spines.values():spine.set_visible(False)
        # Every SNP column maps to its position; faint connectors share one gap.
        for row in mapping.itertuples():
            fig.add_artist(ConnectionPatch(xyA=(row.plot_x,-.25),coordsA=ax.transData,
                xyB=(row.pos/1e6,1.15),coordsB=locus.transData,color='#939da5',
                lw=.25,alpha=.08,clip_on=False,zorder=0))
        mapping.to_csv(out/'compact_column_mapping.tsv',sep='\t',index=False)
    cartoon=fig.add_axes([.14,.10 if compact else .14,.82,.15 if compact else .17]);cartoon.axis('off')
    cartoon.text(0,.93,'Expected genotype pairs',weight='bold',transform=cartoon.transAxes)
    for x,title,left,color in [(0,'Band',0,COLORS['Band']),(.36,'Dennis',1,COLORS['Dennis']),(.70,'Ray',1,COLORS['Ray'])]:
        cartoon.text(x,.65,title,color=color,weight='bold',transform=cartoon.transAxes)
        for yy,side,dose in [(.38,'Left',left),(.12,'Right',.5)]:
            cartoon.text(x,yy,side,va='center',transform=cartoon.transAxes,fontsize=8)
            for j in range(2):
                c='#D7EDF5' if (dose==0 or (dose==.5 and j==0)) else '#662A83'
                cartoon.plot([x+.065+j*.038,x+.095+j*.038],[yy,yy],color=c,lw=7,
                    solid_capstyle='round',transform=cartoon.transAxes)
            label='B/B' if dose==0 else ('DR/DR' if dose==1 else 'B/DR')
            cartoon.text(x+.155,yy,label,va='center',transform=cartoon.transAxes,fontsize=8)
    fig.suptitle('optix: Left and Right genotype evidence',x=.14,ha='left',fontsize=14)
    if subtitle:
        fig.text(.14,.92,subtitle,fontsize=9)
    fig.legend(handles=[Line2D([],[],marker='o',color='#555',ls='',label='Left'),
        Line2D([],[],marker='v',markerfacecolor='white',color='#555',ls='',label='Right')],
        loc='upper right',bbox_to_anchor=(.965,.96),frameon=False,ncol=2)
    fig.text(.14,.045,'Small upper tracks: comparative balanced accuracy. Grey ticks: unavailable focal calls.\n'
        'Individual SNPs passing the stated comparative filters; conflicting genotype pairs retained.\n'
        'B-like / DR-like describes independently oriented markers. Cartoon shows expectations; phase is not inferred.',fontsize=8)
    stem='optix_simple_evidence_condensed' if compact else 'optix_simple_evidence'
    for suffix in ['pdf','svg','png']:fig.savefig(out/f'{stem}.{suffix}',dpi=200)
    if combined is not None:combined.savefig(fig)
    plt.close(fig)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--vcf',type=Path,default=base.PROJECT/'external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf')
    a=p.parse_args()
    # Retain existing cohort/scheme eligibility, remove only the BA threshold.
    base.MIN_BA=-np.inf
    sites=base.load_ba_sites()
    sites=sites.loc[sites.chrom.eq(full.CHROM)&sites.pos.between(full.START,full.END)&sites.ref.str.len().eq(1)&sites.alt.str.len().eq(1)].copy()
    groups=full.cohorts()
    samples=sorted(set().union(*(set(g['sample']) for g in groups.values())))
    calls,seen=full.read_calls(a.vcf,sites,samples)
    values,qc=full.build_values(sites,base.load_recalled_pairs(),groups,calls)
    if values.empty:raise RuntimeError('No comparative-QC SNPs')
    a.output_dir.mkdir(parents=True,exist_ok=True)
    values.to_csv(a.output_dir/'plot_values.tsv.gz',sep='\t',index=False)
    qc.to_csv(a.output_dir/'candidate_sites_and_callability.tsv.gz',sep='\t',index=False)
    with PdfPages(a.output_dir/'optix_layout_comparison.pdf') as combined:
        draw(values,a.output_dir,combined=combined)
        draw(values,a.output_dir,compact=True,combined=combined)
    audit={'BA_cutoff':None,'focal_state_selection':False,'phase_inferred':False,
        'comparative_QC':'DP>=2 GQ>=20; >=80% callable in each phenotype group',
        'focal_QC':'reciprocal-mixture pair posterior>=0.95; no original GT fallback',
        'scope':'All scored biallelic SNPs within Hmr18 1060000-1240000 passing comparative QC; unscored sites cannot be polarized',
        'sources':[str(base.BAND_BA_PATH),str(base.DR_BA_PATH),str(base.LR_PATH),str(a.vcf)],
        'counts':{m:len(values.loc[values.module.eq(m),base.KEY].drop_duplicates()) for m in full.MODULES},
        'absent_from_vcf':len(set(map(tuple,sites[base.KEY].values))-seen)}
    (a.output_dir/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('PDF\t'+str(a.output_dir/'optix_simple_evidence.pdf'))
    print('CONDENSED_PDF\t'+str(a.output_dir/'optix_simple_evidence_condensed.pdf'))
    print('SITES\t'+json.dumps(audit['counts']))
if __name__=='__main__':main()
