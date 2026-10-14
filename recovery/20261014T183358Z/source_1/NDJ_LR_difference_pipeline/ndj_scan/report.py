"""Narrative figures made from complete eligible AF counts and window summaries."""
import json
from pathlib import Path
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap, LogNorm
from .common import *

COLORS = ['#D55E00','#0072B2','#009E73','#CC79A7']


def page(title, subtitle):
    fig = plt.figure(figsize=(12,8))
    fig.suptitle(title, x=.06, y=.96, ha='left', fontsize=20, weight='bold')
    fig.text(.06,.902,subtitle, fontsize=10, va='top', color='#444444')
    fig.text(.06,.025,'NDJ | exploratory paired-library scan | all AFs use ALT / (REF + ALT) reads',fontsize=8,color='#666666')
    return fig


def save(pdf,fig):
    pdf.savefig(fig)
    plt.close(fig)


def af_axis(ax, hist, title, vmax=None):
    centres,sx,sy = hex_geometry()
    offsets = np.array([[.5,-.5],[.5,.5],[0,1],[-.5,.5],[-.5,-.5],[0,-1]])*[sx,sy/3]
    # Fixed reversed batlowW table; pure white remains reserved for empty bins.
    palette_path = Path(__file__).with_name('batlowW_r.json')
    colors = np.asarray(json.loads(palette_path.read_text()), dtype=float) / 255
    cmap = ListedColormap(colors[10:])
    cmap.set_bad('white')
    coll = PolyCollection(centres[:,None,:]+offsets[None,:,:],edgecolors='none',cmap=cmap,
                          norm=LogNorm(vmin=1,vmax=max(2,int(vmax if vmax is not None else hist.max()))), rasterized=True)
    coll.set_array(np.ma.masked_where(hist==0,hist))
    ax.add_collection(coll)
    ax.set(xlim=(-.6*sx,1+.6*sx),ylim=(-.4*sy,1+.4*sy),
           xlabel='Left ALT-read fraction',ylabel='Right ALT-read fraction',title=title)
    ax.set_xticks(np.linspace(0,1,6)); ax.set_yticks(np.linspace(0,1,6))
    ax.set_aspect('equal')
    return coll


def genomic_tracks(fig, win, cfg):
    offsets, ticks, labels = {},[],[]
    end = 0
    for c in cfg['chromosomes']:
        offsets[c['contig']] = end
        ticks.append((end+c['length']/2)/1e6)
        labels.append(c['label'])
        end += c['length'] + 300000
    axes = fig.subplots(3,1,sharex=True,gridspec_kw=dict(hspace=.13))
    fig.subplots_adjust(left=.09,right=.96,bottom=.13,top=.82)
    for i,c in enumerate(cfg['chromosomes']):
        t = win.loc[win.chrom==c['contig']]
        x = (t.start+t.end)/2e6 + offsets[c['contig']]/1e6
        good = t.sites >= cfg['min_window_sites']
        color = '#147D92' if i%2==0 else '#7761A7'
        axes[0].plot(x[good],t.loc[good,'mean_abs_delta'],'.',ms=2,color=color)
        axes[1].plot(x[good],t.loc[good,'right_balanced_left_fixed_fraction']-t.loc[good,'left_balanced_right_fixed_fraction'],'.',ms=2,color=color)
        depth = t.usable_depth
        axes[2].plot(x[depth],t.loc[depth,'log2_LR_depth'],'.',ms=2,color=color)
    axes[0].set_ylabel('Mean |Left AF - Right AF|')
    axes[1].set_ylabel('Balanced / fixed\nRight excess')
    axes[2].set_ylabel('Normalized\nlog2(Left / Right depth)')
    for ax in axes[1:]: ax.axhline(0,color='grey',lw=.8)
    axes[2].set_xticks(ticks,labels,fontsize=8)
    axes[2].set_xlabel('Chromosome (physical coordinates within each chromosome)')
    for ax in axes: ax.grid(axis='y',alpha=.2)


def report(run):
    run,cfg=load_config(run)
    bg=json.loads((run/'background.json').read_text())
    win=pd.read_csv(run/'windows.tsv',sep='\t')
    pieces=[]
    audits=[]
    hist=np.zeros(len(hex_geometry()[0]),dtype=np.int64)
    for i,c in enumerate(cfg['chromosomes']):
        out=chrom_dir(run,i)
        if not (out/'qc.complete.json').exists(): raise RuntimeError(f'Missing QC task {i}')
        pieces.append(pd.read_csv(out/'candidate_read_QC.tsv',sep='\t'))
        audits.append(json.loads((out/'scan.complete.json').read_text()))
        hist+=np.load(out/'AF_hex_counts.npy')
    nonempty = [piece for piece in pieces if not piece.empty]
    candidates = (
        pd.concat(nonempty, ignore_index=True)
        if nonempty
        else pd.DataFrame(columns=pieces[0].columns)
    )
    candidates.to_csv(run/'candidate_sites.with_read_QC.tsv',sep='\t',index=False)
    figures=run/'figures'; figures.mkdir(exist_ok=True)
    # Export exact hex occupancy for reproducibility and cross-project comparisons.
    centres=hex_geometry()[0]
    pd.DataFrame(dict(left_AF=centres[:,0],right_AF=centres[:,1],sites=hist)).to_csv(run/'AF_hex_counts.tsv',sep='\t',index=False)
    delta_hist=np.zeros(100,dtype=np.int64)
    edges=np.linspace(-1,1,101)
    total_sites=0
    for i,c in enumerate(cfg['chromosomes']):
        with np.load(chrom_dir(run,i)/'counts.npz') as d:
            _,_,l,r=fractions(d)
            delta_hist+=np.histogram(l-r,bins=edges)[0]
            total_sites+=len(l)
    if int(hist.sum()) != total_sites: raise RuntimeError('AF count conservation failed')
    output=figures/'NDJ_LR_difference_narrative.pdf'
    with PdfPages(output) as pdf:
        fig=page('1  What is actually in the two libraries?',
            f'{total_sites:,} eligible control-panel SNPs; at least {cfg["min_depth"]} allele-assigned reads on BOTH sides.\nNo genotype, posterior, anchor or apparent-difference selection in these plots.')
        ax=fig.add_axes([.08,.18,.43,.61]); h=af_axis(ax,hist,'Every eligible site')
        fig.colorbar(h,ax=ax,shrink=.65,pad=.02,label='Sites per hexagon (log scale)')
        ax=fig.add_axes([.64,.2,.31,.57]); ax.stairs(delta_hist,edges,color='#0072B2',fill=True,alpha=.7)
        ax.set_yscale('symlog',linthresh=10); ax.set(xlabel='Left AF - Right AF',ylabel='Sites (symlog; linear below 10)',title='All-site differences')
        ax.axvline(0,color='black',lw=.8)
        save(pdf,fig)

        fig=page('2  Do differences form genomic regions?',
            f'{cfg["window"]:,}-bp nonoverlapping windows; AF points need {cfg["min_window_sites"]} eligible SNPs.\nDepth uses all BAM positions, including zero coverage. Positive depth values mean relatively more Left coverage.')
        genomic_tracks(fig,win,cfg); save(pdf,fig)

        fig=page('3  How unusual is a site against the observed background?',
            'Background estimation leaves the scored chromosome out. It can include real biological differences.\nThese are descriptive variance scores and ranks, not probabilities of mutation or founder-model support.')
        ax=fig.add_axes([.08,.2,.4,.58])
        rho=[bg['fits'][str(i)]['rho'] for i in range(len(cfg['chromosomes']))]
        ax.bar(range(len(rho)),rho,color='#147D92')
        ax.set(xticks=range(len(rho)),xticklabels=[c['label'] for c in cfg['chromosomes']],ylabel='Background overdispersion rho',xlabel='Chromosome scored',title='Beta-binomial variance sensitivity')
        ax.tick_params(axis='x',labelsize=8)
        fig.text(.56,.78,'What the score means',weight='bold',fontsize=14)
        text=(r'$\Delta_i=f_{L,i}-f_{R,i}$'+'\n\n'+r'$V_i=p_i(1-p_i)[1/n_L+1/n_R$'+'\n'+r'$\quad+\rho(2-1/n_L-1/n_R)]$'+'\n\n'+r'$z_i=\Delta_i/\sqrt{V_i}$'+'\n\n'
              'The empirical tail rank compares |z| with\nother chromosomes at similar pooled AF\nand minimum read depth.\n\nFisher scores are supplied for shortlisted\nsites as read-sampling diagnostics only.\n\nNo multiple-testing or biological\nsignificance claim is made.')
        fig.text(.56,.72,text,fontsize=12,va='top',linespacing=1.3)
        save(pdf,fig)

        # Choose up to three nonoverlapping top AF windows; all points in them are highlighted.
        top=win.loc[win.sites>=cfg['min_window_sites']].nlargest(3,'mean_abs_delta')
        region_hist=[]
        for row in top.itertuples():
            idx=next(i for i,c in enumerate(cfg['chromosomes']) if c['contig']==row.chrom)
            with np.load(chrom_dir(run,idx)/'counts.npz') as d:
                use=(d['pos']>row.start)&(d['pos']<=row.end)
                sub={k:d[k][use] for k in COUNT_KEYS}
                _,_,l,r=fractions(sub)
                region_hist.append(np.bincount(hex_ids(l,r),minlength=len(hist)))
        fig=page('4  Where do candidate regions sit in AF space?',
            'Background = every eligible site. Outlines = occupied bins from the three strongest AF-difference windows.\nThe shortlist changes the overlay, never the underlying full-data distribution.')
        ax=fig.add_axes([.07,.15,.54,.65]); h=af_axis(ax,hist,'All sites + genomic window overlays')
        fig.colorbar(h,ax=ax,shrink=.6,pad=.02,label='All sites per hexagon')
        yy=.76
        for j,(row,hh) in enumerate(zip(top.itertuples(),region_hist)):
            use=hh>0
            ax.scatter(centres[use,0],centres[use,1],s=32,marker='h',facecolors='none',edgecolors=COLORS[j],linewidths=.6)
            label=f'{row.chrom}\n{int(row.start)+1:,}-{int(row.end):,}\n{int(row.sites):,} sites; mean |delta| {row.mean_abs_delta:.3f}'
            fig.text(.65,yy,label,color=COLORS[j],fontsize=10,va='top'); yy-=.16
        fig.text(.65,.22,'Window outlines can overlap.\nCoordinates in tables remain exact.\nNo posterior anchor selection.',fontsize=10,color='#555555')
        pdf.savefig(fig)
        fig.savefig(figures/'NDJ_AF_candidate_regions.pdf')
        fig.savefig(figures/'NDJ_AF_candidate_regions.png',dpi=180)
        plt.close(fig)

        fig=page('5  Inspect the leading single-site candidates',
            'Ranked by background variance score, with a minimum AF difference. QC flags do not establish causality.\nRead checks collapse overlapping mates; counts may differ from caller AD because BAQ/filtering differ.')
        ax=fig.add_axes([.05,.15,.90,.64]);ax.axis('off')
        show=candidates.sort_values('abs_background_z',ascending=False).head(10)
        rows=[]
        for row in show.to_dict('records'):
            rows.append([row['candidate_id'],f"{row['chrom']}\n{int(row['pos']):,}",f"{row['left_AF']:.3f} / {row['right_AF']:.3f}",
                         f"{int(row['left_total'])} / {int(row['right_total'])}",
                         '\n'.join(textwrap.wrap(str(row['QC_flags']),48))[:220]])
        if rows:
            t=ax.table(cellText=rows,colLabels=['ID','Site (1-based)','AF L / R','AD depth L / R','Read-level flags'],colWidths=[.06,.22,.16,.16,.40],loc='upper center',cellLoc='left')
            t.auto_set_font_size(False);t.set_fontsize(8)
            for (r,c),cell in t.get_celld().items():
                cell.set_height(.065 if r==0 else .077)
                if r==0: cell.set_facecolor('#E7F0F3');cell.set_text_props(weight='bold')
        else: ax.text(.1,.7,'No sites met the configured candidate effect-size threshold.',fontsize=14)
        fig.text(.06,.09,'Complete counts, strand balance, mapping/base quality, distinct starts, clipping and nearby-indel metrics are in candidate_sites.with_read_QC.tsv.',fontsize=9)
        save(pdf,fig)

        fig=page('6  Check technical patterns before interpreting depth',
            'GC and coverage relationships help flag technical explanations. These plots do not automatically correct them.\nMedian depth normalization removes global library size differences and also hides whole-genome dosage shifts.')
        axes=fig.subplots(1,2);fig.subplots_adjust(left=.08,right=.95,bottom=.18,top=.8,wspace=.3)
        valid=win.usable_depth & win.log2_LR_depth.notna()
        axes[0].scatter(win.loc[valid,'gc'],win.loc[valid,'log2_LR_depth'],s=8,alpha=.35,color='#147D92')
        axes[0].set(xlabel='Window GC fraction',ylabel='Normalized log2(Left / Right depth)',title='Depth versus sequence composition')
        axes[1].scatter(win.loc[valid,'mean_abs_delta'],win.loc[valid,'log2_LR_depth'],s=8,alpha=.35,color='#7761A7')
        axes[1].set(xlabel='Window mean absolute AF difference',ylabel='Normalized log2(Left / Right depth)',title='AF shifts and depth shifts')
        for ax in axes: ax.axhline(0,color='grey',lw=.8)
        save(pdf,fig)

        fig=page('7  What this run can and cannot establish',
                 'The aim is to localize unexpected Left/Right differences without assuming the known biological answer.')
        lines=[
          f'Input: {cfg["left"]} and {cfg["right"]}; one library per side.',
          f'All-site AF plot: {total_sites:,} SNPs, with technical filters audited per chromosome.',
          'Site scan: effect size + depth-aware empirical ranking; read-level inspection of the shortlist.',
          'Regional scan: absolute AF differences, balanced/fixed contrasts, and independent BAM depth.',
          'Signed ALT differences are unphased: REF/ALT orientation can cancel a regional signal.',
          'Control-panel ascertainment: NDJ-private variants absent from the BCF are not tested.',
          'Z is reported separately from the default autosomal background; no sex/ploidy assignment is made.',
          'No replicate validation, absolute-ploidy inference, or founder-model selection is claimed.',
          'Repeat annotation is optional; without an indexed BED it is explicitly marked NOT_ASSESSED.',
          'Next: inspect candidate regions and read evidence before choosing a biological model.',
          'Per-chromosome raw count NPZ/TSV outputs are reusable by subsequent founder-model analyses.'
        ]
        y=.80
        for line in lines:
            fig.text(.08,y,'• '+line,fontsize=11,va='top');y-=.061
        save(pdf,fig)
    # Full physical-coordinate tracks for every chromosome, including windows with no SNPs.
    with PdfPages(figures/'NDJ_chromosome_tracks.pdf') as pdf:
        for c in cfg['chromosomes']:
            t=win.loc[win.chrom==c['contig']]
            fig=page(f'Chromosome {c["label"]}  |  {c["contig"]}',
                'Every depth window is retained; AF summaries use eligible control-panel SNPs. Blank AF windows are not zero differences.')
            axes=fig.subplots(3,1,sharex=True);fig.subplots_adjust(top=.8,bottom=.13,left=.1,right=.95,hspace=.2)
            x=(t.start+t.end)/2e6
            good=t.sites>=cfg['min_window_sites']
            axes[0].plot(x[good],t.loc[good,'mean_abs_delta'],color='#147D92');axes[0].set_ylabel('Mean |delta AF|')
            axes[1].plot(x[good],t.loc[good,'left_balanced_fraction'],label='Left',color=COLORS[1],alpha=.75)
            axes[1].plot(x[good],t.loc[good,'right_balanced_fraction'],label='Right',color=COLORS[0],alpha=.75)
            axes[1].set_ylabel('Balanced AF fraction');axes[1].legend(fontsize=8)
            axes[2].plot(x,t.left_depth_relative,label='Left',color=COLORS[1],alpha=.75)
            axes[2].plot(x,t.right_depth_relative,label='Right',color=COLORS[0],alpha=.75)
            axes[2].set(xlabel='Genomic position (Mb)',ylabel='Relative mean depth');axes[2].legend(fontsize=8)
            save(pdf,fig)
    write_json(run/'COMPLETE.json',dict(version=VERSION,eligible_sites=total_sites,
        candidate_sites=len(candidates), per_chromosome=audits,
        narrative_pdf=str(output), status='exploratory_scan_complete_not_biological_validation'))
    print(f'COMPLETE\tsites={total_sites:,}\nNARRATIVE_PDF\t{output}',flush=True)
