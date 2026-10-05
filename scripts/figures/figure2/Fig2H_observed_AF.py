#!/usr/bin/env python3
"""Figure 2H — observed Left/Right AF density, model overlay and genotype marginals."""

from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
from common import GT_COLOURS, MODEL_COLOUR, add_panel_letter, add_predictive_region, expected_mixed_af, load_luka_bleute, loess_like, panel_output, project_root, set_publication_style

BASE=Path("analysis/figures/pairwise_AF_cli")
DEFAULT_BINS=BASE/"TN_merged_control_sites_hex40_batlowW.AF_bins.tsv.gz"
DEFAULT_HIST_ROOT=BASE/"AF_narrative_frames/original_GT_histograms"
DEFAULT_COUNTS=DEFAULT_HIST_ROOT/"original_GT_counts.npz"
DEFAULT_AUDIT=DEFAULT_HIST_ROOT/"original_GT_counts.audit.json"
EXPECTED_N=31_126_758

def draw_marginal(axis,counts,vertical,cap=.015):
    counts=np.asarray(counts,dtype=float)
    total=float(counts.sum())
    if counts.shape[0]<3 or total<=0: raise RuntimeError("Invalid marginal histogram")
    frac=counts[:3]/total; x=(np.arange(counts.shape[1])+.5)/counts.shape[1]
    for gt in (0,1,2):
        smooth=loess_like(x[1:-1],frac[gt,1:-1],span=.22); colour=GT_COLOURS[gt]
        if vertical:
            axis.plot(smooth,x[1:-1],color=colour,lw=.85); axis.fill_betweenx(x[1:-1],0,smooth,color=colour,alpha=.10)
            for idx in (0,len(x)-1): axis.plot([0,min(float(frac[gt,idx]),cap)],[x[idx],x[idx]],color=colour,lw=.8)
        else:
            axis.plot(x[1:-1],smooth,color=colour,lw=.85); axis.fill_between(x[1:-1],0,smooth,color=colour,alpha=.10)
            for idx in (0,len(x)-1): axis.plot([x[idx],x[idx]],[0,min(float(frac[gt,idx]),cap)],color=colour,lw=.8)
    if vertical:
        axis.set_ylim(-.02,1.02); axis.set_xlim(cap*1.03,0); axis.set_xticks([0,cap]); axis.xaxis.set_major_formatter(PercentFormatter(1.,decimals=1)); axis.set_yticks([])
        for s in ("left","top","right"): axis.spines[s].set_visible(False)
    else:
        axis.set_xlim(-.02,1.02); axis.set_ylim(0,cap*1.03); axis.set_yticks([0,cap]); axis.yaxis.set_major_formatter(PercentFormatter(1.,decimals=1)); axis.set_xticks([])
        for s in ("top","right","bottom"): axis.spines[s].set_visible(False)

def main():
    p=argparse.ArgumentParser(); p.add_argument("--project-root",default=None); p.add_argument("--bins",default=str(DEFAULT_BINS)); p.add_argument("--counts",default=str(DEFAULT_COUNTS)); p.add_argument("--audit",default=str(DEFAULT_AUDIT)); p.add_argument("--output",default=None); a=p.parse_args()
    root=project_root(a.project_root)
    paths=[Path(a.bins),Path(a.counts),Path(a.audit)]
    paths=[x if x.is_absolute() else root/x for x in paths]
    bins_path,counts_path,audit_path=paths
    for path in paths:
        if not path.is_file(): raise FileNotFoundError(path)
    output=panel_output(root,"H",a.output); set_publication_style(7.2); cmap=load_luka_bleute(root)
    bins=pd.read_csv(bins_path,sep="\t"); required={"left_centre","right_centre","width","height","count"}; missing=required-set(bins.columns)
    if missing: raise RuntimeError(f"Hex-bin table missing columns: {sorted(missing)}")
    with np.load(counts_path,allow_pickle=False) as cache: hist=cache["hist"].copy()
    metadata=json.loads(audit_path.read_text()); samples=metadata["signature"]["samples"]
    li=samples.index("TN_Left_merged"); ri=samples.index("TN_Right_merged"); left_counts=hist[li]; right_counts=hist[ri]
    n=int(bins["count"].sum())
    if n!=int(left_counts.sum()) or n!=int(right_counts.sum()): raise RuntimeError("Observed density/marginal totals disagree")
    if n!=EXPECTED_N: raise RuntimeError(f"Expected {EXPECTED_N:,} retained sites; found {n:,}")
    occupied=bins.loc[bins["count"]>0].copy(); unit=np.array([[0,-.5],[.5,-.25],[.5,.25],[0,.5],[-.5,.25],[-.5,-.25]])
    centres=occupied[["left_centre","right_centre"]].to_numpy(float); sizes=occupied[["width","height"]].to_numpy(float); vertices=centres[:,None,:]+unit[None,:,:]*sizes[:,None,:]
    fig=plt.figure(figsize=(5.25,5.55)); add_panel_letter(fig,"H",x=.01,y=.995,size=14); fig.text(.51,.977,"observed allelic frequencies",ha="center",va="top",fontsize=9.2,fontweight="bold")
    ax=fig.add_axes([.225,.215,.61,.61]); top=fig.add_axes([.225,.835,.61,.075]); side=fig.add_axes([.115,.215,.075,.61]); cax=fig.add_axes([.30,.085,.46,.025])
    collection=PolyCollection(vertices,array=occupied["count"].to_numpy(float),cmap=cmap,norm=LogNorm(vmin=1.,vmax=float(occupied["count"].max())),edgecolors="none",linewidths=0,rasterized=False); ax.add_collection(collection)
    pad=max(float(bins["width"].max()),float(bins["height"].max()))/2+.003; limits=(-pad,1+pad); ax.set_xlim(*limits); ax.set_ylim(*limits); ax.set_aspect("equal",adjustable="box"); ax.plot([0,1],[0,1],color="white",lw=.5,alpha=.40,zorder=4)
    for pl in (0.,.5,1.):
        for pr in (0.,.5,1.):
            ql,qr=expected_mixed_af(pl,pr,c_left=.20,c_right=.20,error=.001); add_predictive_region(ax,ql,qr,depth_left=40,depth_right=40,concentration=22,mass=.50,colour=MODEL_COLOUR,linewidth=1.3,alpha=.75)
            ax.plot(ql,qr,marker="o",linestyle="none",markersize=3.2,markerfacecolor=MODEL_COLOUR,markeredgecolor="white",markeredgewidth=.5,zorder=10)
    ax.set_xticks([0,.25,.5,.75,1]); ax.set_yticks([0,.25,.5,.75,1]); ax.set_xlabel("Left ALT-read fraction",labelpad=2); ax.set_ylabel("Right ALT-read fraction",labelpad=2); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    draw_marginal(top,left_counts,False); draw_marginal(side,right_counts,True)
    cb=fig.colorbar(collection,cax=cax,orientation="horizontal"); cb.set_label("observed sites / hexagon",fontsize=6.8,labelpad=2); cb.ax.tick_params(labelsize=5.8,length=2,pad=1); cb.outline.set_linewidth(.5)
    legend=[Line2D([0],[0],color=GT_COLOURS[0],lw=1.3,label="original 0/0"),Line2D([0],[0],color=GT_COLOURS[1],lw=1.3,label="original 0/1"),Line2D([0],[0],color=GT_COLOURS[2],lw=1.3,label="original 1/1"),Line2D([0],[0],marker="o",color=MODEL_COLOUR,lw=1.,markersize=3.5,markerfacecolor=MODEL_COLOUR,label="20% model")]
    fig.legend(handles=legend,loc="lower center",bbox_to_anchor=(.52,.005),ncol=4,frameon=False,fontsize=5.8,columnspacing=.8,handlelength=1.4)
    fig.text(.83,.925,f"n = {n:,}",ha="right",va="center",fontsize=6.5,color="#555555")
    fig.savefig(output,bbox_inches="tight"); plt.close(fig)
    print(f"BINS\t{bins_path}"); print(f"COUNTS\t{counts_path}"); print(f"SITES\t{n}"); print(f"WROTE\t{output}")

if __name__=="__main__": main()
