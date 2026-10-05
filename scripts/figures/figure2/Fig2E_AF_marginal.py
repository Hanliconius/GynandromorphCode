#!/usr/bin/env python3
"""Figure 2E — one-dimensional allele-frequency distributions by original genotype."""

from __future__ import annotations
import argparse, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
from common import GT_COLOURS, add_panel_letter, loess_like, panel_output, project_root, set_publication_style

DEFAULT_CACHE=Path("analysis/figures/pairwise_AF_cli/AF_narrative_frames/original_GT_histograms/ERR2298208_vs_TNLeft_original_GT_counts.npz")

def main():
    p=argparse.ArgumentParser(); p.add_argument("--project-root",default=None); p.add_argument("--cache",default=str(DEFAULT_CACHE)); p.add_argument("--output",default=None); a=p.parse_args()
    root=project_root(a.project_root); cache_path=Path(a.cache)
    if not cache_path.is_absolute(): cache_path=root/cache_path
    if not cache_path.is_file(): raise FileNotFoundError(cache_path)
    output=panel_output(root,"E",a.output); set_publication_style(7.6)
    with np.load(cache_path,allow_pickle=False) as cache:
        hist=(cache["hist"] if "hist" in cache else cache["counts"]).copy()
        edges=cache["edges"].copy() if "edges" in cache else np.linspace(0.0,1.0,hist.shape[-1]+1)
    if hist.ndim!=3 or hist.shape[0]<2 or hist.shape[1]<3: raise RuntimeError(f"Unexpected AF histogram shape: {hist.shape}")
    control=hist[0,:3,:].astype(float); focal=hist[1,:3,:].astype(float); control/=control.sum(); focal/=focal.sum()
    centres=0.5*(edges[:-1]+edges[1:]); bw=float(edges[1]-edges[0]); bar_width=bw*0.27; offsets={0:-bar_width,2:0.0,1:+bar_width}; order=[0,2,1]
    fig=plt.figure(figsize=(5.25,3.10)); add_panel_letter(fig,"E",x=0.01,y=0.99,size=14); ax=fig.add_axes([0.12,0.20,0.84,0.68])
    for gt in order:
        ax.bar(centres+offsets[gt],focal[gt],width=bar_width*0.95,color=GT_COLOURS[gt],edgecolor="none",alpha=0.68,zorder=3 if gt==1 else 2)
    for gt in order:
        x=centres[1:-1]; smooth=loess_like(x,control[gt,1:-1],span=0.26)
        ax.plot(np.r_[0.0,x,1.0],np.r_[control[gt,0],smooth,control[gt,-1]],color=GT_COLOURS[gt],linestyle=(0,(3.5,2.5)),linewidth=1.25,alpha=0.92,zorder=5)
    ax.set_xlim(0,1); ymax=max(float(np.nanmax(focal)),float(np.nanmax(control))); ax.set_ylim(0,max(0.018,min(0.032,ymax*1.04)))
    ax.yaxis.set_major_formatter(PercentFormatter(1.0,decimals=1)); ax.set_xlabel("Allelic Frequency (ALT allele)"); ax.set_ylabel("SNP fraction"); ax.set_xticks([0,.25,.5,.75,1])
    ax.axvline(.5,color="#888888",linewidth=.55,linestyle=(0,(2,2)),zorder=1); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.annotate("excess of\nintermediate frequencies",xy=(.18,ax.get_ylim()[1]*.37),xytext=(.49,ax.get_ylim()[1]*.68),ha="center",va="center",fontsize=7.3,color="#777777",arrowprops=dict(arrowstyle="->",color="#777777",lw=.9))
    ax.annotate("",xy=(.84,ax.get_ylim()[1]*.38),xytext=(.51,ax.get_ylim()[1]*.68),arrowprops=dict(arrowstyle="->",color="#777777",lw=.9))
    handles=[Patch(facecolor=GT_COLOURS[0],edgecolor="none",label="hom. REF"),Patch(facecolor=GT_COLOURS[1],edgecolor="none",label="heterozygous"),Patch(facecolor=GT_COLOURS[2],edgecolor="none",label="hom. ALT"),Patch(facecolor="#777777",edgecolor="none",alpha=.55,label="gynandromorph – left leg"),Line2D([0],[0],color="#444444",linestyle=(0,(3.5,2.5)),linewidth=1.25,label="diploid control")]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.56,.99),ncol=3,frameon=False,fontsize=6.7,handlelength=1.3,columnspacing=1.0)
    fig.savefig(output,bbox_inches="tight"); plt.close(fig)
    print(f"INPUT\t{cache_path}"); print(f"FOCAL_SITES\t{int(hist[1].sum())}"); print(f"CONTROL_SITES\t{int(hist[0].sum())}"); print(f"WROTE\t{output}")

if __name__=="__main__": main()
