#!/usr/bin/env python3
"""Figure 2D — final Smudgeplot focal arrangements."""

from __future__ import annotations
import argparse, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
from common import add_panel_letter, draw_smudge_axis, panel_output, project_root, read_smu, set_publication_style, shared_smudge_scale, smudge_density

YMAX=120.0
INPUTS={
 "diploid":Path("analysis/smudgeplot_standard_rerun_20261002/ERR1051012_diploid/L6/kmerpairs.smu"),
 "pseudo":Path("analysis/smudgeplot_lowerL_20261002/pseudoTetraploid_ERR1307111_ERR1307112/L6/kmerpairs.smu"),
 "left":Path("analysis/smudgeplot_arrangements/TN_Left/kmerpairs.smu"),
 "right":Path("analysis/smudgeplot_arrangements/TN_Right/kmerpairs.smu"),
 "t2":Path("analysis/smudgeplot_lowerL_20261002/TN_T2_merge/L6/kmerpairs.smu"),
 "t3":Path("analysis/smudgeplot_lowerL_20261002/TN_T3_merge/L6/kmerpairs.smu"),
}

def prepare(root):
    out={}; densities=[]
    for key,rel in INPUTS.items():
        path=root/rel
        data=read_smu(path)
        density,xe,ye=smudge_density(data,ymax=YMAX)
        out[key]=(path,data,density,xe,ye); densities.append(density)
    return out,*shared_smudge_scale(densities)

def main():
    p=argparse.ArgumentParser(); p.add_argument("--project-root",default=None); p.add_argument("--output",default=None); a=p.parse_args()
    root=project_root(a.project_root); output=panel_output(root,"D",a.output); set_publication_style(7.5)
    prepared,vmin,vmax,levels=prepare(root)
    fig,axes=plt.subplots(2,2,figsize=(5.15,5.25),sharex=True,sharey=True); add_panel_letter(fig,"D",x=0.01,y=0.995,size=14)
    specs=[("left","left legs combined\nT2 + T3; L=12"),("right","right legs combined\nT2 + T3; L=12"),("t2","T2 legs combined\nleft + right; L=6"),("t3","T3 legs combined\nleft + right; L=6")]
    artist=None
    for ax,(key,title) in zip(axes.ravel(),specs):
        path,data,density,xe,ye=prepared[key]
        artist=draw_smudge_axis(ax,density,xe,ye,vmin=vmin,vmax=vmax,levels=levels,ymax=YMAX,title=title)
        print("PANEL",key,path,data["rows"],data["total_pairs"],sep="\t")
    for ax in axes[-1,:]: ax.set_xlabel("k-mer minor fraction")
    for ax in axes[:,0]: ax.set_ylabel("k-mer coverage")
    cb=fig.colorbar(artist,ax=axes.ravel().tolist(),fraction=0.03,pad=0.025); cb.set_label("k-mer pairs",fontsize=7); cb.ax.tick_params(labelsize=6,length=2)
    fig.subplots_adjust(left=0.11,right=0.88,bottom=0.10,top=0.93,hspace=0.27,wspace=0.18)
    fig.savefig(output,bbox_inches="tight"); plt.close(fig); print(f"SHARED_DENSITY_SCALE\t{vmin}\t{vmax}"); print(f"WROTE\t{output}")

if __name__=="__main__": main()
