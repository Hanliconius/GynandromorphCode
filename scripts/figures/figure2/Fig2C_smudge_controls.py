#!/usr/bin/env python3
"""Figure 2C — final Smudgeplot control comparisons."""

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
    root=project_root(a.project_root); output=panel_output(root,"C",a.output); set_publication_style(7.5)
    prepared,vmin,vmax,levels=prepare(root)
    fig,axes=plt.subplots(2,1,figsize=(2.85,5.55),sharex=True,sharey=True); add_panel_letter(fig,"C",x=0.01,y=0.995,size=14)
    specs=[("diploid","single diploid (2n)\nERR1051012; L=6"),("pseudo","2 diploids combined\n(2n + 2n); L=6")]
    artist=None
    for ax,(key,title) in zip(axes,specs):
        path,data,density,xe,ye=prepared[key]
        artist=draw_smudge_axis(ax,density,xe,ye,vmin=vmin,vmax=vmax,levels=levels,ymax=YMAX,title=title)
        print("PANEL",key,path,data["rows"],data["total_pairs"],sep="\t")
    axes[-1].set_xlabel("k-mer minor fraction")
    for ax in axes: ax.set_ylabel("k-mer coverage")
    cb=fig.colorbar(artist,ax=axes.tolist(),fraction=0.045,pad=0.035); cb.set_label("k-mer pairs",fontsize=7); cb.ax.tick_params(labelsize=6,length=2)
    fig.subplots_adjust(left=0.21,right=0.84,bottom=0.09,top=0.94,hspace=0.25)
    fig.savefig(output,bbox_inches="tight"); plt.close(fig); print(f"SHARED_DENSITY_SCALE\t{vmin}\t{vmax}"); print(f"WROTE\t{output}")

if __name__=="__main__": main()
