#!/usr/bin/env python3
"""Figure 2G — simulated paired allele-frequency distributions for four founder configurations."""

from __future__ import annotations
import argparse, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
from common import MODEL_COLOUR, add_panel_letter, add_predictive_region, expected_mixed_af, load_luka_bleute, load_simulation_pair, panel_output, project_root, set_publication_style, validate_fraction_arrays

NARRATIVE=Path("analysis/figures/pairwise_AF_cli/AF_narrative_frames")
SCENARIOS=[
 {"name":"AB / A","source":NARRATIVE/"AB_A_digametic/AB_A_simulated_AF.npz","prefix":"mix20","states":[(0.,0.),(.5,0.),(.5,1.),(1.,1.)]},
 {"name":"AA / BC","source":NARRATIVE/"gyno_three_founder/gyno_simulated_AF.npz","prefix":"mix20","states":[(l,r) for l in (0.,1.) for r in (0.,.5,1.)]},
 {"name":"AB / AC","source":NARRATIVE/"simulated_founder_empirical_depth/simulated_AF.npz","prefix":"mix20_model3","states":[(l,r) for l in (0.,.5,1.) for r in (0.,.5,1.) if (l,r) not in {(0.,1.),(1.,0.)}]},
 {"name":"AB / CD","source":NARRATIVE/"simulated_founder_empirical_depth/simulated_AF.npz","prefix":"mix20_model4","states":[(l,r) for l in (0.,.5,1.) for r in (0.,.5,1.)]},
]

def main():
    p=argparse.ArgumentParser(); p.add_argument("--project-root",default=None); p.add_argument("--output",default=None); p.add_argument("--gridsize",type=int,default=40); a=p.parse_args()
    root=project_root(a.project_root); output=panel_output(root,"G",a.output); set_publication_style(7.3); cmap=load_luka_bleute(root)
    loaded=[]; sizes=set()
    for s in SCENARIOS:
        path=root/s["source"]; left,right=load_simulation_pair(path,s["prefix"]); validate_fraction_arrays(left,right); loaded.append((s,path,left,right)); sizes.add(len(left))
    if len(sizes)!=1: raise RuntimeError(f"Simulation panel sizes differ: {sorted(sizes)}")
    fig,axes=plt.subplots(2,2,figsize=(5.15,5.25),sharex=True,sharey=True); add_panel_letter(fig,"G",x=.01,y=.995,size=14)
    fig.text(.50,.975,"simulated allelic frequencies",ha="center",va="top",fontsize=9.2,fontweight="bold")
    collections=[]
    for ax,(s,path,left,right) in zip(axes.ravel(),loaded):
        collection=ax.hexbin(left,right,gridsize=a.gridsize,extent=(0,1,0,1),mincnt=1,cmap=cmap,linewidths=0,rasterized=False); collections.append(collection)
        for p_left,p_right in s["states"]:
            ql,qr=expected_mixed_af(p_left,p_right,c_left=.20,c_right=.20,error=.001)
            add_predictive_region(ax,ql,qr,depth_left=40,depth_right=40,concentration=22,mass=.50,colour=MODEL_COLOUR,linewidth=1.25,alpha=.72)
            ax.plot(ql,qr,marker="o",linestyle="none",markersize=3.2,markerfacecolor=MODEL_COLOUR,markeredgecolor="white",markeredgewidth=.5,zorder=10)
        ax.plot([0,1],[0,1],color="white",linewidth=.5,alpha=.35,zorder=4); ax.set_xlim(-.018,1.018); ax.set_ylim(-.018,1.018); ax.set_aspect("equal",adjustable="box")
        ax.set_xticks([0,.5,1]); ax.set_yticks([0,.5,1]); ax.set_title(s["name"],fontsize=7.8,pad=3); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        print(f"SOURCE\t{s['name']}\t{path}\t{len(left)}")
    maximum=max(float(c.get_array().max()) for c in collections); norm=LogNorm(vmin=1.,vmax=max(1.,maximum))
    for c in collections: c.set_norm(norm)
    for ax in axes[-1,:]: ax.set_xlabel("Left ALT-read fraction")
    for ax in axes[:,0]: ax.set_ylabel("Right ALT-read fraction")
    cb=fig.colorbar(collections[0],ax=axes.ravel().tolist(),fraction=.03,pad=.025); cb.set_label("simulated sites / hexagon",fontsize=6.8); cb.ax.tick_params(labelsize=6,length=2)
    fig.text(.50,.025,"20% reciprocal read mixing; magenta points = expected centroids; rings = 50% predictive regions",ha="center",va="bottom",fontsize=6.2,color="#555555")
    fig.subplots_adjust(left=.11,right=.88,bottom=.09,top=.92,hspace=.16,wspace=.13); fig.savefig(output,bbox_inches="tight"); plt.close(fig)
    print(f"SHARED_MAX_SITES_PER_HEX\t{maximum}"); print(f"WROTE\t{output}")

if __name__=="__main__": main()
