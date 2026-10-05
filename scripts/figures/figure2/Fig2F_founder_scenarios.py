#!/usr/bin/env python3
"""Figure 2F — conceptual founder/haplotype configurations with reciprocal mixing."""

from __future__ import annotations
import argparse, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0,str(HERE))
from common import HAPLOTYPE_COLOURS, add_panel_letter, panel_output, project_root, set_publication_style

SCENARIOS=[
 {"model":"AB / A","title":"two haplotypes","subtitle":"diploid / haploid","left":("A","B"),"right":("A",)},
 {"model":"AA / BC","title":"three haplotypes","subtitle":"one duplicated","left":("A","A"),"right":("B","C")},
 {"model":"AB / AC","title":"three haplotypes","subtitle":"one shared","left":("A","B"),"right":("A","C")},
 {"model":"AB / CD","title":"four haplotypes","subtitle":"no shared haplotype","left":("A","B"),"right":("C","D")},
]

def draw_cell(ax,x0,y0,width,height,haps,label):
    ax.add_patch(FancyBboxPatch((x0,y0),width,height,boxstyle="round,pad=0.018,rounding_size=0.035",facecolor="white",edgecolor="#8A8A8A",linewidth=.8))
    xs=[x0+width/2] if len(haps)==1 else [x0+width*.34,x0+width*.66]
    for xpos,hap in zip(xs,haps):
        ax.add_patch(FancyBboxPatch((xpos-width*.055,y0+height*.19),width*.11,height*.60,boxstyle="round,pad=0.008,rounding_size=0.018",facecolor=HAPLOTYPE_COLOURS[hap],edgecolor="none"))
        ax.text(xpos,y0+height*.49,hap,ha="center",va="center",fontsize=7.2,color="white",fontweight="bold")
    ax.text(x0+width/2,y0-height*.17,label,ha="center",va="top",fontsize=7.2,color="#555555")

def mixing(ax,cx,y):
    ax.add_patch(FancyArrowPatch((cx+.095,y),(cx-.095,y),connectionstyle="arc3,rad=.55",arrowstyle="-|>",mutation_scale=7,linewidth=.8,color="#333333"))
    ax.add_patch(FancyArrowPatch((cx-.095,y-.005),(cx+.095,y-.005),connectionstyle="arc3,rad=.55",arrowstyle="-|>",mutation_scale=7,linewidth=.8,color="#333333"))
    ax.text(cx,y+.008,"20%",ha="center",va="center",fontsize=6.7)

def main():
    p=argparse.ArgumentParser(); p.add_argument("--project-root",default=None); p.add_argument("--output",default=None); a=p.parse_args()
    root=project_root(a.project_root); output=panel_output(root,"F",a.output); set_publication_style(7.5)
    fig,axes=plt.subplots(2,2,figsize=(5.2,3.95)); add_panel_letter(fig,"F",x=.01,y=.99,size=14)
    for ax,s in zip(axes.ravel(),SCENARIOS):
        ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off"); ax.text(.5,.96,s["model"],ha="center",va="top",fontsize=9,fontweight="bold"); ax.text(.5,.83,s["title"],ha="center",va="top",fontsize=7.5); ax.text(.5,.72,s["subtitle"],ha="center",va="top",fontsize=6.6,color="#666666")
        draw_cell(ax,.08,.30,.36,.31,s["left"],"left"); draw_cell(ax,.56,.30,.36,.31,s["right"],"right"); mixing(ax,.50,.18)
    fig.subplots_adjust(left=.06,right=.99,bottom=.05,top=.94,wspace=.10,hspace=.13); fig.savefig(output,bbox_inches="tight"); plt.close(fig); print(f"WROTE\t{output}")

if __name__=="__main__": main()
