#!/usr/bin/env python3
"""Figure 2B — mirrored chromosome-by-chromosome relative coverage."""

from __future__ import annotations
import argparse
from pathlib import Path
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from common import LEFT_T2, LEFT_T3, RIGHT_T2, RIGHT_T3, add_panel_letter, panel_output, project_root, set_publication_style

DEFAULT_INPUT = Path("analysis/chromosome_coverage_bars/run_20260923T162015902039Z/chromosome_coverage.tsv")

def chromosome_rank(label: str):
    text=str(label).strip().upper()
    if text.isdigit(): return (0,int(text),text)
    if text=="W": return (1,0,text)
    if text=="Z": return (1,1,text)
    return (2,0,text)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--project-root",default=None)
    parser.add_argument("--input",default=str(DEFAULT_INPUT))
    parser.add_argument("--output",default=None)
    args=parser.parse_args()
    root=project_root(args.project_root)
    input_path=Path(args.input)
    if not input_path.is_absolute(): input_path=root/input_path
    if not input_path.is_file(): raise FileNotFoundError(input_path)
    output=panel_output(root,"B",args.output)
    set_publication_style(7.5)
    df=pd.read_csv(input_path,sep="\t")
    required={"label","left_t2_relative_depth","left_t3_relative_depth","right_t2_relative_depth","right_t3_relative_depth"}
    missing=required-set(df.columns)
    if missing: raise RuntimeError(f"Coverage table missing columns: {sorted(missing)}")
    df=df.copy(); df["_rank"]=df["label"].map(chromosome_rank); df=df.sort_values("_rank").reset_index(drop=True)
    labels=df["label"].astype(str).tolist(); y=np.arange(len(df),dtype=float)
    values=np.concatenate([df[c].to_numpy(float) for c in ["left_t2_relative_depth","left_t3_relative_depth","right_t2_relative_depth","right_t3_relative_depth"]])
    finite=values[np.isfinite(values)]
    if not len(finite): raise RuntimeError("No finite relative-depth values")
    xmax=max(1.10,float(np.nanmax(finite))*1.06)
    fig=plt.figure(figsize=(3.55,5.35)); add_panel_letter(fig,"B",x=0.015,y=0.988,size=14)
    fig.text(0.19,0.965,"Chromosome",ha="left",va="top",fontsize=8.5)
    ax=fig.add_axes([0.11,0.10,0.82,0.82]); bar_h=0.15; offsets=(-0.10,+0.10)
    ax.barh(y+offsets[0],-df["left_t2_relative_depth"].to_numpy(float),height=bar_h,color=LEFT_T2,edgecolor="none",label="Left T2",zorder=2)
    ax.barh(y+offsets[1],-df["left_t3_relative_depth"].to_numpy(float),height=bar_h,color=LEFT_T3,edgecolor="none",label="Left T3",zorder=2)
    ax.barh(y+offsets[0],df["right_t2_relative_depth"].to_numpy(float),height=bar_h,color=RIGHT_T2,edgecolor="none",label="Right T2",zorder=2)
    ax.barh(y+offsets[1],df["right_t3_relative_depth"].to_numpy(float),height=bar_h,color=RIGHT_T3,edgecolor="none",label="Right T3",zorder=2)
    ax.axvspan(-0.105*xmax,0.105*xmax,color="white",zorder=3); ax.axvline(0,color="#444444",linewidth=0.65,zorder=4)
    for yi,label in zip(y,labels):
        colour="#B2182B" if label.upper() in {"W","Z"} else "#333333"
        weight="bold" if label.upper() in {"W","Z"} else "normal"
        ax.text(0,yi,label,ha="center",va="center",fontsize=6.6,color=colour,fontweight=weight,zorder=5)
    ax.axvline(-1.0,color="#8A8A8A",linewidth=0.55,linestyle=(0,(2,2)),zorder=1)
    ax.axvline(+1.0,color="#8A8A8A",linewidth=0.55,linestyle=(0,(2,2)),zorder=1)
    ax.set_xlim(-xmax,xmax); ax.set_ylim(-0.6,len(df)-0.4); ax.invert_yaxis(); ax.set_yticks([])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value,_: f"{abs(value):.1f}"))
    ax.set_xlabel("relative depth",labelpad=3); ax.tick_params(axis="x",length=2.5,pad=2)
    for s in ("top","right","left"): ax.spines[s].set_visible(False)
    fig.text(0.23,0.055,"left",ha="center",va="center",fontsize=7.5,color=LEFT_T2)
    fig.text(0.80,0.055,"right",ha="center",va="center",fontsize=7.5,color=RIGHT_T2)
    handles,legend_labels=ax.get_legend_handles_labels()
    fig.legend(handles,legend_labels,loc="lower center",bbox_to_anchor=(0.52,0.002),ncol=4,frameon=False,fontsize=5.7,handlelength=1.2,columnspacing=0.9)
    fig.savefig(output,bbox_inches="tight"); plt.close(fig)
    print(f"INPUT\t{input_path}"); print(f"WROTE\t{output}")

if __name__=="__main__": main()
