#!/usr/bin/env python3
"""Focused Yb d/v coverage for 20 user-rescored homozygous Yb+/Yb+ PK samples."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

WINDOWS = {"Hmel2.5":{"d":(1403500,1405500),"v":(1429500,1430500)},
           "HmRos":{"d":(1633041,1635646),"v":(1658079,1659072)}}
PRIMER = {"Hmel2.5":(1404112,1405075),"HmRos":(1633657,1635225)}
COLORS = {0:"#C06C9A",1:"#693D99"}
SIDES = {"Left":"#1764AB","Right":"#D9472D"}

def binary(value):
    s=str(value).strip()
    if s in ("0","0.0"):return 0
    if s in ("1","1.0"):return 1
    raise ValueError("Unverified Yb state (expected 0/1): "+repr(value))

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--coverage",required=True,type=Path)
    p.add_argument("--fit-summary",required=True,type=Path)
    p.add_argument("--outdir",required=True,type=Path)
    p.add_argument("--mixing-delta",type=float,default=.05)
    a=p.parse_args()
    if not (0<=a.mixing_delta<=.2):raise ValueError("mixing delta must be between 0 and .2")
    cov=pd.read_csv(a.coverage,sep="\t")
    required={"reference","sample","group","position","scale","observed"}
    if not required.issubset(cov.columns):raise RuntimeError("Coverage table missing "+str(required-set(cov.columns)))
    pk=sorted(cov.loc[cov.group.astype(str).str.startswith("PK "),"sample"].unique())
    if len(pk)!=20:raise RuntimeError("Expected 20 retained PK specimens; found "+str(len(pk)))
    # User rescoring (2026-10-14): every retained PK is homozygous Yb+/Yb+.
    # This is independent of Optix B/D/R and PK dark/pale phenotypes.
    states=pd.DataFrame([{"sample":sample,"Yb_genotype":"Yb+/Yb+","state":1,"source":"user rescoring 2026-10-14"} for sample in pk])
    fit=pd.read_csv(a.fit_summary,sep="\t").set_index("statistic").value
    cL=float(fit["estimated_c_left_from_right"])
    cR=float(fit["estimated_c_right_from_left"])
    if not (0<=cL<.5 and 0<=cR<.5 and cL+cR<1):raise RuntimeError("Invalid fitted mixing fractions")
    a.outdir.mkdir(parents=True,exist_ok=True)
    states.to_csv(a.outdir/"pk_homozygous_Ybplus_states_used.tsv",sep="\t",index=False)
    plt.rcParams.update({"pdf.fonttype":42,"font.size":9})
    outfile=a.outdir/"Yb_d_v_separate_race_tracks_observed_LR.pdf"
    audit={"source_coverage":str(a.coverage),"genotype_provenance":"User rescored all 20 retained PK specimens as homozygous Yb+/Yb+ on 2026-10-14; not Optix B/D/R",
           "source_fit":str(a.fit_summary),"fitted_mixing":[cL,cR],
           "mixing_delta":a.mixing_delta,"mixing_envelope":"Range of inverted latent profiles for independent rates within +/- delta around saved Optix estimates. Not a confidence interval.",
           "states":states.to_dict("records"),"panels":[]}
    race_groups=[
        ("rosina","rosina",4,"#662A83"),
        ("melpomene","melpomene",4,"#555555"),
        ("PK Yb+/Yb+","PK",20,"#693D99"),
        ("amaryllis","amaryllis",2,"#35856D"),
        ("thelxinoe","thelxinoe",2,"#D09537"),
        ("florencia","florencia",2,"#A66B94"),
    ]
    with PdfPages(outfile) as pdf:
        for ref in ("Hmel2.5","HmRos"):
            dataset=cov.loc[(cov["reference"]==ref)&(cov["scale"]=="contig normalized")]
            if dataset.empty:raise RuntimeError("Missing contig-normalized coverage for "+ref)
            fig,axes=plt.subplots(len(race_groups)+1,2,figsize=(13,16),sharey=True)
            for j,interval in enumerate(("d","v")):
                lo,hi=WINDOWS[ref][interval]
                subset=dataset.loc[dataset["position"].between(lo-1000,hi+1000)]
                for i,(label,group,n,color) in enumerate(race_groups):
                    ax=axes[i,j]
                    if group=="PK":
                        traces=subset.loc[subset["sample"].isin(pk)]
                    else:
                        traces=subset.loc[subset["group"]==group]
                    names=sorted(traces["sample"].unique())
                    if len(names)!=n:
                        raise RuntimeError(f"{ref} {interval}: {label}: expected {n}, got {len(names)}")
                    for _,trace in traces.groupby("sample"):
                        trace=trace.sort_values("position")
                        ax.plot(trace["position"]/1e6,trace["observed"],color=color,lw=.7,alpha=.3)
                    mean=traces.groupby("position")["observed"].mean()
                    ax.plot(mean.index/1e6,mean.values,color=color,lw=2,label=f"{label} mean (n={n})")
                    if j==0:ax.set_ylabel(f"{label} (n={n})")
                    if j==1:ax.legend(loc="upper right",fontsize=7,frameon=False)
                lower=axes[-1,j]
                arrays={}
                for side in ("Left","Right"):
                    data=subset.loc[subset["sample"]=="TN_"+side+"_combined"].sort_values("position")
                    if data.empty:raise RuntimeError(f"Missing combined {side}: {ref} {interval}")
                    arrays[side]=data
                L,R=arrays["Left"],arrays["Right"]
                if not np.array_equal(L["position"].to_numpy(),R["position"].to_numpy()):
                    raise RuntimeError("Combined L/R position mismatch")
                x=L["position"].to_numpy()/1e6
                obsL=L["observed"].to_numpy(dtype=float)
                obsR=R["observed"].to_numpy(dtype=float)
                aa=np.linspace(max(0,cL-a.mixing_delta),min(.49,cL+a.mixing_delta),21)
                bb=np.linspace(max(0,cR-a.mixing_delta),min(.49,cR+a.mixing_delta),21)
                latent={"Left":[],"Right":[]}
                for ca in aa:
                    for cb in bb:
                        det=1-ca-cb
                        if det<=.05:continue
                        latent["Left"].append(((1-cb)*obsL-ca*obsR)/det)
                        latent["Right"].append(((1-ca)*obsR-cb*obsL)/det)
                for side,obs in (("Left",obsL),("Right",obsR)):
                    if not latent[side]:raise RuntimeError("No valid mixing grid")
                    arr=np.vstack(latent[side])
                    good=np.isfinite(arr).any(axis=0)
                    bottom=np.full(arr.shape[1],np.nan)
                    top=np.full(arr.shape[1],np.nan)
                    if good.any():
                        selected=arr[:,good]
                        bottom[good]=np.nanmin(selected,axis=0)
                        top[good]=np.nanmax(selected,axis=0)
                    lower.fill_between(x,bottom,top,color=SIDES[side],alpha=.16)
                    lower.plot(x,obs,color=SIDES[side],lw=1.8,label="Observed "+side)
                if j==0:lower.set_ylabel("Combined Left/Right")
                if j==1:lower.legend(loc="upper right",fontsize=7,frameon=False)
                for ax in axes[:,j]:
                    ax.axvspan(lo/1e6,hi/1e6,color="#BDBDBD",alpha=.15)
                    if interval=="d":
                        p0,p1=PRIMER[ref]
                        ax.axvspan(p0/1e6,p1/1e6,color="#693D99",alpha=.08)
                    ax.set_xlim((lo-1000)/1e6,(hi+1000)/1e6)
                    ax.set_ylim(0,2)
                    ax.spines[["top","right"]].set_visible(False)
                    ax.tick_params(labelbottom=True,labelsize=7)
                axes[0,j].set_title("Yb("+interval+")",fontsize=12)
                axes[-1,j].set_xlabel("Position (Mb)")
                audit["panels"].append({"reference":ref,"interval":interval,
                                       "scale":"contig normalized",
                                       "groups":{name:n for name,_,n,_ in race_groups}})
            fig.suptitle(ref+" | Yb(d) and Yb(v) | race-specific tracks",fontsize=14)
            fig.text(.09,.016,
                f"50-bp mean depth / full-contig mean; 1-kb flanks. PK cohort: 20 user-rescored homozygous Yb+/Yb+. "
                f"Observed L/R pooled T2+T3; ribbons: mixing sensitivity (not confidence intervals), "
                f"L<-R={cL:.3f}, R<-L={cR:.3f} +/-{a.mixing_delta:.2f}.",fontsize=8)
            fig.tight_layout(rect=(0,.032,1,.974),h_pad=.65,w_pad=1.2)
            pdf.savefig(fig)
            plt.close(fig)
    audit["normalization"]="contig normalized only; reference-control-relative pages removed"
    audit["race_rows"]=[name for name,_,_,_ in race_groups]
    (a.outdir/"audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    print("PDF",outfile)
    print("PK_YBPLUS_HOMOZYGOUS",len(pk))
if __name__=="__main__":main()
