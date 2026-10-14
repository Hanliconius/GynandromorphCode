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
    outfile=a.outdir/"Yb_d_v_PK_homozygous_Ybplus_observed_LR_mixing_sensitivity.pdf"
    audit={"source_coverage":str(a.coverage),"genotype_provenance":"User rescored all 20 retained PK specimens as homozygous Yb+/Yb+ on 2026-10-14; not Optix B/D/R",
           "source_fit":str(a.fit_summary),"fitted_mixing":[cL,cR],
           "mixing_delta":a.mixing_delta,"mixing_envelope":"Range of inverted latent profiles for independent rates within +/- delta around saved Optix estimates. Not a confidence interval.",
           "states":states.to_dict("records"),"panels":[]}
    with PdfPages(outfile) as pdf:
        for scale in ("contig normalized","reference-control relative"):
            for ref in ("Hmel2.5","HmRos"):
                dataset=cov.loc[(cov.reference==ref)&(cov.scale==scale)]
                fig,axes=plt.subplots(2,2,figsize=(12,7),sharey="row")
                for j,interval in enumerate(("d","v")):
                    lo,hi=WINDOWS[ref][interval]; subset=dataset.loc[dataset.position.between(lo-1000,hi+1000)]
                    upper,lower=axes[:,j]
                    curves=subset.loc[subset["sample"].isin(pk)]
                    if curves.empty:raise RuntimeError("No Yb+/Yb+ PK coverage for "+ref+" "+interval)
                    for _,sample_data in curves.groupby("sample"):
                        sample_data=sample_data.sort_values("position")
                        upper.plot(sample_data.position/1e6,sample_data.observed,color=COLORS[1],alpha=.20,lw=.65)
                    means=curves.groupby("position").observed.mean()
                    upper.plot(means.index/1e6,means.values,color=COLORS[1],lw=2.2,label="PK Yb+/Yb+ (n=20)")
                    for group,color in (("rosina","#662A83"),("melpomene","#555555")):
                        ctrl=subset.loc[subset.group==group]
                        if not ctrl.empty:
                            mean=ctrl.groupby("position").observed.mean()
                            upper.plot(mean.index/1e6,mean.values,color=color,lw=1,ls="--",label=group)
                    arrays={}
                    for side in ("Left","Right"):
                        data=subset.loc[subset["sample"]=="TN_"+side+"_combined"].sort_values("position")
                        if data.empty:raise RuntimeError("Missing focal combined "+side+" for "+ref+" "+scale)
                        arrays[side]=data
                    L,R=arrays["Left"],arrays["Right"]
                    if not np.array_equal(L.position.to_numpy(),R.position.to_numpy()):raise RuntimeError("L/R grid mismatch")
                    x=L.position.to_numpy()/1e6
                    obsL=L.observed.to_numpy(dtype=float);obsR=R.observed.to_numpy(dtype=float)
                    for side in ("Left","Right"):
                        obs=obsL if side=="Left" else obsR
                        lower.plot(x,obs,color=SIDES[side],lw=1.8,label="Observed "+side)
                    aa=np.linspace(max(0,cL-a.mixing_delta),min(.49,cL+a.mixing_delta),21)
                    bb=np.linspace(max(0,cR-a.mixing_delta),min(.49,cR+a.mixing_delta),21)
                    latent={"Left":[],"Right":[]}
                    for ca in aa:
                        for cb in bb:
                            det=1-ca-cb
                            if det<=.05:continue
                            latent["Left"].append(((1-cb)*obsL-ca*obsR)/det)
                            latent["Right"].append(((1-ca)*obsR-cb*obsL)/det)
                    for side in ("Left","Right"):
                        if not latent[side]:raise RuntimeError("Mixing grid contains no valid values")
                        arr=np.vstack(latent[side])
                        valid=np.isfinite(arr)
                        with np.errstate(all="ignore"):
                            bottom=np.nanmin(np.where(valid,arr,np.nan),axis=0)
                            top=np.nanmax(np.where(valid,arr,np.nan),axis=0)
                        lower.fill_between(x,bottom,top,color=SIDES[side],alpha=.16,label=side+" sensitivity envelope")
                    for ax in (upper,lower):
                        ax.axvspan(lo/1e6,hi/1e6,color="grey",alpha=.10)
                        if interval=="d":
                            a0,b0=PRIMER[ref]
                            ax.axvspan(a0/1e6,b0/1e6,color="#693D99",alpha=.07)
                        ax.set_xlim((lo-1000)/1e6,(hi+1000)/1e6)
                        ax.set_ylim(0,2)
                        ax.spines[["top","right"]].set_visible(False)
                        ax.set_xlabel("Position (Mb)")
                    upper.set_title("Yb("+interval+") — homozygous Yb+ PK cohort")
                    if j==0:
                        upper.set_ylabel("PK coverage")
                        lower.set_ylabel("Combined Left/Right coverage")
                    audit["panels"].append({"reference":ref,"interval":interval,"scale":scale,
                                             "pk_Ybplus_homozygous":len(pk)})
                handles,labels=[],[]
                for ax in (axes[0,0],axes[1,0]):
                    h,l=ax.get_legend_handles_labels();handles+=h;labels+=l
                fig.legend(handles,labels,ncol=4,loc="upper center",bbox_to_anchor=(.5,.94),fontsize=8)
                fig.suptitle(ref+" | "+scale,fontsize=13,y=.99)
                fig.text(.08,.035,f"Solid L/R lines: measured T2+T3 combined coverage. Ribbons: linear-mixing sensitivity, not confidence intervals.\nMixing from saved Optix fit: L<-R={cL:.3f}, R<-L={cR:.3f}; varied independently +/-{a.mixing_delta:.2f}. 1-kb flanks. All 20 PK samples: user-rescored homozygous Yb+/Yb+; no Yb-negative PK group.",fontsize=8)
                fig.tight_layout(rect=(0,.085,1,.87))
                pdf.savefig(fig);plt.close(fig)
    (a.outdir/"audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    print("PDF",outfile)
    print("PK_YBPLUS_HOMOZYGOUS",len(pk))
if __name__=="__main__":main()
