#!/usr/bin/env python3
"""Focused Yb(d) comparison using cached coverage; guide anchors are not CRE bounds."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams.update({"pdf.fonttype":42,"font.size":9})
ANCHORS={
"Hmel2.5":{"contig":"Hmel215003o","primers":[1404112,1405075],
    "guides":{"Del2":[1404514,1404537],"Del3":[1404936,1404959],"Del4":[1404754,1404776]}},
"HmRos":{"contig":"NC_140403.1","primers":[1633657,1635225],
    "guides":{"Del1":[1633984,1634007],"Del2":[1634694,1634717]}}}
def main():
    p=argparse.ArgumentParser();p.add_argument("--run-dir",required=True,type=Path);a=p.parse_args()
    cfg=json.loads((a.run_dir/"run.json").read_text())
    out=a.run_dir/"Yb_d_Luca_guide_anchored_rosina_vs_melpomene_both_references.pdf"
    tmp=out.with_name(out.stem+".tmp.pdf")
    with PdfPages(tmp) as pdf:
        for column,title in [("normalized_depth","Normalized coverage (0–2)"),
                             ("raw_bigwig","Raw coverage: ln(1 + coverage)")]:
            fig,axes=plt.subplots(2,2,figsize=(12,7),sharey=True,sharex="col")
            fig.subplots_adjust(left=.10,right=.98,bottom=.18,top=.82,hspace=.3,wspace=.14)
            for col,(ref,anchor) in enumerate(ANCHORS.items()):
                tasks=[e for e in cfg["tasks"] if e["reference"]["name"]==ref]
                focal=[e for e in tasks if e["role"]=="focal"]
                for row,race in enumerate(["rosina","melpomene"]):
                    ax=axes[row,col]
                    controls=[e for e in tasks if e["role"]=="control" and e["species"]=="melpomene" and e["subspecies"]==race]
                    if len(controls)!=2:raise RuntimeError(f"Expected 2 {race} controls")
                    tables={}
                    lo,hi=anchor["primers"];left,right=lo-500,hi+500
                    for e in controls+focal:
                        t=pd.read_csv(a.run_dir/ref/"profiles"/(e["sample"]+".tsv.gz"),sep="\t")
                        tables[e["sample"]]=t.loc[t.position.between(left,right)].copy()
                    x=tables[controls[0]["sample"]].position.to_numpy()
                    if not len(x):raise RuntimeError("No cached bins in focused range")
                    for t in tables.values():
                        if not np.array_equal(t.position.to_numpy(),x):raise RuntimeError("Mismatched bins")
                    transform=np.log1p if column=="raw_bigwig" else lambda v:v
                    vals=np.vstack([tables[e["sample"]][column].to_numpy() for e in controls])
                    # Thin control individuals expose within-race variability; thick line is their mean.
                    for e,v in zip(controls,vals):
                        ax.plot(x/1e6,transform(v),color="#888888",lw=.7,alpha=.5)
                    ax.plot(x/1e6,transform(vals.mean(axis=0)),color="black",lw=1.8,label="Control mean (n=2)")
                    for e in focal:
                        color="#1764AB" if e["side"].lower()=="left" else "#D04B35"
                        ax.plot(x/1e6,transform(tables[e["sample"]][column].to_numpy()),color=color,
                            ls="--" if e["leg"]=="T3" else "-",lw=1.1,label=e["side"]+" "+e["leg"])
                    ax.axvspan((lo+.5)/1e6,(hi+.5)/1e6,color="#E6E1D6",alpha=.55,zorder=0)
                    for name,(gs,ge) in anchor["guides"].items():
                        mid=(gs+ge)/2e6
                        ax.axvline(mid,color="#787878",ls=":",lw=.6)
                        if row==0:
                            ax.text(mid,1.01,name,transform=ax.get_xaxis_transform(),ha="center",va="bottom",fontsize=7)
                    ax.set_xlim(left/1e6,right/1e6)
                    if column=="normalized_depth":
                        ax.set_ylim(0,2);ax.set_yticks([0,.5,1,1.5,2])
                        ax.axhline(1,color="#999999",ls=":",lw=.6)
                    ax.set_title(f"{ref} | H. m. {race} | yellow bar "+("present" if race=="rosina" else "absent"),loc="left",pad=25)
                    ax.ticklabel_format(axis="x",style="plain",useOffset=False)
                    ax.spines[["top","right"]].set_visible(False)
                    ax.text(.01,.025,", ".join(e["sample"] for e in controls),transform=ax.transAxes,fontsize=7)
                axes[-1,col].set_xlabel(anchor["contig"]+" position (Mb)")
            handles,labels=axes[0,0].get_legend_handles_labels()
            fig.legend(handles,labels,ncol=5,loc="upper center",bbox_to_anchor=(.54,.935))
            fig.suptitle("Yb(d): Luca guide-anchored region | rosina versus melpomene\n"+title,y=.99,fontsize=12)
            fig.supylabel("50-bp median / scaffold mean" if column=="normalized_depth" else "ln(1 + raw bigWig coverage)",x=.025)
            fig.text(.10,.065,"Shading = primer-bounded reference span, not a verified ATAC peak or exact functional CRE boundary; 500 bp flanks.\n"
                "Dotted lines = exact published guide/PAM matches. Unmatched guides are omitted, not inferred. Grey curves = individual controls.\n"
                "Blue/red = focal Left/Right; solid/dashed = T2/T3 libraries from one butterfly. Each race has two control individuals.\n"+
                ("Normalized: MAPQ ≥30, duplicate/secondary/QC-failed/unmapped exclusions; paired overlaps counted once; axis clipped at 2."
                 if column=="normalized_depth" else "Raw: existing unfiltered 50-bp bamCoverage values; mean control depth transformed after averaging."),fontsize=8)
            pdf.savefig(fig);plt.close(fig)
    tmp.replace(out)
    (a.run_dir/"Yb_d_Luca_focus.audit.json").write_text(json.dumps({
        "anchors_0based":ANCHORS,"functional_CRE_bounds_verified":False,
        "guide_source":"Livraghi2021 Figure4 source data2; exact match output job138462070",
        "coverage_recomputed":False,"pdf":str(out)},indent=2)+"\n")
    print("PDF\t"+str(out))
if __name__=="__main__":main()
