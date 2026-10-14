#!/usr/bin/env python3
"""Recalculate Yb coverage using Luca's 50-bp per-base median / scaffold mean.

Sources: the completed Yb all-samples run, its frozen sample manifests,
and existing BAMs. No downloading or alignment. Two vector PDFs:
race-specific tracks and phenotype-pooled non-PK control tracks.
"""
import argparse, json, subprocess, sys, shlex, shutil, hashlib
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd

ROOT=Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
RUN=ROOT/"analysis/yb_all_samples_coverage/20261014T194126Z"
REGIONS={"Hmel2.5":{"d":("Hmel215003o",1403500,1405500),"v":("Hmel215003o",1429500,1430500)},
         "HmRos":{"d":("NC_140403.1",1633041,1635646),"v":("NC_140403.1",1658079,1659072)}}
STATES={"rosina":1,"melpomene":0,"amaryllis":1,"thelxinoe":1,"florencia":0}
RACES=["rosina","melpomene","PK Yb+/Yb+","amaryllis","thelxinoe","florencia"]
PALETTE={"rosina":"#662A83","melpomene":"#555555","PK Yb+/Yb+":"#693D99",
         "amaryllis":"#3A8A70","thelxinoe":"#CC9332","florencia":"#AE668C",
         "Yb+ controls":"#24775E","Yb- controls":"#A35B47"}
SIDES={"Left":"#1764AB","Right":"#D9472D"}
FRAG_FLAGS="1796"

def need(path):
    p=Path(path)
    if not p.is_file() or p.stat().st_size<=0:
        raise RuntimeError("Missing or empty: "+str(p))
    return p

def records(run):
    audit=json.loads(need(run/"audit.json").read_text())
    cfg=json.loads(need(run/"run.json").read_text())
    origins=cfg["sources"]
    source_cfg=[json.loads(need(Path(p)/"run.json").read_text()) for p in origins.values()]
    bam_index={}
    for source in source_cfg:
        for entry in source.get("tasks",[]):
            ref=entry.get("reference",{}).get("name")
            sample=entry.get("sample")
            path=entry.get("bam")
            if path and ref and sample:
                bam_index[(ref,sample)]=Path(path)
    samples=cfg["samples"]
    if len(samples)!=76:raise RuntimeError("Expected 38 per reference in retained-PK analysis")
    tasks=[]
    for e in samples:
        ref=e["reference"]; sample=e["sample"]
        group=e["group"]
        if group.startswith("PK "):
            group="PK Yb+/Yb+"
            base=Path(origins["pk_run"])/ref/sample
            hits=sorted(p for p in base.rglob("*.bam") if not p.name.endswith(".bai"))
            # Some PK manifests point to a BAM outside the per-sample directory.
            if not hits:
                pc=source_cfg[1]
                for entry in pc.get("tasks",[]):
                    if entry.get("sample")==sample and entry.get("reference",{}).get("name")==ref:
                        hits.append(Path(entry["bam"]))
            if len(hits)!=1:
                raise RuntimeError(f"Need exactly one PK BAM {ref} {sample}; found {len(hits)} in {base}")
            bam=hits[0]
        else:
            bam=bam_index.get((ref,sample))
            if bam is None:
                base=Path(origins["mean_run"])/ref
                hits=sorted(base.rglob(sample+"*.bam"))
                if len(hits)==1:bam=hits[0]
            if bam is None:raise RuntimeError(f"Cannot resolve BAM for {ref} {sample}")
        need(bam)
        profile=pd.read_csv(need(e["profile"]),sep="\t",usecols=["mean_depth","normalized_depth"])
        valid=np.isfinite(profile["mean_depth"]) & np.isfinite(profile["normalized_depth"]) & (profile["normalized_depth"]>0)
        ratio=(profile.loc[valid,"mean_depth"]/profile.loc[valid,"normalized_depth"]).to_numpy(dtype=float)
        if not len(ratio):raise RuntimeError(f"No valid scaffold denominator: {ref} {sample}")
        den=float(np.median(ratio))
        if not np.allclose(ratio,den,rtol=1e-5,atol=1e-5):
            raise RuntimeError(f"Inconsistent scaffold denominator: {ref} {sample}")
        tasks.append(dict(ref=ref,sample=sample,group=group,bam=str(bam),denominator=den,
                          provenance="user-rescored homozygous Yb+/Yb+" if group=="PK Yb+/Yb+" else "repository race Yb phenotype"))
    # Focal side means must be reconstructed by pooling base-level depth, not averaging per-library medians.
    for ref in REGIONS:
        for side in ("Left","Right"):
            es=[e for e in samples if e["reference"]==ref and e["group"]=="focal" and e.get("side")==side]
            if len(es)!=2 or {e.get("leg") for e in es}!={"T2","T3"}:
                raise RuntimeError(f"Missing focal T2/T3 for {ref} {side}")
            orig=[t for t in tasks if t["ref"]==ref and t["sample"] in {e["sample"] for e in es}]
            if len(orig)!=2:raise RuntimeError("Focal profile mapping mismatch")
            tasks.append(dict(ref=ref,sample="TN_"+side+"_combined",group="focal combined",
                              bam=[t["bam"] for t in orig],denominator=sum(t["denominator"] for t in orig),
                              provenance="per-base pooled T2+T3 depth; summed scaffold means"))
    for ref in REGIONS:
        rt=[t for t in tasks if t["ref"]==ref]
        assert len(rt)==40,(ref,len(rt))
        for group,n in (("rosina",4),("melpomene",4),("PK Yb+/Yb+",20),("amaryllis",2),("thelxinoe",2),("florencia",2),("focal combined",2)):
            got=sum(t["group"]==group for t in rt)
            if got!=n:raise RuntimeError(f"{ref}: {group}: expected {n}, got {got}")
    return tasks, origins

def submit(a):
    base=ROOT/"analysis/yb_luca_median_replication"/datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    tasks,origins=records(a.run)
    base.mkdir(parents=True)
    (base/"logs").mkdir()
    script=base/Path(__file__).name
    shutil.copy2(__file__,script)
    cfg={"run_dir":str(base),"source_run":str(a.run),"sources":origins,"samtools":shutil.which("samtools"),
         "method":"samtools depth -aa -s -q 0 -Q 30 -G1796; per-base depth median in 50bp absolute reference bins; divide by original own full-scaffold mean",
         "comparison_to_Luca":"Luca used mosdepth; using samtools with equivalent intended MAPQ>=30 and per-base window medians, but not byte-identical mosdepth output",
         "phenotype_labels":STATES,"pk_genotype":"20 PK are homozygous Yb+/Yb+ by user rescoring",
         "tasks":tasks}
    if not cfg["samtools"]:raise RuntimeError("samtools must be on PATH")
    (base/"run.json").write_text(json.dumps(cfg,indent=2)+"\n")
    py=shlex.quote(sys.executable)
    runner=" ".join((py,shlex.quote(str(script)),"--config",shlex.quote(str(base/"run.json"))))
    common=["sbatch","--parsable","--partition=cpu","--chdir="+str(base)]
    def job(name,wrap,array=None,dependency=None):
        flags=["--job-name="+name,"--time=04:00:00","--mem=4G","--cpus-per-task=1"]
        if array:flags+=["--array="+array]
        if dependency:flags+=["--dependency=afterok:"+dependency]
        flags+=["--output="+str(base/"logs"/(name+".%A_%a.out" if array else name+".%j.out")),
                "--error="+str(base/"logs"/(name+".%A_%a.err" if array else name+".%j.err"))]
        return subprocess.check_output(common+flags+["--wrap="+wrap],text=True).strip().split(";")[0]
    work=job("YbLucaMedian",runner+" --stage worker --task \"$SLURM_ARRAY_TASK_ID\"",f"0-{len(tasks)-1}%20")
    plot=job("YbLucaPlot",runner+" --stage report",dependency=work)
    (base/"jobs.json").write_text(json.dumps({"median_array":work,"report":plot},indent=2)+"\n")
    print("MEDIAN_ARRAY",work)
    print("REPORT_JOB",plot)
    print("RUN_DIR",base)

def worker(cfg,i):
    t=cfg["tasks"][i]; ref=t["ref"]
    samtools=cfg["samtools"]
    bams=t["bam"] if isinstance(t["bam"],list) else [t["bam"]]
    out=Path(cfg["run_dir"])/ref/"median"/(t["sample"]+".tsv.gz")
    out.parent.mkdir(parents=True,exist_ok=True)
    rows=[]
    for key,(contig,lo,hi) in REGIONS[ref].items():
        start=max(0,lo-1000);end=hi+1000
        # The original analysis uses 50-bp windows on an absolute contig grid.
        first=(start//50)*50
        last=((end+49)//50)*50
        cmd=[samtools,"depth","-aa","-s","-q","0","-Q","30","-G",FRAG_FLAGS,
             "-r",f"{contig}:{first+1}-{last}",*map(str,bams)]
        proc=subprocess.Popen(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        per_bin={}
        last_pos=None
        for line in proc.stdout:
            fields=line.rstrip("\n").split("\t")
            pos=int(fields[1])-1
            if not first<=pos<last:continue
            if last_pos is not None and pos<=last_pos:raise RuntimeError("Unsorted/duplicated samtools depth")
            last_pos=pos
            depth=sum(map(int,fields[2:]))
            per_bin.setdefault(pos//50,[]).append(depth)
        stderr=proc.stderr.read()
        if proc.wait():raise RuntimeError(f"samtools depth failed {t['sample']}: {stderr[-1200:]}")
        for bs in range(first,last,50):
            values=per_bin.get(bs//50,[])
            if len(values)!=50:
                raise RuntimeError(f"Incorrect zero-inclusive depth for {ref} {t['sample']} {bs}: {len(values)} bases")
            median=float(np.median(values))
            rows.append({"reference":ref,"sample":t["sample"],"group":t["group"],"interval":key,
                         "position":bs+25.0,"median_depth":median,
                         "normalized_depth":median/t["denominator"]})
    pd.DataFrame(rows).to_csv(out,sep="\t",index=False)
    print("DONE",ref,t["sample"])

def report(cfg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams.update({"pdf.fonttype":42,"font.size":9})
    root=Path(cfg["run_dir"])
    allrows=[]
    for t in cfg["tasks"]:
        path=root/t["ref"]/"median"/(t["sample"]+".tsv.gz")
        tab=pd.read_csv(need(path),sep="\t")
        allrows.append(tab)
    data=pd.concat(allrows,ignore_index=True)
    data.to_csv(root/"median_coverage_values.tsv.gz",sep="\t",index=False)
    fit=ROOT/"analysis/original_502/reciprocal_L_R_mixture/optix_reciprocal_L_R_mixture.summary.tsv"
    params=pd.read_csv(need(fit),sep="\t").set_index("statistic").value
    cl=float(params["estimated_c_left_from_right"])
    cr=float(params["estimated_c_right_from_left"])
    for mode in ("separate_races","pooled_nonPK_phenotype"):
        path=root/("Yb_Luca50bpMedian_"+mode+".pdf")
        with PdfPages(path) as pdf:
            for ref in REGIONS:
                groups=(RACES if mode=="separate_races" else ["Yb+ controls","Yb- controls","PK Yb+/Yb+"])
                nrow=len(groups)+1
                fig,axes=plt.subplots(nrow,2,figsize=(13,max(8,nrow*2.1)),sharey=True)
                for j,key in enumerate(("d","v")):
                    subset=data[(data["reference"]==ref)&(data["interval"]==key)]
                    contig,lo,hi=REGIONS[ref][key]
                    for i,g in enumerate(groups):
                        ax=axes[i,j]
                        if g=="Yb+ controls":ids=["rosina","amaryllis","thelxinoe"]
                        elif g=="Yb- controls":ids=["melpomene","florencia"]
                        else:ids=[g]
                        traces=subset[subset["group"].isin(ids)]
                        names=traces["sample"].unique()
                        if not len(names):raise RuntimeError("Empty "+g)
                        for _,d in traces.groupby("sample"):
                            ax.plot(d["position"]/1e6,d["normalized_depth"],color=PALETTE[g],lw=.65,alpha=.3)
                        avg=traces.groupby("position")["normalized_depth"].mean()
                        ax.plot(avg.index/1e6,avg.values,color=PALETTE[g],lw=2,label=f"{g} (n={len(names)})")
                        if j==0:ax.set_ylabel(g+f" (n={len(names)})")
                        if j==1:ax.legend(loc="upper right",fontsize=7,frameon=False)
                    ax=axes[-1,j]
                    lr={}
                    for side in SIDES:
                        d=subset[subset["sample"]=="TN_"+side+"_combined"].sort_values("position")
                        lr[side]=d["normalized_depth"].to_numpy()
                        if side=="Left":x=d["position"].to_numpy()/1e6
                        elif not np.array_equal(d["position"].to_numpy()/1e6,x):raise RuntimeError("LR grid mismatch")
                    aa=np.linspace(max(0,cl-.05),min(.49,cl+.05),15)
                    bb=np.linspace(max(0,cr-.05),min(.49,cr+.05),15)
                    samples={"Left":[],"Right":[]}
                    L,R=lr["Left"],lr["Right"]
                    for a in aa:
                        for b in bb:
                            if 1-a-b<=.05:continue
                            samples["Left"].append(((1-b)*L-a*R)/(1-a-b))
                            samples["Right"].append(((1-a)*R-b*L)/(1-a-b))
                    for side,color in SIDES.items():
                        arr=np.array(samples[side])
                        ok=np.isfinite(arr).any(axis=0)
                        lower=np.full(len(x),np.nan);upper=lower.copy()
                        if ok.any():
                            lower[ok]=np.nanmin(arr[:,ok],axis=0)
                            upper[ok]=np.nanmax(arr[:,ok],axis=0)
                        ax.fill_between(x,lower,upper,color=color,alpha=.16)
                        ax.plot(x,lr[side],color=color,lw=1.7,label="Observed "+side)
                    if j==0:ax.set_ylabel("Combined Left / Right")
                    ax.legend(loc="upper right",fontsize=7,frameon=False)
                    for axi in axes[:,j]:
                        axi.axvspan(lo/1e6,hi/1e6,color="#BDBDBD",alpha=.14)
                        axi.set_xlim((lo-1000)/1e6,(hi+1000)/1e6)
                        axi.spines[["top","right"]].set_visible(False)
                        axi.tick_params(labelbottom=True,labelsize=7)
                    axes[0,j].set_title("Yb("+key+")")
                    axes[-1,j].set_xlabel(contig+" position (Mb)")
                fig.suptitle(ref+" | Luca-style 50-bp median / full-scaffold mean | "+mode.replace("_"," "),fontsize=12)
                fig.text(.09,.012,f"MAPQ>=30, 50-bp median of zero-inclusive per-base depth; full-scaffold mean from existing validated coverage. "
                         f"20 PK Yb+/Yb+ (user-rescored). L/R ribbons: Optix mixing fit {cl:.3f}/{cr:.3f} ±0.05, sensitivity only.",fontsize=8)
                fig.tight_layout(rect=(0,.035,1,.96),h_pad=.6)
                pdf.savefig(fig);plt.close(fig)
        print("PDF",path)
    (root/"audit.json").write_text(json.dumps({"method":cfg["method"],"comparison_to_Luca":cfg["comparison_to_Luca"],
                    "controls":cfg["phenotype_labels"],"pk":cfg["pk_genotype"],"fit":str(fit),
                    "sources":cfg["sources"],"rows":len(data)},indent=2)+"\n")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--run",type=Path,default=RUN)
    p.add_argument("--config",type=Path)
    p.add_argument("--stage",choices=["worker","report"])
    p.add_argument("--task",type=int)
    a=p.parse_args()
    if not a.config:submit(a)
    else:
        cfg=json.loads(need(a.config).read_text())
        if a.stage=="worker":worker(cfg,a.task)
        elif a.stage=="report":report(cfg)
        else:raise RuntimeError("--stage required")
if __name__=="__main__":main()
