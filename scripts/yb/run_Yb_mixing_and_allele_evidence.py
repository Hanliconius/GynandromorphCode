#!/usr/bin/env python3
"""Run linear-depth mixing sensitivity and local CRE allele evidence preparation."""
import argparse,csv,importlib.util,json,os,shlex,shutil,subprocess,sys,uuid
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
ANCHORS={"Hmel2.5":("Hmel215003o",1404112,1405075),"HmRos":("NC_140403.1",1633657,1635225)}
def run(cmd,**kwargs):subprocess.run([str(x) for x in cmd],check=True,**kwargs)
def submit(source):
    cfg=json.loads((source/"run.json").read_text())
    for e in cfg["tasks"]:
        for path in (e["bam"],e["bam"]+".bai"):
            if not Path(path).is_file():raise RuntimeError("Missing "+path)
    root=source.parent.parent/"yb_mixing_and_allele_evidence"
    base=root/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8])
    (base/"logs").mkdir(parents=True)
    code=base/Path(__file__).name;shutil.copy2(__file__,code)
    here=Path(__file__).resolve().parent
    mean=(here/"Yb_control_overlays.py").read_text().replace("np.median(depth[","np.mean(depth[")
    mean=mean.replace("median_depth","mean_depth")
    (base/"mean_worker.py").write_text(mean)
    shutil.copy2(here/"livraghi2021_fig8_sanger_sequences.fasta",base/"paper_sanger.fasta")
    cfg.update(run_dir=str(base),source_run=str(source),method="50-bp mean depth / full-contig mean; samtools depth MAPQ30",
        mixing_assumption="symmetric 20% in both legs; sensitivity 15%,25%; NOT fitted estimates")
    config=base/"run.json";config.write_text(json.dumps(cfg,indent=2)+"\n")
    worker=base/"worker.sh"
    worker.write_text("#!/usr/bin/env bash\nset -euo pipefail\nexec "+shlex.join([sys.executable,str(code),"--config",str(config)])+' --stage "$1" --task "$SLURM_ARRAY_TASK_ID"\n')
    jobs={}
    try:
        for stage,array,cpus,mem,limit in [("mean","0-27%28",2,"8G","04:00:00"),("alleles","0-1%2",4,"8G","02:00:00"),("report","0-0%1",1,"8G","00:30:00")]:
            cmd=["sbatch","--parsable","--partition=cpu",f"--job-name=Yb_{stage}",f"--array={array}",
                f"--cpus-per-task={cpus}",f"--mem={mem}",f"--time={limit}",f"--chdir={base}",
                f"--output={base}/logs/{stage}.%A_%a.out",f"--error={base}/logs/{stage}.%A_%a.err"]
            if stage=="report":cmd.append("--dependency=afterok:"+jobs["mean"])
            jobs[stage]=subprocess.check_output(cmd+[str(worker),stage],text=True).strip().split(";")[0]
            (base/"jobs.json").write_text(json.dumps(jobs,indent=2)+"\n")
            print(stage.upper()+"_ARRAY\t"+jobs[stage],flush=True)
    except Exception:
        for job in jobs.values():subprocess.run(["scancel",job],check=False)
        raise
    print("RUN_DIR\t"+str(base))
def mean(cfg,task):
    spec=importlib.util.spec_from_file_location("mean_worker",Path(cfg["run_dir"])/"mean_worker.py")
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);mod.coverage(cfg,task)
def report(cfg,relative=False,bundle=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams.update({"pdf.fonttype":42,"font.size":9})
    base=Path(cfg["run_dir"]);out=base/("Yb_d_control_relative_mean_coverage.pdf" if relative else "Yb_d_both_references_observed_vs_assumed20pct_unmixed_mean_coverage.pdf")
    allrows=[]
    with PdfPages(out) as pdf:
        fig,axes=plt.subplots(2,2,figsize=(12,7),sharey=True,sharex="col")
        fig.subplots_adjust(left=.09,right=.98,top=.83,bottom=.2,hspace=.25)
        for col,(ref,(contig,lo,hi)) in enumerate(ANCHORS.items()):
            entries=[e for e in cfg["tasks"] if e["reference"]["name"]==ref]
            profiles={}
            for e in entries:
                t=pd.read_csv(base/ref/"profiles"/(e["sample"]+".tsv.gz"),sep="\t")
                profiles[e["sample"]]=t.loc[t.position.between(lo-500,hi+500)].copy()
            x=next(iter(profiles.values())).position.to_numpy()
            for t in profiles.values():
                if not np.array_equal(x,t.position.to_numpy()):raise RuntimeError("Coordinate mismatch")
            if relative:
                baseline="melpomene" if ref=="Hmel2.5" else "rosina"
                ids=[e["sample"] for e in entries if e["role"]=="control" and e.get("species")=="melpomene" and e.get("subspecies")==baseline]
                denominator=np.vstack([profiles[s].normalized_depth.to_numpy() for s in ids]).mean(axis=0)
                valid=np.isfinite(denominator)&(denominator>=.25)
                for table in profiles.values():
                    table["normalized_depth"]=np.divide(table.normalized_depth.to_numpy(),denominator,
                        out=np.full(denominator.shape,np.nan),where=valid)
            for race,color in [("rosina","#662A83"),("melpomene","#555555")]:
                samples=[e["sample"] for e in entries if e["role"]=="control" and e["species"]=="melpomene" and e["subspecies"]==race]
                y=np.vstack([profiles[s].normalized_depth.to_numpy() for s in samples]).mean(axis=0)
                for ax in axes[:,col]:ax.plot(x/1e6,y,color=color,lw=1.5,label=race+" controls (n=2)")
            for leg,style in [("T2","-"),("T3","--")]:
                pair={e["side"].lower():e for e in entries if e["role"]=="focal" and e["leg"]==leg}
                L=profiles[pair["left"]["sample"]].normalized_depth.to_numpy()
                R=profiles[pair["right"]["sample"]].normalized_depth.to_numpy()
                for side,color,obs in [("Left","#1764AB",L),("Right","#D04B35",R)]:
                    unmix=lambda m:((1-m)*obs-m*(R if side=="Left" else L))/(1-2*m)
                    central=unmix(.2);a=unmix(.15);b=unmix(.25)
                    axes[0,col].plot(x/1e6,obs,color=color,ls=style,lw=1,label=side+" "+leg)
                    axes[1,col].plot(x/1e6,central,color=color,ls=style,lw=1,label=side+" "+leg)
                    axes[1,col].fill_between(x/1e6,np.minimum(a,b),np.maximum(a,b),color=color,alpha=.10)
                    negative=central<0
                    axes[1,col].scatter(x[negative]/1e6,np.zeros(negative.sum()),marker="v",s=12,color=color,zorder=5)
                    for pos,o,u,v,w in zip(x,obs,central,a,b):
                        allrows.append(dict(reference=ref,leg=leg,side=side,position=pos,observed=o,
                            assumed20pct_unmixed=u,sensitivity15pct=v,sensitivity25pct=w,negative=u<0))
            for row in range(2):
                ax=axes[row,col];ax.axvspan(lo/1e6,hi/1e6,color="#E6E1D6",alpha=.4)
                ax.set_ylim(0,2);ax.set_xlim((lo-500)/1e6,(hi+500)/1e6)
                ax.spines[["top","right"]].set_visible(False)
                ax.set_title(ref+" | "+("Observed" if row==0 else "Assumed 20% symmetric mixing corrected"),loc="left")
                ax.ticklabel_format(axis="x",style="plain",useOffset=False)
            axes[-1,col].set_xlabel(contig+" position (Mb)")
        h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,ncol=3,loc="upper center",bbox_to_anchor=(.54,.96))
        fig.suptitle("Yb(d): linear mean coverage and mixing sensitivity",fontsize=12)
        fig.supylabel("Coverage / reference-matched control mean" if relative else "50-bp mean / full-contig mean",x=.02)
        fig.text(.09,.05,"20% is a user-specified exploratory assumption, not the fitted per-library mixture. Shading around corrected curves spans 15–25%, not a confidence interval.\n"
            "Each library scaled by its own contig mean before inversion; assumes comparable underlying contig dosage and mapping response.\n"
            "MAPQ >=30; duplicate/secondary/QC-failed/unmapped excluded; paired overlap counted once. No log transform or median inversion.\n"
            "Negative estimates retained in TSV and marked by downward triangles at zero; display range 0–2. Primer span shaded; exact CRE boundaries pending.",fontsize=8)
        if relative:fig.text(.09,.015,"Per-bin baseline: Hmel2.5 = melpomene mean; HmRos = rosina mean. Baseline <0.25 scaffold-normalized depth masked; retained baseline = 1.",fontsize=8)
        pdf.savefig(fig)
        if bundle is not None:bundle.savefig(fig)
        plt.close(fig)
    pd.DataFrame(allrows).to_csv(base/("control_relative_unmixing_sensitivity.tsv.gz" if relative else "mean_coverage_unmixing_sensitivity.tsv.gz"),sep="\t",index=False)
    (base/("control_relative.audit.json" if relative else "unmixing.audit.json")).write_text(json.dumps({"mixing":cfg["mixing_assumption"],
        "negative_bins":int(sum(bool(r["negative"]) for r in allrows)),"method":cfg["method"],
        "scale":"own contig mean then per-bin reference-matched race mean" if relative else "own full-contig mean; requires equal underlying contig dosage","low_control_mask_threshold":.25 if relative else None,"uncertainty":"mixing sensitivity only; no statistical confidence interval",
        "pdf":str(out)},indent=2)+"\n")
    print("PDF\t"+str(out))
def alleles(cfg,task):
    base=Path(cfg["run_dir"]);ref=list(ANCHORS)[task];contig,lo,hi=ANCHORS[ref]
    focal=Path(cfg["focal_run"]);fc=json.loads((focal/"run.json").read_text())
    tools=fc["tools"];directory=base/"allele_evidence"/ref;directory.mkdir(parents=True)
    templates=directory/"reference_span.fasta"
    fasta=focal/ref/"reference/genome.fa"
    with templates.open("w") as h:
        run([tools["samtools"],"faidx",fasta,f"{contig}:{lo-2000+1}-{hi+2000}"],stdout=h)
    # Align paper Sanger sequences to the complete genome, retaining alternative mappings.
    with (directory/"paper_sanger_to_genome.sam").open("w") as h, (directory/"sanger_mapping.log").open("w") as err:
        run([tools["bwa"],"mem","-a","-t","4",fasta,base/"paper_sanger.fasta"],stdout=h,stderr=err)
    manifest=[]
    for e in cfg["tasks"]:
        if e["reference"]["name"]!=ref:continue
        if e["role"]!="focal" and not (e.get("species")=="melpomene" and e.get("subspecies") in ("rosina","melpomene")):continue
        raw=directory/(e["sample"]+".region_and_mates.unsorted.bam")
        out=directory/(e["sample"]+".region_and_mates.bam")
        run([tools["samtools"],"view","--fetch-pairs","-b","-o",raw,e["bam"],f"{contig}:{lo-2000+1}-{hi+2000}"])
        run([tools["samtools"],"sort","-@", "2","-o",out,raw]);run([tools["samtools"],"index",out])
        raw.unlink()
        manifest.append({"sample":e["sample"],"role":e["role"],"bam":str(out),"source_bam":e["bam"]})
    (directory/"evidence_manifest.json").write_text(json.dumps({
        "reference":ref,"primer_span0":[lo,hi],"flanks":2000,"samples":manifest,
        "read_assignment_performed":False,"genotypes_inferred":False,
        "next":"Align reference spans and paper Sanger sequences; validate unique allele junctions before classifying reads"},indent=2)+"\n")
    print("ALLELE_PREPARATION_COMPLETE\t"+str(directory))
def main():
    p=argparse.ArgumentParser();p.add_argument("--source-run",type=Path);p.add_argument("--config",type=Path)
    p.add_argument("--stage",choices=["mean","alleles","report"]);p.add_argument("--task",type=int)
    a=p.parse_args()
    if a.source_run:submit(a.source_run)
    else:
        cfg=json.loads(a.config.read_text())
        if a.stage=="mean":mean(cfg,a.task)
        elif a.stage=="alleles":alleles(cfg,a.task)
        else:
            from matplotlib.backends.backend_pdf import PdfPages
            out=Path(cfg["run_dir"])/"Yb_d_mean_coverage_and_reference_matched_control_ratios.pdf"
            with PdfPages(out) as bundle:
                report(cfg,bundle=bundle)
                report(cfg,relative=True,bundle=bundle)
            print("BUNDLE_PDF\\t"+str(out))
if __name__=="__main__":
    try:main()
    except Exception as e:print("ERROR: "+str(e),file=sys.stderr);sys.exit(2)
