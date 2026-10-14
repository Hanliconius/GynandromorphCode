#!/usr/bin/env python3
"""Resume the existing Optix/Yb Rosina BAMs with the established bcftools caller.

No GATK and no gVCFs. Reuse all BWA/Picard outputs except known-failed
P4003 / SRR11243351. The calls come from joint multi-BAM mpileup followed
by bcftools call, using the exact previous MAPQ/BQ/calling options.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shlex
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

EXCLUDED="SRR11243351"
ALIGN_ARRAY="138468792"
LOCI={"Optix":("NC_140406.1",1060000,1240000),
      "Yb":("NC_140403.1",1600000,1700000)}
SHARD_BP=5000


def make_shards(loci=LOCI, width=SHARD_BP):
    """Partition 1-based inclusive intervals with no overlaps or gaps.

    Absorb a last <=width remainder into the preceding shard to avoid tiny
    one-base jobs. Never split the *samples* between joint calls.
    """
    if width < 1:
        raise ValueError("Shard width must be positive")
    result=[]
    for locus,(contig,start,end) in loci.items():
        if start < 1 or end < start:
            raise ValueError("Invalid locus "+locus)
        pos=start
        while pos<=end:
            stop=min(end,pos+width-1)
            if end-stop<=width:
                stop=end
            result.append(dict(locus=locus,chrom=contig,start=pos,end=stop))
            pos=stop+1
    return result
MP_OPTIONS=["-Ou","-a","FORMAT/DP,FORMAT/AD","-q","30","-Q","20"]
CALL_OPTIONS=["-m","-A","-f","GQ","-Oz"]


def load(path):
    return json.loads(Path(path).read_text())


def dump(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True)+"\n")


def manifest(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def require(path):
    path=Path(path)
    if not path.is_file() or path.stat().st_size==0:
        raise RuntimeError("Missing required file: "+str(path))
    return path


def cohort(src):
    original=manifest(src/"manifest.tsv")
    config=load(src/"run.json")
    if len(original)!=48 or len(config["comparative"])!=44:
        raise RuntimeError("Original run expected 44 comparisons + 4 focal")
    if len([r for r in original if r["sample"]==EXCLUDED])!=1:
        raise RuntimeError("Expected exactly one excluded sample")
    dropped=next(r for r in original if r["sample"]==EXCLUDED)
    if (dropped["label"],dropped["role"])!=("P4003","comparative"):
        raise RuntimeError("Excluded sample identity/role unexpected")
    selected=[r for r in original if r["sample"]!=EXCLUDED]
    comparisons=[r for r in selected if r["role"]=="comparative"]
    focal=[r for r in selected if r["role"]=="focal"]
    if (len(comparisons),len(focal))!=(43,4):
        raise RuntimeError("Expected 43 comparisons and four focal libraries")
    if len({r["specimen"] for r in comparisons})!=43:
        raise RuntimeError("Comparative specimen duplication")
    groups=Counter(r["group"] for r in comparisons)
    if (groups["cythera"],groups["vulcanus"])!=(5,5):
        raise RuntimeError("Race cohort not preserved")
    if Counter(r["side"] for r in focal)!={"Left":2,"Right":2}:
        raise RuntimeError("Focal side assignment changed")
    indices=[i for i,r in enumerate(original) if r["sample"]!=EXCLUDED]
    if len(indices)!=47 or sorted(indices)!=[i for i in range(48) if i!=14]:
        raise RuntimeError("Excluded Slurm array index is not 14")
    return selected,comparisons,focal,indices


def prepare(src):
    src=src.resolve()
    selected,comparisons,focal,indices=cohort(src)
    prev=load(src/"jobs.json")
    if prev.get("align")!=ALIGN_ARRAY:
        raise RuntimeError("Original alignment job ID changed")
    # Reuse exactly the project's established bcftools executable if present.
    # This helps preserve its caller version as well as command-line settings.
    historical=Path(
        "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/"
        "analysis/pk_dark_pale_optix_yb/20261013T220116Z_0b6e0845/run.json"
    )
    previous_tool=None
    if historical.is_file():
        previous_tool=load(historical).get("tools",{}).get("bcftools")
    preferred=Path(previous_tool) if previous_tool else src/"toolenv/bin/bcftools"
    bcftools=require(preferred)
    require(src/"ref/complete")
    original=load(src/"run.json")
    require(original["fasta"])
    require(original["fasta"]+".fai")
    out=src/"bcftools_joint_recall"
    cfgpath=out/"run.json"
    if out.exists():
        if not cfgpath.exists() or load(cfgpath)["comparative"]!=[
                r["sample"] for r in comparisons]:
            raise RuntimeError("Existing recovery run conflicts with expected samples")
        return out
    out.mkdir()
    for d in ["logs","joint","shards","focal","audit"]:
        (out/d).mkdir()
    shutil.copy2(__file__,out/"resume_bcftools_optix_joint.py")
    (out/"bam").symlink_to(src/"bam",target_is_directory=True)
    with (out/"comparative_bams.list").open("w") as h:
        for r in comparisons:
            h.write(str(src/"bam"/(r["sample"]+".dedup.bam"))+"\n")
    with (out/"manifest.tsv").open("w",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(selected[0]),delimiter="\t",extrasaction="ignore")
        w.writeheader()
        w.writerows(selected)
    planned=make_shards()
    dump(out/"shards.json",planned)
    cfg=dict(root=str(out),source_run=str(src),fasta=original["fasta"],
        bcftools=str(bcftools),bam_list=str(out/"comparative_bams.list"),
        comparative=[r["sample"] for r in comparisons],
        focal=[r["sample"] for r in focal],
        samples=selected,loci=LOCI,shards=planned,
        shard_bp=SHARD_BP,selected_alignment_indices=indices,
        excluded_sample=EXCLUDED,excluded_task=14,
        expected_method="bcftools mpileup -a FORMAT/DP,FORMAT/AD -q30 -Q20 | bcftools call -m -A -f GQ",
        previous_manifest_sha256=hashlib.sha256((src/"manifest.tsv").read_bytes()).hexdigest(),
        new_script_sha256=hashlib.sha256(
            (out/"resume_bcftools_optix_joint.py").read_bytes()).hexdigest())
    dump(cfgpath,cfg)
    worker=out/"worker.sh"
    worker.write_text("#!/usr/bin/env bash\nset -euo pipefail\nexec "+
        shlex.join([sys.executable,str(out/"resume_bcftools_optix_joint.py"),
                    "--config",str(cfgpath),"--stage"])+
        ' "$1" --task "$SLURM_ARRAY_TASK_ID"\n')
    worker.chmod(0o755)
    return out


def submit(src):
    out=prepare(src)
    cfg=load(out/"run.json")
    submitted={}
    jobs_path=out/"jobs.json"
    if jobs_path.exists():
        raise RuntimeError("This recovery run has already submitted jobs")
    # Specific array elements are an established Slurm dependency syntax.
    # Do not depend on whole 138468792: element 14 is known to fail.
    dep="afterok:"+":".join(ALIGN_ARRAY+"_"+str(i) for i in
                           cfg["selected_alignment_indices"])
    def sb(stage,task_count,cpus,mem,time,depends):
        cmd=["sbatch","--parsable","--partition=cpu",
             "--array=0-"+str(task_count-1),"--cpus-per-task="+str(cpus),
             "--mem="+mem,"--time="+time,"--dependency="+depends,
             "--job-name=OxBCF_"+stage,
             "--output="+str(out/"logs"/(stage+".%A_%a.out")),
             "--error="+str(out/"logs"/(stage+".%A_%a.err")),
             str(out/"worker.sh"),stage]
        jid=subprocess.check_output(cmd,text=True).strip().split(";")[0]
        submitted[stage]=jid
        dump(jobs_path,submitted)
        print("JOB\t"+stage+"\t"+jid,flush=True)
        return jid
    try:
        shard=sb("shard",len(cfg["shards"]),1,"8G","08:00:00",dep)
        combined=sb("concat",2,2,"12G","03:00:00","afterok:"+shard)
        focaljob=sb("focal",2,2,"16G","06:00:00","aftercorr:"+combined)
        sb("audit",1,2,"12G","02:00:00","afterok:"+focaljob)
    except Exception:
        for jid in submitted.values():
            subprocess.run(["scancel",jid],stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
        raise
    print("METHOD\tbcftools mpileup -> bcftools call (established parameters)")
    print("ARRAYS\t"+str(len(cfg["shards"]))+
          " parallel genomic shards; 2 concat; 2 focal; 1 audit")
    print("COHORT\t43 comparative, 4 focal, excluded P4003")
    print("RUN_DIR\t"+str(out))
    print("GATK_SUBMISSIONS\t0")



def check_bams(cfg):
    import pysam
    old=Path(cfg["source_run"])
    for name in cfg["comparative"]+cfg["focal"]:
        bam=require(old/"bam"/(name+".dedup.bam"))
        require(str(bam)+".bai")
        meta=load(old/"bam"/(name+".complete.json"))
        if meta.get("sample")!=name or not meta.get("full_genome_alignment"):
            raise RuntimeError("Bad BAM provenance "+name)
        with pysam.AlignmentFile(str(bam),"rb") as f:
            tags={rg.get("SM") for rg in f.header.to_dict().get("RG",[])}
            if tags!={name}:
                raise RuntimeError("Incorrect BAM read-group sample "+name)


def joint(cfg,task):
    import pysam
    out=Path(cfg["root"])
    locus=list(cfg["loci"])[task]
    ch,start,end=cfg["loci"][locus]
    dest=out/"joint"/(locus+".all_sites.vcf.gz")
    done=out/"joint"/(locus+".complete.json")
    if dest.exists() and done.exists():
        return
    check_bams(cfg)
    tmp=out/"joint"/(locus+".partial.vcf.gz")
    tmp.unlink(missing_ok=True)
    Path(str(tmp)+".tbi").unlink(missing_ok=True)
    mp=[cfg["bcftools"],"mpileup","-Ou",
        "-f",cfg["fasta"],"-r",f"{ch}:{start}-{end}",
        "-a","FORMAT/DP,FORMAT/AD","-q","30","-Q","20",
        "-b",cfg["bam_list"]]
    call=[cfg["bcftools"],"call","-m","-A","-f","GQ",
          "-Oz","-o",str(tmp)]
    print("CALL\t"+locus+"\t43 comparative BAMs",flush=True)
    p1=subprocess.Popen(mp,stdout=subprocess.PIPE)
    try:
        p2=subprocess.Popen(call,stdin=p1.stdout)
        p1.stdout.close()
        status2=p2.wait()
        status1=p1.wait()
    finally:
        if p1.poll() is None:
            p1.kill()
            p1.wait()
    if status1 or status2:
        raise RuntimeError("bcftools mpileup/call exit statuses: "+
                           str((status1,status2)))
    require(tmp)
    subprocess.run([cfg["bcftools"],"index","-t",str(tmp)],check=True)
    require(str(tmp)+".tbi")
    with pysam.VariantFile(str(tmp)) as f:
        names=list(f.header.samples)
        sites=sum(1 for _ in f)
    if len(names)!=43 or set(names)!=set(cfg["comparative"]) or sites==0:
        raise RuntimeError("Incorrect bcftools VCF sample set or zero sites")
    tmp.replace(dest)
    Path(str(tmp)+".tbi").replace(str(dest)+".tbi")
    version=subprocess.check_output([cfg["bcftools"],"--version"],text=True).splitlines()[0]
    dump(done,dict(caller="bcftools",bcftools_version=version,
                   locus=locus,sites=sites,samples=names,
                   mpileup_options=MP_OPTIONS,call_options=CALL_OPTIONS,
                   source_method="scripts/figures/figure3/extend_PK_optix_Yb.py",
                   legacy_VCF_used=False,gvcf_used=False))
    print("JOINT_COMPLETE\t"+locus+"\t"+str(sites),flush=True)


def get_previous_focal_and_audit(cfg):
    # The previous pipeline used exactly the expected pysam focal pileup/QC.
    # Reuse only those two functions. NEVER execute its GATK stages.
    source=Path(cfg["source_run"])/"optix_joint_recall.py"
    spec=importlib.util.spec_from_file_location("original_focal_counts",source)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.focal,mod.audit


def focal(cfg,task):
    method,_=get_previous_focal_and_audit(cfg)
    method(cfg,task)


def audit(cfg,task):
    _,method=get_previous_focal_and_audit(cfg)
    method(cfg,task)
    # The old report contains a literal 44; this recovery excludes P4003.
    marker=Path(cfg["root"])/"audit.complete.json"
    record=load(marker)
    record["comparative"]=43
    record["caller"]="bcftools mpileup|call"
    record["excluded_comparative"]=EXCLUDED
    record["gvcf_used"]=False
    dump(marker,record)
    print("BCFTOOLS_AUDIT_COMPLETE\t43 comparative; 4 focal",flush=True)


STAGES={"joint":joint,"focal":focal,"audit":audit}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--submit",action="store_true")
    p.add_argument("--preview",action="store_true")
    p.add_argument("--source-run",type=Path)
    p.add_argument("--config",type=Path)
    p.add_argument("--stage",choices=list(STAGES))
    p.add_argument("--task",type=int)
    a=p.parse_args()
    if a.submit or a.preview:
        if a.source_run is None:
            p.error("Missing --source-run")
        _,comparisons,focal,_=cohort(a.source_run)
        print("COHORT_OK\t"+str(len(comparisons))+
              " comparative + "+str(len(focal))+" focal; "+EXCLUDED+" excluded")
        print("CALLER\tbcftools mpileup -q30 -Q20 | bcftools call -m -A -f GQ")
        if a.submit:
            submit(a.source_run)
    elif a.stage and a.config and a.task is not None:
        STAGES[a.stage](load(a.config),a.task)
    else:
        p.error("Specify --preview/--submit or --config/--stage/--task")


if __name__=="__main__":
    try:
        main()
    except Exception as e:
        print("BCFTOOLS_RESUME_ERROR\t"+str(e),file=sys.stderr,flush=True)
        sys.exit(2)
