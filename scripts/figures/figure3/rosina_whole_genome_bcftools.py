#!/usr/bin/env python3
"""Genome-wide ROSINA bcftools callset shared by Optix and Yb.

Uses the established bcftools mpileup/call pipeline, genomic window Slurm
arrays, and the existing BWA/Picard alignments. Calls every contig in the
Rosina reference. NEVER invokes Hmel2.5, GATK or new BAM alignments.
Focal samples remain four separate BAMs, not naive diploid GT samples.
"""
import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

BASE=Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
DEST_BASE=BASE/"analysis/rosina_genomewide_joint_recall"
ALIGN_ARRAY="138468792"
EXCLUDED="SRR11243351"
YB="NC_140403.1"
OPTIX="NC_140406.1"
YB_RANGE=(1600000,1700000)
OPTIX_RANGE=(1060000,1240000)
WINDOW_BP=2_000_000
PADDING=200
SLURM_BATCH=450


def read_json(path):
    return json.loads(Path(path).read_text())


def put_json(path,obj):
    Path(path).write_text(json.dumps(obj,sort_keys=True,indent=2)+"\n")


def read_tsv(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def write_tsv(path,rows,fields):
    with Path(path).open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=fields,extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def exists(path):
    p=Path(path)
    if not p.is_file() or p.stat().st_size==0:
        raise RuntimeError("Missing or empty input: "+str(p))
    return p


def exec_cmd(cmd,stdout=None):
    return subprocess.run([str(x) for x in cmd],check=True,stdout=stdout)


def original_samples(source):
    data=read_tsv(exists(source/"manifest.tsv"))
    prior=read_json(source/"run.json")
    if len(data)!=48 or len(prior["comparative"])!=44 or len(prior["focal"])!=4:
        raise RuntimeError("Unexpected frozen original cohort")
    if [r["sample"] for r in data]!=[r["sample"] for r in prior["samples"]]:
        raise RuntimeError("Manifest order differs from worker array order")
    removed=[r for r in data if r["sample"]==EXCLUDED]
    if len(removed)!=1 or removed[0]["label"]!="P4003":
        raise RuntimeError("Known P4003 exclusion does not match")
    out=[dict(row,array_task=i) for i,row in enumerate(data)
         if row["sample"]!=EXCLUDED]
    comp=[r for r in out if r["role"]=="comparative"]
    focal=[r for r in out if r["role"]=="focal"]
    if len(out)!=47 or len(comp)!=43 or len(focal)!=4:
        raise RuntimeError("Cohort after P4003 exclusion is not 43+4")
    if len({r["sample"] for r in out})!=47 or len({r["specimen"] for r in comp})!=43:
        raise RuntimeError("Duplicated sample or biological specimen")
    groups=Counter(r["group"] for r in comp)
    if groups["cythera"]!=5 or groups["vulcanus"]!=5:
        raise RuntimeError("Five cythera/five vulcanus samples not retained")
    if Counter(r["side"] for r in focal)!={"Left":2,"Right":2}:
        raise RuntimeError("Four focal tissues no longer form L/R pairs")
    if not all(re.fullmatch(r"[A-Za-z0-9_.-]+",r["sample"]) for r in out):
        raise RuntimeError("Unsafe BAM sample ID")
    return out,comp,focal


def reference_contigs(path):
    refs=[]
    for line in exists(str(path)+".fai").read_text().splitlines():
        terms=line.split("\t")
        if len(terms)<2 or not terms[1].isdigit():
            raise RuntimeError("Malformed Rosina .fai")
        if not terms[0] or int(terms[1])<1:
            raise RuntimeError("Empty FASTA contig")
        refs.append(dict(chrom=terms[0],length=int(terms[1])))
    if len(refs)<10 or len({x["chrom"] for x in refs})!=len(refs):
        raise RuntimeError("Invalid whole-genome FASTA index")
    lengths={r["chrom"]:r["length"] for r in refs}
    if lengths.get(YB,0)<YB_RANGE[1] or lengths.get(OPTIX,0)<OPTIX_RANGE[1]:
        raise RuntimeError("Rosina assembly lacks required native Yb/Optix contigs")
    if any("Hmel2.5" in r["chrom"] for r in refs):
        raise RuntimeError("Hmel2.5 reference used instead of Rosina")
    return refs


def windows_for(refs,width=WINDOW_BP):
    if width<1:raise ValueError("Window size must be positive")
    out=[]
    for i,r in enumerate(refs):
        for start in range(1,r["length"]+1,width):
            out.append(dict(task=len(out),reference_index=i,chrom=r["chrom"],
                            start=start,end=min(r["length"],start+width-1),
                            contig_length=r["length"]))
    return out


def verified_prior_bams(path,samples):
    rows=read_tsv(exists(path))
    if not rows or "long_contigs_mapped" not in rows[0]:
        raise RuntimeError("Rerun the updated Rosina whole-genome BAM reuse audit")
    wanted={r["sample"] for r in samples}
    adopted={}
    for row in rows:
        if row["sample"] not in wanted:continue
        good=(row["status"]=="POTENTIAL_REUSE" and
              row["source_ref_explicit_match"]=="YES" and
              int(row["long_contigs_mapped"] or 0)>=10 and
              (row.get("picard_metrics_present")=="YES"
               or row.get("has_picard_pg")=="YES"))
        if not good:continue
        proposed=exists(row["candidate"]).resolve(strict=True)
        if (row["sample"] in adopted and
            Path(adopted[row["sample"]]).resolve()!=proposed):
            raise RuntimeError("Ambiguous reuse BAM for "+row["sample"])
        adopted[row["sample"]]=str(proposed)
    return adopted


def build_plan(source,audit,width):
    source=source.resolve(strict=True)
    rows,comp,focal=original_samples(source)
    prior=read_json(source/"run.json")
    if read_json(source/"jobs.json")["align"]!=ALIGN_ARRAY:
        raise RuntimeError("Frozen alignment array mismatch")
    fasta=exists(prior["original_fasta"]).resolve(strict=True)
    refs=reference_contigs(fasta)
    tasks=windows_for(refs,width)
    if len(tasks)>500000:raise RuntimeError("Window size too small for assembly")
    prior_bams=verified_prior_bams(audit,rows)
    for r in rows:
        source_bam=(Path(prior_bams[r["sample"]])
                    if r["sample"] in prior_bams else
                    source/"bam"/(r["sample"]+".dedup.bam"))
        r["input_bam"]=str(source_bam)
        r["origin"]="reused_old_rosina" if r["sample"] in prior_bams else "new_oct14"
    needed=[int(r["array_task"]) for r in rows if r["origin"]=="new_oct14"]
    config=dict(source=str(source),fasta=str(fasta),
        samples=rows,comparative=[r["sample"] for r in comp],
        focal=[r["sample"] for r in focal],refs=refs,windows=tasks,
        window_bp=width,padding=PADDING,
        verified_reuse=sorted(prior_bams),pending_alignment_ids=needed,
        old_alignment_job=ALIGN_ARRAY,
        audit_path=str(Path(audit).resolve()),
        original_manifest_sha256=hashlib.sha256(
            (source/"manifest.tsv").read_bytes()).hexdigest(),
        reference_fai_sha256=hashlib.sha256(
            Path(str(fasta)+".fai").read_bytes()).hexdigest(),
        standard_call="bcftools mpileup -Ou -a FORMAT/DP,FORMAT/AD -q30 -Q20"
                      " | bcftools call -m -A -v -f GQ",
        variant_only_flag="bcftools call -v changes site emission, not genotype model",
        default_ploidy_note="Sex-linked hemizygous regions need sex-aware downstream validation",
        no_hmel25=True)
    return config


def print_plan(cfg):
    print("REFERENCE\tROSINA WHOLE GENOME; NO HMEL2.5")
    print("CONTIGS\t"+str(len(cfg["refs"])))
    print("REFERENCE_BASES\t"+str(sum(c["length"] for c in cfg["refs"])))
    print("WINDOWS\t"+str(len(cfg["windows"])))
    print("WINDOW_SIZE_BP\t"+str(cfg["window_bp"]))
    print("COMPARATIVE\t43; 5 cythera; 5 vulcanus")
    print("FOCAL\t4; T2+T3 left/right maintained individually")
    print("REUSABLE_OLD_ROSINA_BAMS\t"+str(len(cfg["verified_reuse"])))
    print("ORIGINAL_BAM_ALIGNMENTS_STILL_NEEDED\t"+str(len(cfg["pending_alignment_ids"])))
    print("SNP_ENGINE\tbcftools mpileup | bcftools call")
    print("SHARED_OUTPUT\t"+str(DEST_BASE))
    print("YB_HANDOFF\tChromosome VCF, Yb interval VCF, BAM manifest, focal depths")


def find_bcftools(source):
    history=BASE/"analysis/pk_dark_pale_optix_yb/20261013T220116Z_0b6e0845/run.json"
    if history.is_file():
        bcftools=read_json(history).get("tools",{}).get("bcftools")
        if bcftools and Path(bcftools).is_file():
            return str(Path(bcftools).resolve())
    candidate=source/"toolenv/bin/bcftools"
    return str(exists(candidate).resolve())


def prepare(cfg):
    DEST_BASE.mkdir(parents=True,exist_ok=True)
    output=DEST_BASE/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
                    +"_"+cfg["original_manifest_sha256"][:8])
    if output.exists():raise RuntimeError("No run overwrite: "+str(output))
    output.mkdir()
    for d in ("logs","bams","windows","chromosomes","joint","qc",
              "integration_Yb","integration_Optix"):
        (output/d).mkdir()
    script=output/"rosina_whole_genome_bcftools.py"
    shutil.copy2(__file__,script)
    cfg["output"]=str(output)
    cfg["bcftools"]=find_bcftools(Path(cfg["source"]))
    cfg["frozen_script_sha256"]=hashlib.sha256(script.read_bytes()).hexdigest()
    for row in cfg["samples"]:
        name=row["sample"]
        upstream=Path(row["input_bam"])
        bam=output/"bams"/(name+".dedup.bam")
        bam.symlink_to(upstream)
        Path(str(bam)+".bai").symlink_to(Path(str(upstream)+".bai"))
        row["bam"]=str(bam)
    write_tsv(output/"bam_manifest.tsv",cfg["samples"],list(cfg["samples"][0]))
    (output/"comparative_bams.list").write_text(
        "".join(next(r["bam"] for r in cfg["samples"] if r["sample"]==sample)+"\n"
                for sample in cfg["comparative"]))
    write_tsv(output/"genomic_windows.tsv",cfg["windows"],
              ["task","reference_index","chrom","start","end","contig_length"])
    write_tsv(output/"reference_contigs.tsv",cfg["refs"],["chrom","length"])
    put_json(output/"run.json",cfg)
    worker=output/"worker.sh"
    worker.write_text("#!/usr/bin/env bash\nset -euo pipefail\nexec "+
        shlex.join([sys.executable,str(script),"--run",str(output),"--stage"])+
        ' "$1" --offset "$2" --task "$SLURM_ARRAY_TASK_ID"\n')
    worker.chmod(0o755)
    return output


def submit(cfg):
    output=prepare(cfg)
    submitted={}
    jobsfile=output/"jobs.json"
    def launch(name,total,cpus,mem,limit,dep):
        new=[]
        for offset in range(0,total,SLURM_BATCH):
            size=min(SLURM_BATCH,total-offset)
            args=["sbatch","--parsable","--partition=cpu",
                  "--array=0-"+str(size-1),
                  "--cpus-per-task="+str(cpus),"--mem="+mem,
                  "--time="+limit,"--job-name=RosJoint_"+name,
                  "--output="+str(output/"logs"/(name+".%A_%a.out")),
                  "--error="+str(output/"logs"/(name+".%A_%a.err"))]
            if dep:args+=["--dependency=afterok:"+":".join(dep)]
            args.extend([str(output/"worker.sh"),name,str(offset)])
            jid=subprocess.check_output(args,text=True).strip().split(";")[0]
            new.append(jid)
            submitted[name]=new
            put_json(jobsfile,submitted)
            print("JOB\t"+name+"\t"+jid+"\t"+str(offset)+"-"+
                  str(offset+size-1),flush=True)
        return new
    try:
        deps=[ALIGN_ARRAY+"_"+str(i) for i in cfg["pending_alignment_ids"]]
        gate=launch("validate",1,2,"12G","03:00:00",deps)
        call=launch("call",len(cfg["windows"]),2,"14G","12:00:00",gate)
        ctgs=launch("chrom",len(cfg["refs"]),2,"14G","06:00:00",call)
        full=launch("combine",1,4,"24G","12:00:00",ctgs)
        derived=launch("extract",2,2,"12G","04:00:00",full)
        foc=launch("focal",2,2,"16G","08:00:00",derived)
        launch("publish",1,4,"20G","08:00:00",foc)
    except Exception:
        for ids in submitted.values():
            for jid in ids:
                subprocess.run(["scancel",jid],stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)
        raise
    print("ORIGINAL_BAM_ARRAY_RETAINED\t"+ALIGN_ARRAY)
    print("HMEDL2.5_CALL_JOBS\t0")
    print("GATK_CALL_JOBS\t0")
    print("RUN_DIR\t"+str(output))


def validate(cfg,index):
    import pysam
    if index!=0:raise RuntimeError("Bad validation task")
    output=Path(cfg["output"])
    if hashlib.sha256(Path(cfg["fasta"]+".fai").read_bytes()).hexdigest()!=(
            cfg["reference_fai_sha256"]):
        raise RuntimeError("Reference index drift since submission")
    expected=[(r["chrom"],r["length"]) for r in cfg["refs"]]
    verified=[]
    for row in cfg["samples"]:
        sample=row["sample"]
        bam=exists(row["bam"])
        exists(row["bam"]+".bai")
        with pysam.AlignmentFile(str(bam),"rb") as f:
            head=f.header.to_dict()
            sq=[(e["SN"],int(e["LN"])) for e in head.get("SQ",[])]
            rg={e.get("SM") for e in head.get("RG",[])}
            if sq!=expected or rg!={sample}:
                raise RuntimeError("BAM SQ/RG mismatch: "+sample)
            coverage={e.contig:e.mapped for e in f.get_index_statistics()}
            if sum(coverage.get(name,0)>=20 for name,length in expected
                   if length>=100000)<10:
                raise RuntimeError("Not genome-wide mapped BAM: "+sample)
        verified.append(dict(sample=sample,bam=str(bam),
                             reused=row["origin"]=="reused_old_rosina"))
    write_tsv(output/"qc"/"verified_bams.tsv",verified,list(verified[0]))
    version=subprocess.check_output([cfg["bcftools"],"--version"],
                                    text=True).splitlines()[0]
    put_json(output/"qc"/"validated.json",
             dict(sample_count=47,comparative_count=43,reference=cfg["fasta"],
                  bcftools_version=version))
    print("ROSINA_47_BAMS_VALIDATED\t"+version,flush=True)


def call_window(cfg,index):
    import pysam
    output=Path(cfg["output"])
    exists(output/"qc"/"validated.json")
    x=cfg["windows"][index]
    c,s,e=x["chrom"],x["start"],x["end"]
    stem=str(index).zfill(7)
    dest=output/"windows"/(stem+".bcf")
    done=output/"windows"/(stem+".complete.json")
    if dest.is_file() and done.is_file() and Path(str(dest)+".csi").is_file():
        if read_json(done)["region"]!=[c,s,e]:
            raise RuntimeError("Cached window region mismatch")
        return
    tmp=output/"windows"/(stem+".partial.bcf")
    tmp.unlink(missing_ok=True)
    Path(str(tmp)+".csi").unlink(missing_ok=True)
    padded=c+":"+str(max(1,s+(-cfg["padding"])))+"-"+str(
        min(x["contig_length"],e+cfg["padding"]))
    pile=[cfg["bcftools"],"mpileup","-Ou","-f",cfg["fasta"],
          "-r",padded,"-a","FORMAT/DP,FORMAT/AD","-q","30","-Q","20",
          "-b",str(output/"comparative_bams.list")]
    # Previous model retained; -v changes only variant site emission.
    call=[cfg["bcftools"],"call","-m","-A","-v","-f","GQ","-Ou"]
    crop=[cfg["bcftools"],"view","-i",
          "POS>="+str(s)+" && POS<="+str(e),"-Ob","-o",str(tmp)]
    pipeline=" | ".join(shlex.join(p) for p in [pile,call,crop])
    subprocess.run(["bash","-o","pipefail","-c",pipeline],check=True)
    exists(tmp)
    exec_cmd([cfg["bcftools"],"index","-c",str(tmp)])
    exists(str(tmp)+".csi")
    n=0
    with pysam.VariantFile(str(tmp)) as vf:
        if set(vf.header.samples)!=set(cfg["comparative"]) or len(vf.header.samples)!=43:
            raise RuntimeError("Wrong joint cohort in window "+str(index))
        last=s-1
        for v in vf:
            if v.contig!=c or v.pos<last or not s<=v.pos<=e:
                raise RuntimeError("Unsorted or off-window variant: "+str(index))
            last=v.pos
            n+=1
    tmp.replace(dest)
    Path(str(tmp)+".csi").replace(str(dest)+".csi")
    put_json(done,dict(region=[c,s,e],padded=padded,records=n,
                       comparative=43))
    print("CALLED_WINDOW\t"+str(index)+"\t"+str(n),flush=True)


def concat_bcfs(bcftools,inputs,output):
    for item in inputs:
        exists(item)
        exists(str(item)+".csi")
    plist=Path(str(output)+".list")
    plist.write_text("".join(str(p)+"\n" for p in inputs))
    output.unlink(missing_ok=True)
    Path(str(output)+".csi").unlink(missing_ok=True)
    exec_cmd([bcftools,"concat","-f",str(plist),
              "-Ob" if output.suffix==".bcf" else "-Oz",
              "-o",str(output)])
    exec_cmd([bcftools,"index","-c",str(output)])


def combine_contig(cfg,index):
    output=Path(cfg["output"])
    r=cfg["refs"][index]
    target=r["chrom"]
    items=[x for x in cfg["windows"] if x["chrom"]==target]
    if not items or items[0]["start"]!=1 or items[-1]["end"]!=r["length"]:
        raise RuntimeError("Contig window coverage incomplete "+target)
    if any(a["end"]+1!=b["start"] for a,b in zip(items,items[1:])):
        raise RuntimeError("Contig intervals have a gap "+target)
    files=[]
    total=0
    for w in items:
        name=str(w["task"]).zfill(7)
        marker=read_json(exists(output/"windows"/(name+".complete.json")))
        if marker["region"]!=[target,w["start"],w["end"]]:
            raise RuntimeError("Mismatched window metadata "+name)
        total+=marker["records"]
        files.append(output/"windows"/(name+".bcf"))
    name=str(index).zfill(5)
    dest=output/"chromosomes"/(name+".bcf")
    done=output/"chromosomes"/(name+".complete.json")
    if dest.is_file() and done.is_file():return
    tmp=output/"chromosomes"/(name+".partial.bcf")
    concat_bcfs(cfg["bcftools"],files,tmp)
    tmp.replace(dest)
    Path(str(tmp)+".csi").replace(str(dest)+".csi")
    put_json(done,dict(chrom=target,length=r["length"],records=total,
                       interval_count=len(items)))
    print("CONTIG_COMPLETE\t"+target+"\t"+str(total),flush=True)


def combine_genome(cfg,index):
    import pysam
    if index!=0:raise RuntimeError("Bad combine task")
    output=Path(cfg["output"])
    files=[]
    n=0
    for i,r in enumerate(cfg["refs"]):
        stem=str(i).zfill(5)
        marker=read_json(exists(output/"chromosomes"/(stem+".complete.json")))
        if marker["chrom"]!=r["chrom"]:
            raise RuntimeError("Reference contig order changed")
        n+=marker["records"]
        files.append(output/"chromosomes"/(stem+".bcf"))
    dest=output/"joint"/"Rosina.genomewide.43comparative.variants.vcf.gz"
    complete=output/"joint"/"combined.complete.json"
    if dest.is_file() and complete.is_file():return
    tmp=output/"joint"/"Rosina.genomewide.partial.vcf.gz"
    concat_bcfs(cfg["bcftools"],files,tmp)
    with pysam.VariantFile(str(tmp)) as vf:
        if set(vf.header.samples)!=set(cfg["comparative"]) or len(vf.header.samples)!=43:
            raise RuntimeError("Final joint VCF sample set mismatch")
    tmp.replace(dest)
    Path(str(tmp)+".csi").replace(str(dest)+".csi")
    snps=output/"joint"/"Rosina.genomewide.43comparative.snps.vcf.gz"
    exec_cmd([cfg["bcftools"],"view","-v","snps","-Oz",
              "-o",str(snps),str(dest)])
    exec_cmd([cfg["bcftools"],"index","-c",str(snps)])
    put_json(complete,dict(reference="Rosina",contigs=len(cfg["refs"]),
                           windows=len(cfg["windows"]),variants=n,
                           cohort=43,vcf=str(dest),snps=str(snps)))
    print("WHOLE_GENOME_VCFS_READY\t"+str(dest),flush=True)


def export_locus(cfg,index):
    if index not in (0,1):raise RuntimeError("Bad export task")
    output=Path(cfg["output"])
    full=read_json(exists(output/"joint"/"combined.complete.json"))
    label,chrom,limits=(("Yb",YB,YB_RANGE) if index==0
                        else ("Optix",OPTIX,OPTIX_RANGE))
    folder=output/("integration_"+label)
    for typ,region in [("chromosome",chrom),
                       ("window",chrom+":"+str(limits[0])+"-"+str(limits[1]))]:
        result=folder/(label+"_"+typ+".vcf.gz")
        exec_cmd([cfg["bcftools"],"view","-r",region,"-Oz",
                  "-o",str(result),full["vcf"]])
        exec_cmd([cfg["bcftools"],"index","-c",str(result)])
    bam_manifest=folder/"rosina_bam_manifest.tsv"
    write_tsv(bam_manifest,cfg["samples"],list(cfg["samples"][0]))
    put_json(folder/"complete.json",dict(label=label,reference="Rosina",
             chrom=chrom,window=list(limits),genomewide_variants=full["vcf"],
             genomewide_snps=full["snps"],bam_manifest=str(bam_manifest),
             chrom_vcf=str(folder/(label+"_chromosome.vcf.gz")),
             window_vcf=str(folder/(label+"_window.vcf.gz")),
             no_hmel25_snp_calls=True))
    print("SHARED_HANDOFF_READY\t"+label+"\t"+str(bam_manifest),flush=True)


def focal_counts(cfg,index):
    import gzip
    import pysam
    if index not in (0,1):raise RuntimeError("Bad focal task")
    output=Path(cfg["output"])
    label="Yb" if index==0 else "Optix"
    folder=output/("integration_"+label)
    locus=read_json(exists(folder/"complete.json"))
    ch=locus["chrom"]
    first,last=locus["window"]
    sites={}
    with pysam.VariantFile(str(locus["window_vcf"])) as vf:
        for rec in vf.fetch(ch,first-1,last):
            if (len(rec.ref)==1 and rec.ref.upper() in "ACGT" and rec.alts
                    and len(rec.alts)==1 and len(rec.alts[0])==1
                    and rec.alts[0].upper() in "ACGT"):
                sites[rec.pos]=(rec.ref.upper(),rec.alts[0].upper())
    if not sites:raise RuntimeError("No biallelic SNPs in "+label+" window")
    samples={r["sample"]:r for r in cfg["samples"]}
    counts={}
    for name in cfg["focal"]:
        found=defaultdict(lambda:[0,0])
        with pysam.AlignmentFile(str(exists(samples[name]["bam"])),"rb") as f:
            for col in f.pileup(ch,first-1,last,truncate=True,
                    stepper="samtools",min_base_quality=20,
                    min_mapping_quality=30,ignore_overlaps=True,
                    flag_filter=1796,max_depth=100000):
                pos=col.reference_pos+1
                if pos not in sites:continue
                ref,alt=sites[pos]
                for aln in col.pileups:
                    if aln.is_del or aln.is_refskip or aln.query_position is None:
                        continue
                    base=aln.alignment.query_sequence[aln.query_position].upper()
                    if base==ref:found[pos][0]+=1
                    elif base==alt:found[pos][1]+=1
        counts[name]=found
    dest=folder/(label+"_focal_Left_Right_allele_counts.tsv.gz")
    with gzip.open(dest,"wt",newline="") as h:
        writer=csv.writer(h,delimiter="\t")
        writer.writerow(["chrom","pos","ref","alt","L_ref","L_alt",
                         "R_ref","R_alt","L_AF","R_AF","per_library_json"])
        for pos,(ref,alt) in sorted(sites.items()):
            side={"Left":[0,0],"Right":[0,0]}
            per={}
            for name in cfg["focal"]:
                q=counts[name].get(pos,[0,0])
                per[name]=q
                dest_side=side[samples[name]["side"]]
                dest_side[0]+=q[0]
                dest_side[1]+=q[1]
            L,R=side["Left"],side["Right"]
            af=lambda z:z[1]/sum(z) if sum(z) else ""
            writer.writerow([ch,pos,ref,alt,*L,*R,af(L),af(R),
                             json.dumps(per,sort_keys=True)])
    put_json(folder/"focal.complete.json",
             dict(biallelic_snps=len(sites),files=str(dest),
                  focal_libs=cfg["focal"],per_tissue_counts_retained=True,
                  pooled_T2_T3_Left_Right=True,min_mapq=30,min_bq=20))
    print("FOCAL_COUNTS_READY\t"+label+"\t"+str(len(sites)),flush=True)


def publish(cfg,index):
    if index!=0:raise RuntimeError("Bad publish task")
    output=Path(cfg["output"])
    full=read_json(exists(output/"joint"/"combined.complete.json"))
    for locus in ("Yb","Optix"):
        where=output/("integration_"+locus)
        exists(where/"complete.json")
        exists(where/"focal.complete.json")
    stats=output/"qc"/"joint_bcftools_stats.txt"
    with stats.open("w") as handle:
        exec_cmd([cfg["bcftools"],"stats","-s","-",full["vcf"]],
                 stdout=handle)
    put_json(output/"complete.json",dict(
        status="COMPLETE",reference="Rosina",no_Hmel25=True,
        comparative=43,focal=4,contigs=len(cfg["refs"]),
        calling_windows=len(cfg["windows"]),reuse=len(cfg["verified_reuse"]),
        genomewide_vcf=full["vcf"],genomewide_snps=full["snps"],
        yb=read_json(output/"integration_Yb/complete.json"),
        optix=read_json(output/"integration_Optix/complete.json"),
        original_alignment_job=ALIGN_ARRAY,
        default_ploidy_note=cfg["default_ploidy_note"]))
    (DEST_BASE/"latest_complete.path").write_text(str(output)+"\n")
    print("ROSINA_GENOMEWIDE_COMPLETE\t"+str(output),flush=True)


STAGES={"validate":validate,"call":call_window,"chrom":combine_contig,
        "combine":combine_genome,"extract":export_locus,
        "focal":focal_counts,"publish":publish}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run",type=Path)
    parser.add_argument("--reuse-audit",type=Path)
    parser.add_argument("--window-bp",type=int,default=WINDOW_BP)
    parser.add_argument("--preview",action="store_true")
    parser.add_argument("--submit",action="store_true")
    parser.add_argument("--run",type=Path)
    parser.add_argument("--stage",choices=STAGES)
    parser.add_argument("--offset",type=int,default=0)
    parser.add_argument("--task",type=int)
    a=parser.parse_args()
    if a.preview or a.submit:
        if a.source_run is None or a.reuse_audit is None:
            parser.error("--source-run and --reuse-audit required")
        config=build_plan(a.source_run,a.reuse_audit,a.window_bp)
        print_plan(config)
        if a.submit:submit(config)
    elif a.run and a.stage and a.task is not None:
        config=read_json(a.run/"run.json")
        STAGES[a.stage](config,a.offset+a.task)
    else:
        parser.error("--preview/--submit or --run --stage --task required")


if __name__=="__main__":
    try:
        main()
    except Exception as exc:
        print("ROSINA_WHOLE_GENOME_ERROR\t"+str(exc),
              file=sys.stderr,flush=True)
        sys.exit(2)
