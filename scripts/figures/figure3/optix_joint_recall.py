#!/usr/bin/env python3
"""Slurm-orchestrated joint Optix/Yb recall: 44 comparative WGS, four focal libraries.

Fresh Rosina whole-genome alignments; diploid comparative GVCFs only.
The gynandromorph is not diploid-genotyped: its four BAMs feed pooled L/R
allele counts at newly discovered biallelic SNPs.
"""
import argparse,csv,gzip,hashlib,json,os,re,shlex,shutil,subprocess,sys,tempfile
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT=Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
LOCI={"Optix":("NC_140406.1",1060000,1240000),
      "Yb":("NC_140403.1",1600000,1700000)}
FIELDS=["sample","specimen","role","group","label","ena_run","ena_sample",
        "r1","r2","source","notes","side","fastq_r1","fastq_r2",
        "md5_r1","md5_r2","bytes_r1","bytes_r2"]


def need(p):
    p=Path(p)
    if not p.is_file() or not p.stat().st_size:
        raise RuntimeError("Missing required file: "+str(p))
    return p


def read(path):
    return json.loads(need(path).read_text())


def write(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True)+"\n")


def table(path):
    with need(path).open() as f:
        return list(csv.DictReader(f,delimiter="\t"))


def tsv(path,rows,fields):
    with Path(path).open("w") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def run(cmd,**kw):
    subprocess.run(list(map(str,cmd)),check=True,**kw)


def md5(path):
    h=hashlib.md5()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b""):
            h.update(chunk)
    return h.hexdigest()


def select_samples(preflight,candidates,sweep):
    original=table(preflight/"manifest.preflight.tsv")
    historical=table(candidates/"verified_unique_wgs_candidates.tsv")
    extra=table(sweep/"additional_cythera_wgs_candidates.tsv")
    audit=table(sweep/"previous_race_taxonomy_audit.tsv")
    if len(original)!=38 or Counter(r["role"] for r in original)!={"comparative":34,"focal":4}:
        raise RuntimeError("Expected exactly 34 comparative + 4 focal inputs")
    if any(x["status"]!="RACE_METADATA_SUPPORTED" for x in audit):
        raise RuntimeError("Historical race audit found unresolved/contradictory metadata")
    cyth=sorted([r for r in historical if r["subspecies"]=="cythera"],
                key=lambda x:x["biosample"])
    if len(cyth)!=2:
        raise RuntimeError("Expected two archival cythera genomes")
    novel={r["biosample"]:r for r in extra if r["subspecies"]=="cythera"
           and r["evidence"]=="TAXON_CONFIRMED" and r["wgs_pair"]=="yes"
           and r["already_in_previous"]=="no"}
    novel=sorted(novel.values(),key=lambda r:(int(r["bytes_r1"])+int(r["bytes_r2"]),
                                               r["biosample"]))
    if len(novel)<3:
        raise RuntimeError("Insufficient new verified cythera genomes")
    newcyth=[novel[round((len(novel)-1)*q)] for q in (0.5,0.7,0.9)]
    vul={r["biosample"]:r for r in historical if r["subspecies"]=="vulcanus"}
    vul=sorted(vul.values(),key=lambda r:(
        0 if r["source_evidence"]=="Martin_2016_archival_seed" else 1,
        -(int(r["bytes_r1"])+int(r["bytes_r2"])),r["biosample"]))
    if len(vul)<5:
        raise RuntimeError("Insufficient verified vulcanus genomes")
    picked=cyth+newcyth+vul[:5]
    if Counter(r["subspecies"] for r in picked)!={"cythera":5,"vulcanus":5}:
        raise RuntimeError("New race balance mismatch")
    result=[]
    for row in original:
        row={**{k:"" for k in FIELDS},**row}
        if row["role"]=="focal":
            s=(row["sample"]+" "+row["label"]).lower()
            row["side"]="Left" if ("left" in s or re.search(r"t[23]l",s)) else (
                "Right" if ("right" in s or re.search(r"t[23]r",s)) else "")
            if not row["side"]:
                raise RuntimeError("No left/right assignment for focal "+row["sample"])
        result.append(row)
    for row in picked:
        result.append({**{k:"" for k in FIELDS},
            "sample":row["run_accession"],"specimen":row["biosample"],
            "role":"comparative","group":row["subspecies"],
            "label":row.get("label") or row["biosample"],
            "ena_run":row["run_accession"],"ena_sample":row["biosample"],
            "source":row.get("source_evidence","ENA confirmed taxon"),
            **{k:row.get(k,"") for k in ("fastq_r1","fastq_r2","md5_r1",
                                         "md5_r2","bytes_r1","bytes_r2")}})
    if len(result)!=48 or len({r["sample"] for r in result})!=48:
        raise RuntimeError("Expected 48 unique libraries")
    cmp=[r for r in result if r["role"]=="comparative"]
    if len(cmp)!=44 or len({r["specimen"] for r in cmp})!=44:
        raise RuntimeError("Expected 44 distinct comparative specimens")
    if Counter(r["side"] for r in result if r["role"]=="focal")!={"Left":2,"Right":2}:
        raise RuntimeError("Expected two libraries on each focal side")
    if any(not re.fullmatch(r"[a-zA-Z0-9_.-]+",r["sample"]) for r in result):
        raise RuntimeError("Unsafe sample name")
    return result


def submit(a):
    pf,cr,sw=(p.resolve() for p in (a.preflight,a.candidates,a.sweep))
    rows=select_samples(pf,cr,sw)
    if a.preview:
        print("COHORT_OK\t44 comparative; 4 focal; 5 cythera; 5 vulcanus")
        for row in rows:
            if row["group"] in ("cythera","vulcanus"):
                print("PICK\t"+row["group"]+"\t"+row["sample"]+"\t"+row["specimen"])
        print("SUBMITTED\t0")
        return
    audited=read(pf/"audit.json")["reference"]
    original_fasta=need(audited["reference"]).resolve()
    import pysam
    with pysam.FastaFile(str(original_fasta)) as fa:
        for contig,beg,end in LOCI.values():
            if contig not in fa.references or fa.get_reference_length(contig)<end:
                raise RuntimeError("Rosina FASTA lacks locus "+contig)
    output=(a.output or ROOT/"analysis/optix_joint_all_samples"/
            datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")).resolve()
    if output.exists():
        raise RuntimeError("Output already exists "+str(output))
    for d in ("logs","raw","bam","gvcf","joint","focal","ref"):
        (output/d).mkdir(parents=True,exist_ok=True)
    shutil.copy2(__file__,output/"optix_joint_recall.py")
    (output/"ref/genome.fa").symlink_to(original_fasta)
    tsv(output/"manifest.tsv",rows,FIELDS)
    cfg=dict(root=str(output),fasta=str(output/"ref/genome.fa"),
             env=str(output/"toolenv"),samples=rows,loci=LOCI,
             comparative=[r["sample"] for r in rows if r["role"]=="comparative"],
             focal=[r["sample"] for r in rows if r["role"]=="focal"],
             source_dirs=dict(preflight=str(pf),candidate_run=str(cr),sweep=str(sw)),
             original_fasta=str(original_fasta),
             original_fasta_size=original_fasta.stat().st_size,
             frozen_ref_allele_audit=audited)
    write(output/"run.json",cfg)
    worker=output/"worker.sh"
    worker.write_text("#!/usr/bin/env bash\nset -euo pipefail\nexec "+
                      shlex.join([sys.executable,str(output/"optix_joint_recall.py"),
                                  "--config",str(output/"run.json"),"--stage"])+
                      ' "$1" --task "$SLURM_ARRAY_TASK_ID"\n')
    worker.chmod(0o755)
    jobs={}
    def job(stage,n,cpus,mem,limit,cap=None,prereq=(),part="cpu",aftercorr=None):
        ar="0-"+str(n-1)+(("%"+str(cap)) if cap else "")
        cmd=["sbatch","--parsable","--partition="+part,"--array="+ar,
             "--cpus-per-task="+str(cpus),"--mem="+mem,"--time="+limit,
             "--job-name=OxJ_"+stage,
             "--output="+str(output/"logs"/(stage+".%A_%a.out")),
             "--error="+str(output/"logs"/(stage+".%A_%a.err"))]
        deps=[]
        if prereq:
            deps.append("afterok:"+":".join(prereq))
        if aftercorr:
            deps.append("aftercorr:"+str(aftercorr))
        if deps:
            cmd.append("--dependency="+",".join(deps))
        jid=subprocess.check_output(cmd+[str(worker),stage],text=True).strip().split(";")[0]
        jobs[stage]=jid
        write(output/"jobs.json",jobs)
        print("JOB\t"+stage+"\t"+jid,flush=True)
        return jid
    try:
        a0=job("setup",1,4,"16G","04:00:00")
        a1=job("download",48,1,"4G","20:00:00",8)
        # Each alignment can start as soon as its own mate pair is validated.
        # setup must be complete for every task, but other downloads need not finish.
        a2=job("align",48,16,"64G","3-00:00:00",4,[a0],aftercorr=a1)
        a3=job("haplotype",88,4,"16G","12:00:00",10,[a2])
        a4=job("joint",2,4,"24G","12:00:00",2,[a3])
        a5=job("focal",2,2,"12G","06:00:00",2,[a4])
        job("audit",1,2,"12G","02:00:00",None,[a5],"cpu")
    except Exception:
        for jid in jobs.values():
            subprocess.run(["scancel",jid],capture_output=True)
        raise
    (output.parent/"latest_submission.path").write_text(str(output)+"\n")
    print("RUN_DIR\t"+str(output))
    print("SCOPE\t44 comparative WGS; 4 focal; 5 cythera; 5 vulcanus")
    print("REFERENCE\t"+str(original_fasta))


def exe(cfg,name):
    x=Path(cfg["env"])/"bin"/name
    if not x.is_file():
        raise RuntimeError("Tool not installed: "+str(x))
    return str(x)


def setup(cfg,_):
    base=Path(cfg["root"])
    env=Path(cfg["env"])
    if not (env/"bin/gatk").is_file():
        run(["micromamba","create","-y","-p",env,"-c","conda-forge",
             "-c","bioconda","gatk4","bwa","samtools","bcftools","picard"])
    for name in ("gatk","bwa","samtools","bcftools","picard"):
        exe(cfg,name)
    ref=cfg["fasta"]
    if not Path(ref+".fai").is_file():
        run([exe(cfg,"samtools"),"faidx",ref])
    d=Path(ref).with_suffix(".dict")
    if not d.is_file():
        run([exe(cfg,"gatk"),"CreateSequenceDictionary","-R",ref,"-O",d])
    if not Path(ref+".bwt").is_file():
        run([exe(cfg,"bwa"),"index",ref])
    write(base/"tool_versions.json",{name:subprocess.run(
        [exe(cfg,name),"--version"],capture_output=True,text=True).stdout[:200]
        for name in ("gatk","samtools","bcftools")})
    (base/"ref/complete").write_text("OK\n")
    print("SETUP_COMPLETE")


def ena_mates(acc):
    params=urlencode(dict(accession=acc,result="read_run",format="tsv",
        fields="run_accession,library_layout,fastq_ftp,fastq_md5,fastq_bytes"))
    with urlopen("https://www.ebi.ac.uk/ena/portal/api/filereport?"+params,timeout=90) as r:
        rows=list(csv.DictReader(r.read().decode().splitlines(),delimiter="\t"))
    matches=[r for r in rows if r.get("run_accession")==acc and r.get("library_layout")=="PAIRED"]
    if len(matches)!=1:
        raise RuntimeError("ENA run missing/ambiguous "+acc)
    r=matches[0]
    urls=(r.get("fastq_ftp") or "").split(";")
    sums=(r.get("fastq_md5") or "").split(";")
    sizes=(r.get("fastq_bytes") or "").split(";")
    if not len(urls)==len(sums)==len(sizes):
        raise RuntimeError("ENA FASTQ fields mismatch "+acc)
    mates={}
    for url,checksum,size in zip(urls,sums,sizes):
        m=re.search(r"_([12])\.(?:fastq|fq)\.gz$",url)
        if not m:continue
        mate=int(m.group(1))
        if mate in mates or not re.fullmatch(r"[0-9a-fA-F]{32}",checksum) or not size.isdigit():
            raise RuntimeError("ENA mate metadata invalid "+acc)
        mates[mate]=(url,checksum,int(size))
    if set(mates)!={1,2}:
        raise RuntimeError("No unique ENA pair for "+acc)
    return mates


def download(cfg,index):
    row=cfg["samples"][index]
    name=row["sample"]
    d=Path(cfg["root"])/"raw"/name
    d.mkdir(exist_ok=True)
    done=d/"complete.json"
    if done.exists() and all((d/(name+"_{}.fastq.gz".format(i))).is_file() for i in (1,2)):
        print("READS_REUSED\t"+name)
        return
    prior=[Path(row["r1"]),Path(row["r2"])] if row["r1"] and row["r2"] else [None,None]
    acc=row["ena_run"] or (name if re.fullmatch(r"[SED]RR\d+",name) else "")
    meta={}
    if row["fastq_r1"] and row["fastq_r2"]:
        for i in (1,2):
            meta[i]=(row["fastq_r"+str(i)],row["md5_r"+str(i)],
                     int(row["bytes_r"+str(i)]))
    elif acc:
        meta=ena_mates(acc)
    notes=[]
    for i in (1,2):
        target=d/(name+"_{}.fastq.gz".format(i))
        source=prior[i-1]
        vals=meta.get(i)
        def valid(path):
            return bool(path and path.is_file() and path.stat().st_size>0 and
                        (not vals or (path.stat().st_size==vals[2] and md5(path)==vals[1])))
        if valid(target):
            pass
        elif valid(source):
            target.unlink(missing_ok=True)
            target.symlink_to(source.resolve())
        elif vals:
            url=vals[0]
            if not url.startswith("https://"):
                url="https://"+url.removeprefix("ftp://")
            if not url.startswith("https://ftp.sra.ebi.ac.uk/"):
                raise RuntimeError("Unexpected ENA download domain: "+url)
            tmp=Path(str(target)+".partial")
            run(["curl","--fail","--location","--retry","5",
                 "--silent","--show-error","--output",tmp,url])
            if not valid(tmp):
                tmp.unlink(missing_ok=True)
                raise RuntimeError("FASTQ checksum mismatch "+name)
            tmp.replace(target)
        else:
            raise RuntimeError("No available FASTQ for "+name+" mate "+str(i))
        if not valid(target):
            raise RuntimeError("FASTQ validation failed "+name)
        notes.append(dict(mate=i,bytes=target.stat().st_size,
                          md5_verified=bool(vals),path=str(target.resolve())))
    write(done,dict(sample=name,ena_run=acc,mates=notes))
    print("READS_READY\t"+name)


def align(cfg,index):
    row=cfg["samples"][index]
    name=row["sample"]
    base=Path(cfg["root"])
    dst=base/"bam"/(name+".dedup.bam")
    done=base/"bam"/(name+".complete.json")
    if dst.is_file() and Path(str(dst)+".bai").exists() and done.exists():
        print("BAM_REUSED\t"+name)
        return
    d=base/"raw"/name
    need(d/"complete.json")
    with tempfile.TemporaryDirectory(prefix=name+".",dir=base/"bam") as tmp:
        t=Path(tmp)
        sortedbam=t/"sorted.bam"
        bam=t/"dedup.bam"
        rg="@RG\\tID:"+name+"\\tSM:"+name+"\\tPL:ILLUMINA"
        cmd=[exe(cfg,"bwa"),"mem","-M","-t","16","-R",rg,
             cfg["fasta"],d/(name+"_1.fastq.gz"),d/(name+"_2.fastq.gz")]
        srt=[exe(cfg,"samtools"),"sort","-@","4","-m","2G","-o",sortedbam,"-"]
        with subprocess.Popen(list(map(str,cmd)),stdout=subprocess.PIPE) as p1:
            with subprocess.Popen(list(map(str,srt)),stdin=p1.stdout) as p2:
                p1.stdout.close()
                status2=p2.wait()
            status1=p1.wait()
        if status1 or status2:
            raise RuntimeError("bwa/samtools failed for "+name)
        run([exe(cfg,"picard"),"MarkDuplicates","I="+str(sortedbam),
             "O="+str(bam),"M="+str(t/"metrics.txt"),
             "REMOVE_DUPLICATES=false","VALIDATION_STRINGENCY=LENIENT"])
        run([exe(cfg,"samtools"),"index","-@","4",bam])
        run([exe(cfg,"samtools"),"quickcheck",bam])
        bam.replace(dst)
        Path(str(bam)+".bai").replace(str(dst)+".bai")
        (t/"metrics.txt").replace(base/"bam"/(name+".duplicate_metrics.txt"))
    write(done,dict(sample=name,reference=cfg["original_fasta"],
                    full_genome_alignment=True,duplicate_marking=True))
    print("BAM_READY\t"+name)


def haplotype(cfg,index):
    base=Path(cfg["root"])
    names=cfg["comparative"]
    loci=list(cfg["loci"])
    locus=loci[index//len(names)]
    sample=names[index%len(names)]
    d=base/"gvcf"/locus
    d.mkdir(exist_ok=True)
    dest=d/(sample+".g.vcf.gz")
    if dest.is_file() and Path(str(dest)+".tbi").is_file():
        return
    chrom,start,end=cfg["loci"][locus]
    bam=need(base/"bam"/(sample+".dedup.bam"))
    temp=d/(sample+".partial.g.vcf.gz")
    run([exe(cfg,"gatk"),"--java-options","-Xmx12g","HaplotypeCaller",
         "-R",cfg["fasta"],"-I",bam,"-O",temp,"-ERC","GVCF",
         "-L",f"{chrom}:{start}-{end}","--sample-ploidy","2",
         "--native-pair-hmm-threads","4"])
    need(temp)
    indexfile=need(str(temp)+".tbi")
    temp.replace(dest)
    indexfile.replace(str(dest)+".tbi")
    print("GVCF_READY\t"+locus+"\t"+sample)


def joint(cfg,index):
    base=Path(cfg["root"])
    locus=list(cfg["loci"])[index]
    chrom,beg,end=cfg["loci"][locus]
    d=base/"joint"/locus
    d.mkdir(exist_ok=True)
    outfile=d/(locus+".all_sites.vcf.gz")
    if outfile.is_file() and (d/"complete.json").is_file():
        return
    merged=d/(locus+".merged.g.vcf.gz")
    paths=[need(base/"gvcf"/locus/(n+".g.vcf.gz")) for n in cfg["comparative"]]
    if not merged.is_file():
        cmd=[exe(cfg,"gatk"),"--java-options","-Xmx20g","CombineGVCFs",
             "-R",cfg["fasta"],"-O",merged,"-L",f"{chrom}:{beg}-{end}"]
        for path in paths:
            cmd.extend(["--variant",path])
        run(cmd)
    run([exe(cfg,"gatk"),"--java-options","-Xmx20g","GenotypeGVCFs",
         "-R",cfg["fasta"],"-V",merged,"-O",outfile,
         "-L",f"{chrom}:{beg}-{end}","--include-non-variant-sites","true"])
    import pysam
    with pysam.VariantFile(str(need(outfile))) as v:
        samples=set(v.header.samples)
    if samples!=set(cfg["comparative"]):
        raise RuntimeError("Joint VCF sample list mismatch "+locus)
    write(d/"complete.json",dict(locus=locus,samples=len(samples),
                                   ref_sites_retained=True))
    print("JOINT_READY\t"+locus+"\t"+str(len(samples)))


def focal(cfg,index):
    import pysam
    base=Path(cfg["root"])
    locus=list(cfg["loci"])[index]
    chrom,start,end=cfg["loci"][locus]
    var=need(base/"joint"/locus/(locus+".all_sites.vcf.gz"))
    sites={}
    with pysam.VariantFile(str(var)) as v:
        for rec in v.fetch(chrom,start-1,end):
            if len(rec.ref)==1 and rec.alts and len(rec.alts)==1 and (
                rec.ref.upper() in "ACGT" and rec.alts[0].upper() in "ACGT"):
                sites[rec.pos]=(rec.ref.upper(),rec.alts[0].upper())
    if not sites:
        raise RuntimeError("No biallelic SNPs at "+locus)
    per={}
    for name in cfg["focal"]:
        counts=defaultdict(lambda:[0,0])
        with pysam.AlignmentFile(str(need(base/"bam"/(name+".dedup.bam"))),"rb") as bam:
            for col in bam.pileup(chrom,start-1,end,truncate=True,
                stepper="samtools",min_mapping_quality=30,min_base_quality=20,
                flag_filter=1796,ignore_overlaps=True,max_depth=100000):
                pos=col.reference_pos+1
                if pos not in sites:continue
                ref,alt=sites[pos]
                for r in col.pileups:
                    if r.is_del or r.is_refskip or r.query_position is None:continue
                    b=r.alignment.query_sequence[r.query_position].upper()
                    if b==ref:counts[pos][0]+=1
                    elif b==alt:counts[pos][1]+=1
        per[name]=counts
    sides={r["sample"]:r["side"] for r in cfg["samples"] if r["role"]=="focal"}
    outfile=base/"focal"/(locus+".allele_counts.tsv.gz")
    with gzip.open(outfile,"wt") as fh:
        w=csv.writer(fh,delimiter="\t")
        w.writerow(["chrom","pos","ref","alt","Left_ref","Left_alt",
                    "Right_ref","Right_alt","Left_AF","Right_AF",
                    "per_library_counts_json"])
        for pos,(ref,alt) in sorted(sites.items()):
            counts={"Left":[0,0],"Right":[0,0]}
            details={}
            for name in cfg["focal"]:
                pair=per[name].get(pos,[0,0])
                details[name]=pair
                side=sides[name]
                counts[side][0]+=pair[0]
                counts[side][1]+=pair[1]
            l=counts["Left"];r=counts["Right"]
            frac=lambda a:(a[1]/sum(a)) if sum(a) else ""
            w.writerow([chrom,pos,ref,alt,*l,*r,frac(l),frac(r),
                        json.dumps(details,sort_keys=True)])
    write(base/"focal"/(locus+".complete.json"),
          dict(locus=locus,snps=len(sites),focal_gt_not_called=True,
               combination="sum T2+T3 per side"))
    print("FOCAL_COUNTS_READY\t"+locus+"\t"+str(len(sites)))


def audit(cfg,_):
    import pysam
    base=Path(cfg["root"])
    rows=[]
    for locus in cfg["loci"]:
        need(base/"joint"/locus/"complete.json")
        need(base/"focal"/(locus+".complete.json"))
        vcf=base/"joint"/locus/(locus+".all_sites.vcf.gz")
        totals={n:Counter() for n in cfg["comparative"]}
        with pysam.VariantFile(str(vcf)) as v:
            if set(v.header.samples)!=set(cfg["comparative"]):
                raise RuntimeError("Incorrect cohort in "+locus)
            for rec in v:
                for n in totals:
                    s=rec.samples[n]; t=totals[n]
                    t["sites"]+=1
                    gt=s.get("GT")
                    if not gt or any(x is None for x in gt):
                        t["missing"]+=1
                        continue
                    t["GT_"+"/".join(map(str,gt))]+=1
                    dp=s.get("DP");gq=s.get("GQ")
                    if dp is None:t["missing_DP"]+=1
                    if gq is None:t["missing_GQ"]+=1
                    if dp is not None and gq is not None and dp>=2 and gq>=20:
                        t["callable"]+=1
        for name,t in totals.items():
            rows.append(dict(locus=locus,sample=name,sites=t["sites"],
                GT_00=t["GT_0/0"],GT_01=t["GT_0/1"]+t["GT_1/0"],
                GT_11=t["GT_1/1"],missing=t["missing"],
                missing_DP=t["missing_DP"],missing_GQ=t["missing_GQ"],
                callable_DP2_GQ20=t["callable"]))
    tsv(base/"joint"/"genotype_QC.tsv",rows,list(rows[0]))
    if any(r["sites"]==0 for r in rows):
        raise RuntimeError("Empty joint VCF")
    write(base/"audit.complete.json",dict(
        joint_call_completed=True,comparative=44,focal_libraries=4,
        note="Next: locus-specific SNP ascertainment, reciprocal mixing, PDF figures"))
    print("AUDIT_COMPLETE\t"+str(base/"joint"/"genotype_QC.tsv"))


STAGES=dict(setup=setup,download=download,align=align,
            haplotype=haplotype,joint=joint,focal=focal,audit=audit)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--submit",action="store_true")
    p.add_argument("--preview",action="store_true",
                   help="Validate/present the 44-specimen cohort without submitting")
    p.add_argument("--preflight",type=Path)
    p.add_argument("--candidates",type=Path)
    p.add_argument("--sweep",type=Path)
    p.add_argument("--output",type=Path)
    p.add_argument("--config",type=Path)
    p.add_argument("--stage",choices=STAGES)
    p.add_argument("--task",type=int)
    a=p.parse_args()
    if a.submit or a.preview:
        if not (a.preflight and a.candidates and a.sweep):
            p.error("--submit/--preview needs --preflight, --candidates, --sweep")
        submit(a)
    elif a.config and a.stage and a.task is not None:
        STAGES[a.stage](read(a.config),a.task)
    else:
        p.error("Choose --submit or --config/--stage/--task")


if __name__=="__main__":
    try:main()
    except Exception as e:
        print("OPTIX_JOINT_ERROR\t"+str(e),file=sys.stderr)
        sys.exit(2)
