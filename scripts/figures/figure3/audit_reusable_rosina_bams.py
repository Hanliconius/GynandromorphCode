#!/usr/bin/env python3
"""Read-only reuse audit: existing whole-genome Rosina BAMs for Optix/Yb.

Search *known established pipeline output layouts*, not a recursive sweep of
the entire GPFS scratch filesystem. No Slurm changes, no BAM changes.
"""
import argparse
import csv
import json
from collections import Counter,defaultdict
from pathlib import Path
import pysam

ROLES={"comparative","focal"}
LOCI=("NC_140406.1","NC_140403.1")
PREVIOUS=("yb_local_both_references","yb_control_both_references",
          "yb_extra_controls")
FIELDS=("sample","role","group","candidate","source","status","reason",
        "optix_mapped","yb_mapped","long_contigs_mapped",
        "bam_bytes","has_picard_pg","picard_metrics_present",
        "source_ref_explicit_match","run_ref","original_task_index")


def load(path):
    return json.loads(Path(path).read_text())


def table(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def same_ref(source_path,reference):
    try:
        return Path(source_path).resolve(strict=True)==reference.resolve(strict=True)
    except (OSError,TypeError,ValueError):
        return False


def locate(root,wanted):
    """Yield candidate BAM files only from previous, known pipeline layouts."""
    found={}
    for group in PREVIOUS:
        parent=root/"analysis"/group
        if not parent.is_dir():
            continue
        for job in parent.iterdir():
            if not job.is_dir():
                continue
            config=job/"run.json"
            explicit_reference=False
            refpath=""
            if config.is_file():
                try:
                    content=load(config)
                    options=[x.get("source") or x.get("fasta") or ""
                             for x in content.get("references",[])
                             if x.get("name")=="HmRos"]
                    if options:refpath=options[0]
                except (ValueError,KeyError,TypeError):
                    pass
            folder=job/"HmRos"/"coverage"
            if not folder.is_dir():
                continue
            for bam in folder.glob("*.dedup.bam"):
                sample=bam.name[:-len(".dedup.bam")]
                if sample in wanted:
                    found[str(bam)]=(sample,group,refpath)
    # These earlier coverage runs store explicit BAM paths in their frozen
    # task manifest. Some paths point to chromosome-only 'locus.bam'; we
    # still inspect them, but the two-locus read-content check rejects them.
    cov=root/"analysis"/"yb_all_samples_coverage"
    if cov.is_dir():
        for work in cov.iterdir():
            config=work/"run.json"
            if not config.is_file():
                continue
            try:
                content=load(config)
            except (ValueError,TypeError):
                continue
            for item in content.get("tasks",[]):
                if not isinstance(item,dict):
                    continue
                sample=item.get("sample")
                if sample not in wanted or item.get("reference") not in ("HmRos",None):
                    continue
                bam=item.get("bam")
                if bam:
                    found.setdefault(str(Path(bam)),(sample,"yb_all_samples_coverage",""))
    return sorted(found.items())


def inspect_candidate(path,sample,group,ref_hint,target_sqs,source_fasta):
    p=Path(path)
    answer=dict(sample=sample,candidate=path,source=group,
                status="REJECT",reason="",optix_mapped=0,yb_mapped=0,
                bam_bytes=0,has_picard_pg="NO",
                picard_metrics_present="NO",long_contigs_mapped=0,
                source_ref_explicit_match="YES" if ref_hint and same_ref(ref_hint,source_fasta) else "NO",
                run_ref=ref_hint)
    if not p.is_file() or p.stat().st_size==0:
        answer["reason"]="file_absent"
        return answer
    answer["bam_bytes"]=p.stat().st_size
    # The established Yb BWA/Picard worker writes duplicate metrics to
    # <sample>.picard_metrics.txt and a <sample>.complete.json in coverage/.
    # BAM @PG is not guaranteed to include Picard even when it ran.
    answer["picard_metrics_present"]="YES" if (
        p.with_name(sample+".picard_metrics.txt").is_file()
        and p.with_name(sample+".complete.json").is_file()
    ) else "NO"
    if not Path(str(p)+".bai").is_file() and not p.with_suffix(".bai").is_file():
        answer["reason"]="missing_BAM_index"
        return answer
    try:
        with pysam.AlignmentFile(str(p),"rb") as reader:
            header=reader.header.to_dict()
            sqs={x["SN"]:x["LN"] for x in header.get("SQ",[])}
            if sqs!=target_sqs:
                answer["reason"]="reference_contig_dictionary_not_identical"
                return answer
            rg={x.get("SM") for x in header.get("RG",[])}
            if rg!={sample}:
                answer["reason"]="incorrect_BAM_read_group_sample"
                return answer
            pg=header.get("PG",[])
            has_markdup=any(
                "MarkDuplicates" in " ".join(str(x.get(k,"")) for k in
                   ("ID","PN","CL")) for x in pg)
            answer["has_picard_pg"]="YES" if has_markdup else "NO"
            mapped={s.contig:s.mapped for s in reader.get_index_statistics()}
            optix=mapped.get(LOCI[0],0)
            yb=mapped.get(LOCI[1],0)
            answer["optix_mapped"]=optix
            answer["yb_mapped"]=yb
            answer["long_contigs_mapped"]=sum(
                1 for name,length in target_sqs.items()
                if length>=100000 and mapped.get(name,0)>=20
            )
            if not optix or not yb:
                answer["reason"]="no_mapped_reads_on_both_target_contigs"
                return answer
            if answer["long_contigs_mapped"]<10:
                answer["reason"]="not_proven_genomewide_mapping"
                return answer
            # This is evidence of compatibility, NOT automatic permission to
            # substitute the BAM: explicit source FASTA provenance, read
            # identity and duplicate-marking remain to be reviewed.
            answer["status"]="POTENTIAL_REUSE"
            answer["reason"]="BAM_index_both_contigs_RG_and_dictionary_match"
    except (OSError,ValueError,RuntimeError) as e:
        answer["reason"]="BAM_open_or_index_error:"+str(e)[:130]
    return answer


def audit(source):
    source=source.resolve()
    cfg=load(source/"run.json")
    ref=Path(cfg["original_fasta"]).resolve()
    originals=table(source/"manifest.tsv")
    if len(originals)!=48:
        raise RuntimeError("Expected the frozen 48-row source manifest")
    wanted={r["sample"]:dict(r,index=i) for i,r in enumerate(originals)
            if r["sample"]!="SRR11243351"}
    if len(wanted)!=47:
        raise RuntimeError("Unexpected preserved cohort; P4003 exclusion mismatches")
    with pysam.FastaFile(str(ref)) as fa:
        target_sqs=dict(zip(fa.references,fa.lengths))
    all_candidates=[]
    for pathname,(sample,group,ref_hint) in locate(
        Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"),
        set(wanted)):
        result=inspect_candidate(pathname,sample,group,ref_hint,target_sqs,ref)
        row=wanted[sample]
        result.update(role=row["role"],group=row["group"],
                      original_task_index=row["index"])
        all_candidates.append(result)
    out=source/"existing_bam_reuse_audit"
    out.mkdir(exist_ok=True)
    path=out/"candidate_audit.tsv"
    with path.open("w",newline="") as fh:
        writer=csv.DictWriter(fh,fieldnames=FIELDS,extrasaction="ignore",
                              delimiter="\t")
        writer.writeheader()
        writer.writerows(all_candidates)
    samples=defaultdict(list)
    for row in all_candidates:
        samples[row["sample"]].append(row)
    good={s for s,ls in samples.items()
          if any(row["status"]=="POTENTIAL_REUSE" for row in ls)}
    strict={s for s,ls in samples.items()
            if any(row["status"]=="POTENTIAL_REUSE" and
                   row["source_ref_explicit_match"]=="YES" and
                   (row["has_picard_pg"]=="YES" or
                    row["picard_metrics_present"]=="YES") for row in ls)}
    missing=set(wanted)-good
    summary=dict(
        total_retained=47,
        existing_candidate_files=len(all_candidates),
        potential_reuse_specimens=len(good),
        strong_provenance_specimens=len(strict),
        note="Older frozen HmRos references use references[].source not fasta; Picard metrics can establish duplicate marking without a @PG line.",
        still_need_alignment_or_manual_review=len(missing),
        rejected_candidate_reasons=dict(Counter(
            row["reason"] for row in all_candidates
            if row["status"]=="REJECT")),
        potential_samples=sorted(good),
        missing_samples=sorted(missing),
        output_tsv=str(path),
        method="Metadata, BAM header and index only; no BAM reads modified",
        caution="Candidate reuse must not be automatic without confirmed reference and read provenance",
    )
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("RETAINED\t47")
    print("CANDIDATE_BAMS\t"+str(len(all_candidates)))
    print("POTENTIALLY_REUSABLE\t"+str(len(good)))
    print("STRONG_REFERENCE_AND_PICARD_PROVENANCE\t"+str(len(strict)))
    print("MISSING_OR_UNVERIFIED\t"+str(len(missing)))
    print("EXAMPLES\t"+",".join(sorted(good)[:12]))
    print("REPORT\t"+str(path))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run",type=Path,required=True)
    args=p.parse_args()
    audit(args.run)


if __name__=="__main__":
    main()
