#!/usr/bin/env python3
"""Resolve the SAME completed Rosina callset/BAMs for Yb and Optix analyses.

Read-only: does not launch SNP calls, align reads, or modify the old Yb runs.
By default reads the pointer written ONLY after successful genomewide QC.
"""
import argparse
import csv
import json
from pathlib import Path

ROOT=Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
SHARED=ROOT/"analysis/rosina_genomewide_joint_recall"


def load(path):
    return json.loads(Path(path).read_text())


def resolve(path=None):
    if path is None:
        pointer=SHARED/"latest_complete.path"
        if not pointer.is_file():
            raise RuntimeError("No completed genomewide Rosina recall is published")
        path=Path(pointer.read_text().strip())
    else:
        path=Path(path)
    meta=load(path/"complete.json")
    if meta["status"]!="COMPLETE" or meta["reference"]!="Rosina":
        raise RuntimeError("Shared run is incomplete or reference incorrect")
    yb=meta["yb"]
    if not yb["no_hmel25_snp_calls"]:
        raise RuntimeError("Yb output must be Rosina-only")
    with Path(yb["bam_manifest"]).open(newline="") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    if len(rows)!=47 or len([r for r in rows if r["role"]=="comparative"])!=43:
        raise RuntimeError("Shared Yb BAM manifest has wrong cohort")
    if any(r["sample"]=="SRR11243351" for r in rows):
        raise RuntimeError("Excluded P4003 reintroduced to callset")
    needed=[meta["genomewide_vcf"],meta["genomewide_snps"],
            yb["chrom_vcf"],yb["window_vcf"],yb["bam_manifest"]]
    if any(not Path(p).is_file() for p in needed):
        raise RuntimeError("Some published files are missing")
    return path,meta


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run",type=Path,
                   help="Optional completed run; default is latest_complete.path")
    p.add_argument("--json",action="store_true",
                   help="Print machine-readable complete run JSON")
    a=p.parse_args()
    path,meta=resolve(a.run)
    if a.json:
        print(json.dumps(meta,indent=2))
        return
    yb=meta["yb"]
    print("ROSINA_RUN\t"+str(path))
    print("REFERENCE\tRosina; Hmel2.5 SNP calling: none")
    print("SHARED_GENOMEWIDE_VCF\t"+meta["genomewide_vcf"])
    print("SHARED_GENOMEWIDE_SNPS\t"+meta["genomewide_snps"])
    print("YB_CHROMOSOME_VCF\t"+yb["chrom_vcf"])
    print("YB_WINDOW_VCF\t"+yb["window_vcf"])
    print("YB_BAM_MANIFEST\t"+yb["bam_manifest"])
    print("COMPARATIVE\t43; FOCAL LIBRARIES\t4")
    print("NOTE\tReuse Yb BAMs for coverage/median plots; never re-call SNPs on Hmel2.5")


if __name__=="__main__":
    main()
