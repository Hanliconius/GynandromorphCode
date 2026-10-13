#!/usr/bin/env python3
"""Exact whole-genome placement of published CRE guide/PAM sequences."""
import argparse,csv,json
from pathlib import Path
QUERIES={"Del1":"GCTGTAATATGTTGCATATAGGG","Del2":"GGGTATATTTTTCAATACTAAGG","Del3":"ATAATACAGTATAACGTACGTGG","Del4":"GGTGGCATGGCATGGTTCATGG","primer1":"CGCTTACAGTCGCACTTTGT","primer2":"AATGTCGCCATGTGTACTCG"}
def records(path):
    name=None;parts=[]
    with path.open() as h:
        for line in h:
            if line.startswith(">"):
                if name is not None:yield name,"".join(parts).upper()
                name=line[1:].split()[0];parts=[]
            else:parts.append(line.strip())
        if name is not None:yield name,"".join(parts).upper()
def main():
    p=argparse.ArgumentParser();p.add_argument("--focal-run",type=Path,required=True);a=p.parse_args()
    cfg=json.loads((a.focal_run/"run.json").read_text());rows=[]
    for ref in cfg["references"]:
        for chrom,seq in records(a.focal_run/ref["name"]/"reference/genome.fa"):
            for name,target in QUERIES.items():
                reverse=target.translate(str.maketrans("ACGT","TGCA"))[::-1]
                for strand,needle in (("+",target),("-",reverse)):
                    pos=seq.find(needle)
                    while pos>=0:
                        rows.append([ref["name"],chrom,name,pos,pos+len(target),strand,chrom==ref["contig"]])
                        pos=seq.find(needle,pos+1)
    out=a.focal_run/"luca_CRE_mapping";out.mkdir(exist_ok=True)
    with (out/"exact_matches.tsv").open("w") as h:
        w=csv.writer(h,delimiter="\t");w.writerow(["reference","contig","query","start0","end0","strand","expected_Yb_contig"]);w.writerows(rows)
    for ref in cfg["references"]:
        for name in QUERIES:
            hits=[r for r in rows if r[0]==ref["name"] and r[2]==name]
            print(ref["name"],name,len(hits),hits,sep="\t")
    (out/"source.txt").write_text("Figure 4 source data 2: https://cdn.elifesciences.org/articles/68549/elife-68549-fig4-data2-v1.xlsx\nExperiment pairs: Del1/Del2 and Del3/Del4. Exact placements are targets, not yet functional boundaries. Del4 preserved as published (19-base protospacer).\n")
    print("MATCHES\t"+str(out/"exact_matches.tsv"))
if __name__=="__main__":main()
