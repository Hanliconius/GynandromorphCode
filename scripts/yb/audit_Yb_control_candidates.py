#!/usr/bin/env python3
"""Find local Yb control files by accessions AND manifest-resolved aliases."""
import csv,json,os,re,subprocess,shutil
from pathlib import Path
from collections import defaultdict,Counter
ROOT=Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe")
PROJECT=ROOT/"GynandroCalling"
HERE=Path(__file__).resolve().parent
ACCESSION=re.compile(r"\b(?:ERR|ERS|SRR|SRS|SRX|ERX|SAMN|SAMEA)\d+\b")
SKIP={".git","logs","tmp","cache","node_modules","__pycache__"}
def main():
    rows=list(csv.DictReader((HERE/"livraghi2021_and_krista_control_candidates.tsv").open(),delimiter="\t"))
    identifiers={r["accession"] for r in rows}
    aliases={a:{a} for a in identifiers}
    manifests=[];files=[];errors=[]
    print("YB_CONTROL_AUDIT_START",flush=True)
    for directory,dirs,names in os.walk(ROOT):
        dirs[:]=[d for d in dirs if d not in SKIP]
        for name in names:
            p=Path(directory)/name
            if name.endswith((".fastq",".fastq.gz",".fq",".fq.gz",".bam",".cram")):files.append(p)
            if name.endswith(".tsv") and any(t in name.lower() for t in ("manifest","samples","accession")):
                try:
                    if p.stat().st_size<=20_000_000:manifests.append(p)
                except OSError:pass
    matches=[];evidence=[]
    # Two passes resolve sample -> run -> local library labels across manifests.
    for _ in range(2):
        for p in manifests:
            try:
                with p.open() as handle:
                    reader=csv.DictReader(handle,delimiter="\t")
                    for row in reader:
                        values=[v for v in row.values() if isinstance(v,str)]
                        tokens=set().union(*(set(ACCESSION.findall(v)) for v in values))
                        for key in ("library_id","sample","sample_id","resolved_run","source_accession","accession"):
                            value=row.get(key,"")
                            if value and "/" not in value:tokens.add(value)
                        owners=[a for a in identifiers if tokens&aliases[a]]
                        for a in owners:
                            aliases[a].update(tokens)
                            for value in values:
                                if value.startswith("/") and value.endswith((".fastq",".fastq.gz",".fq",".fq.gz",".bam",".cram")):
                                    f=Path(value)
                                    if f.is_file() and f.stat().st_size:
                                        matches.append((a,str(f.resolve()),"manifest:"+str(p)))
                            evidence.append((a,str(p),json.dumps(row,sort_keys=True)))
            except (OSError,csv.Error,UnicodeError) as e:errors.append(str(p)+": "+str(e))
    for p in files:
        for a,labels in aliases.items():
            if any(re.search(r"(?<![A-Za-z0-9])"+re.escape(label)+r"(?![A-Za-z0-9])",str(p)) for label in labels):
                try:
                    if p.is_file() and p.stat().st_size:matches.append((a,str(p.resolve()),"filename/manifest_alias"))
                except OSError:pass
    matches=sorted(set(matches));byacc=defaultdict(set)
    for a,p,e in matches:byacc[a].add(p)
    out=PROJECT/"analysis/yb_control_audit";out.mkdir(parents=True,exist_ok=True)
    with (out/"local_file_candidates.tsv").open("w") as h:
        h.write("source_accession\tpath\tevidence\treference_status\n")
        for a,p,e in matches:h.write("\t".join((a,p,e,"unverified"))+"\n")
    with (out/"manifest_evidence.tsv").open("w") as h:
        h.write("accession\tmanifest\trow_json\n")
        for r in sorted(set(evidence)):h.write("\t".join(r)+"\n")
    for r in rows:
        paths=byacc[r["accession"]]
        r["local_fastq_files"]=sum(p.endswith((".fq",".fq.gz",".fastq",".fastq.gz")) for p in paths)
        r["local_bam_files"]=sum(p.endswith((".bam",".cram")) for p in paths)
        r["eligible_for_overlay"]=r["yellow_bar"]!="unresolved"
    columns=list(rows[0])
    with (out/"sample_inventory.tsv").open("w") as h:
        w=csv.DictWriter(h,fieldnames=columns,delimiter="\t");w.writeheader();w.writerows(rows)
    (out/"audit.json").write_text(json.dumps({"roots":[str(ROOT)],"manifest_count":len(manifests),
        "scanned_sequence_files":len(files),"source_rows":len(rows),"aliases":{a:sorted(v) for a,v in aliases.items()},
        "errors":errors,"BAM_reference_verified":False,"downloads_performed":False},indent=2)+"\n")
    for phenotype in ("present","absent","unresolved"):
        group=[r for r in rows if r["yellow_bar"]==phenotype]
        print(f"GROUP\t{phenotype}\tlisted_rows={len(group)}\twith_FASTQ={sum(bool(r['local_fastq_files']) for r in group)}\twith_BAM={sum(bool(r['local_bam_files']) for r in group)}")
    print(f"INVENTORY\t{out/'sample_inventory.tsv'}")
    print(f"FILES\t{out/'local_file_candidates.tsv'}")
    print("BAM_REFERENCE_STATUS\tunverified until header check; no automatic reuse")
if __name__=="__main__":main()
