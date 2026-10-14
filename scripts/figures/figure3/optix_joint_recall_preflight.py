#!/usr/bin/env python3
"""Read-only all-sample Optix/Yb cohort preflight; does NOT schedule or download."""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
PK_EXCLUDED = "SRR11243345"  # P4009: user-requested exclusion
EXTRAS = {"ERR2298205", "ERR2298206", "ERR2298204", "ERR2298214"}
SEEDS = [
    ("cythera", "cyth.CJ2856", "ERS977691"),
    ("cythera", "cyth.CJ2857", "ERS977692"),
    ("vulcanus", "vul.CJ14632", "ERS977709"),
    ("vulcanus", "vul.CS519", "ERS977710"),
    ("vulcanus", "vul.CS10", "ERS1030540"),
    ("vulcanus", "vul.CS11", "ERS1030541"),
]
FIELDS = ("sample", "specimen", "role", "group", "label", "ena_run",
          "ena_sample", "r1", "r2", "source", "notes")
VALID_NAME = re.compile(r"^[A-Za-z0-9_.-]+$")


def die(message):
    raise RuntimeError(message)


def read_json(path):
    path = Path(path)
    if not path.is_file():
        die(f"Missing run config: {path}")
    return json.loads(path.read_text())


def read_table(path):
    with Path(path).open(newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def check_name(value, where):
    if not value or not VALID_NAME.fullmatch(value):
        die(f"Unsafe or missing {where}: {value!r}")
    return value


def add(rows, sample, specimen, role, group, label, ena_run, ena_sample,
        r1, r2, source, notes=""):
    sample = check_name(sample, "sample")
    specimen = check_name(specimen, "specimen")
    rows.append(dict(sample=sample, specimen=specimen, role=role,
                     group=group, label=label, ena_run=ena_run,
                     ena_sample=ena_sample, r1=str(r1 or ""), r2=str(r2 or ""),
                     source=source, notes=notes))


def run_accession_metadata(accession):
    """Resolve archived run(s) for a BioSample or run; fail rather than invent."""
    params = urlencode({"accession": accession, "result": "read_run",
                        "fields": ("run_accession,sample_accession,library_strategy,"
                                   "library_layout,fastq_ftp,fastq_md5,fastq_bytes"),
                        "format": "tsv"})
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + params
    with urlopen(url, timeout=90) as response:
        records = list(csv.DictReader(io.StringIO(response.read().decode()),
                                      delimiter="\t"))
    accepted = []
    for r in records:
        if r.get("library_layout") != "PAIRED" or r.get("library_strategy") != "WGS":
            continue
        if not all(r.get(k) for k in ("run_accession", "sample_accession",
                                     "fastq_ftp", "fastq_md5", "fastq_bytes")):
            continue
        paths, md5s, sizes = [r[k].split(";") for k in
                               ("fastq_ftp", "fastq_md5", "fastq_bytes")]
        if len(paths) != 2 or len(md5s) != 2 or len(sizes) != 2:
            continue
        if not all(s.isdigit() for s in sizes):
            continue
        if not paths[0].endswith(("_1.fastq.gz", "_1.fq.gz")) or not paths[1].endswith(
                ("_2.fastq.gz", "_2.fq.gz")):
            continue
        accepted.append(r)
    accepted.sort(key=lambda r: (-sum(map(int, r["fastq_bytes"].split(";"))),
                                 r["run_accession"]))
    return accepted


def audit_reference(ref, allele_path):
    import pysam
    if not ref.is_file() or not Path(str(ref) + ".fai").is_file():
        die(f"Missing reference/FAI: {ref}")
    if not allele_path.is_file():
        die(f"Missing frozen allele file: {allele_path}")
    with pysam.FastaFile(str(ref)) as fa:
        for chrom in ("NC_140406.1", "NC_140403.1"):
            if chrom not in fa.references:
                die(f"Rosina contig not present: {chrom}")
        with gzip.open(allele_path, "rt") as h:
            sites = {}
            for line in h:
                chrom, pos, alleles = line.rstrip("\n").split("\t")[:3]
                if chrom != "NC_140406.1":
                    die("Frozen Optix alleles are not on the native NC_140406.1 contig")
                pos = int(pos)
                refbase = alleles.split(",")[0].upper()
                if len(refbase) != 1:
                    die(f"Non-SNP frozen Optix ref at {pos}")
                if pos in sites and sites[pos] != refbase:
                    die(f"Conflicting frozen REF at {pos}")
                sites[pos] = refbase
        if len(sites) < 100:
            die(f"Only {len(sites)} frozen Optix sites; investigate")
        first, last = min(sites), max(sites)
        sequence = fa.fetch("NC_140406.1", first - 1, last).upper()
        mismatch = [p for p, allele in sites.items()
                    if sequence[p - first] != allele]
        if mismatch:
            die(f"FROZEN_REF_MISMATCH {len(mismatch)} sites, first {mismatch[:5]}")
    digest = hashlib.sha256(allele_path.read_bytes()).hexdigest()
    return dict(reference=str(ref.resolve()), contigs=["NC_140406.1", "NC_140403.1"],
                frozen_optix_sites=len(sites), ref_agreement="all",
                frozen_alleles_sha256=digest)


def locate_candidate_text(root, out):
    """Bounded local search for metadata candidates, not variant/VCF contents."""
    roots = [root / "analysis/original_502",
             root / "analysis/genomewide_mixing_input/metadata"]
    findings = []
    for base in roots:
        if not base.is_dir():
            continue
        for p in base.rglob("*"):
            if not p.is_file() or p.suffix not in {".tsv", ".csv", ".txt", ".json"}:
                continue
            if p.stat().st_size > 8_000_000 or len(p.relative_to(base).parts) > 6:
                continue
            try:
                with p.open(errors="replace") as h:
                    for lineno, line in enumerate(h, 1):
                        if "cythera" in line.lower() or "vulcanus" in line.lower():
                            findings.append((str(p), lineno, line.strip()[:500]))
                            if len(findings) >= 10000:
                                break
            except (OSError, UnicodeError):
                pass
            if len(findings) >= 10000:
                break
    with out.open("w") as h:
        writer = csv.writer(h, delimiter="\t")
        writer.writerow(["file", "line", "evidence"])
        writer.writerows(findings)
    return len(findings)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pk-run", required=True, type=Path)
    p.add_argument("--mean-run", required=True, type=Path)
    p.add_argument("--extra-run", required=True, type=Path)
    p.add_argument("--new-races", type=Path, help="TSV: subspecies, biosample, [label, run_accession]")
    p.add_argument("--resolve-ena", action="store_true",
                   help="Ask ENA for paired-end WGS run/FASTQ metadata for race candidates")
    p.add_argument("--output", type=Path,
                   default=ROOT / "analysis/optix_joint_recall_preflight" /
                   datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    args = p.parse_args()
    out = args.output.resolve()
    if out.exists():
        die(f"Output exists; use new run directory: {out}")
    pk, old, extra = [read_json(x / "run.json") for x in
                      (args.pk_run, args.mean_run, args.extra_run)]
    raw = ROOT / "analysis/yb_extra_control_fastq"
    rows = []
    if len(pk["samples"]) != 21:
        die("Expected 21 old PK manifest entries before P4009 exclusion")
    for e in pk["samples"]:
        sample = e["sample"]
        if sample == PK_EXCLUDED:
            continue
        add(rows, sample, e["biosample"], "comparative",
            "PK_" + e["phenotype"], e["label"], sample, e["biosample"],
            raw / sample / (sample + "_1.fastq.gz"),
            raw / sample / (sample + "_2.fastq.gz"),
            str(args.pk_run / "run.json"))
    if sum(r["group"].startswith("PK_") for r in rows) != 20:
        die("PK cohort is not 20 unique retained specimens")

    def unique_rosina_tasks(cfg):
        matches = {}
        for t in cfg["tasks"]:
            if t.get("reference", {}).get("name") != "HmRos":
                continue
            sample = t["sample"]
            if sample in matches:
                die(f"Multiple Rosina tasks for {sample}")
            matches[sample] = t
        return matches

    old_tasks = unique_rosina_tasks(old)
    if len(old_tasks) != 14:
        die(f"Mean-run expected 10 controls + 4 focal, found {len(old_tasks)}")
    for sample, t in old_tasks.items():
        focal = t.get("role") == "focal"
        group = "focal" if focal else t.get("subspecies", "unclassified")
        add(rows, sample, "TN_gynandromorph" if focal else sample,
            "focal" if focal else "comparative", group,
            t.get("library_id", sample), "", "", t.get("r1"), t.get("r2"),
            str(args.mean_run / "run.json"))
    if sum(r["role"] == "focal" for r in rows) != 4:
        die("Expected exactly four focal tissue libraries")
    extra_tasks = unique_rosina_tasks(extra)
    if not EXTRAS.issubset(extra_tasks):
        die(f"Extra control samples missing: {EXTRAS - set(extra_tasks)}")
    for sample in sorted(EXTRAS):
        t = extra_tasks[sample]
        add(rows, sample, sample, "comparative",
            t.get("subspecies", "unclassified"), sample,
            sample, "", t.get("r1"), t.get("r2"),
            str(args.extra_run / "run.json"))

    candidates = []
    if args.new_races:
        candidates = read_table(args.new_races)
    else:
        candidates = [dict(subspecies=species, label=label, biosample=bio,
                           run_accession="") for species, label, bio in SEEDS]
    race_stats = defaultdict(set)
    resolved_runs = []
    ena_issues = []
    for entry in candidates:
        race = entry.get("subspecies", "").strip().lower()
        if race not in ("cythera", "vulcanus"):
            die(f"Unknown candidate race {race!r}")
        biosample = check_name(entry.get("biosample", ""), "candidate BioSample")
        if biosample in race_stats[race]:
            die(f"Duplicate candidate BioSample {biosample}")
        race_stats[race].add(biosample)
        accession = entry.get("run_accession", "").strip()
        selected = None
        if args.resolve_ena:
            runrows = run_accession_metadata(accession or biosample)
            runrows = [x for x in runrows if x.get("sample_accession") == biosample]
            if not runrows:
                ena_issues.append(f"{race} {biosample}: no validated paired WGS FASTQ run")
                continue
            selected = runrows[0]
            accession = selected["run_accession"]
            resolved_runs.append(dict(race=race, biosample=biosample, selected=accession,
                                      other_runs=[r["run_accession"] for r in runrows[1:]],
                                      ena=selected))
        if not accession:
            ena_issues.append(f"{race} {biosample}: run not yet resolved")
            continue
        check_name(accession, "ENA run")
        if not re.fullmatch(r"[SED]RR\d+", accession):
            die(f"Invalid ENA/SRA sequencing run accession: {accession}")
        d = ROOT / "analysis/optix_joint_raw_fastq" / accession
        add(rows, accession, biosample, "comparative", race,
            entry.get("label") or accession, accession, biosample,
            d / (accession + "_1.fastq.gz"), d / (accession + "_2.fastq.gz"),
            str(args.new_races) if args.new_races else "Martin_2016_archival_seeds",
            "Requires ENA checksum verification before use")

    names = Counter(x["sample"] for x in rows)
    duplicates = [s for s, n in names.items() if n != 1]
    if duplicates:
        die(f"Duplicate sample IDs: {duplicates}")
    specimens = defaultdict(list)
    for row in rows:
        if row["role"] == "comparative":
            specimens[row["specimen"]].append(row["sample"])
    duplicates = {s: ids for s, ids in specimens.items() if len(ids) > 1}
    if duplicates:
        die(f"Multiple comparative runs per BioSample: {duplicates}")
    reference = args.mean_run / "HmRos/reference/genome.fa"
    alleles = args.pk_run / "optix_alleles_rosina.tsv.gz"
    ref_audit = audit_reference(reference, alleles)

    out.mkdir(parents=True)
    with (out / "manifest.preflight.tsv").open("w") as h:
        writer = csv.DictWriter(h, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    with (out / "ena_selected_runs.json").open("w") as h:
        json.dump(resolved_runs, h, indent=2)
    (out / "ena_issues.txt").write_text("\n".join(ena_issues) + "\n")
    candidates_n = locate_candidate_text(ROOT, out / "original142_race_metadata_hits.tsv")
    available = sum(all(Path(r[k]).is_file() and Path(r[k]).stat().st_size > 0
                        for k in ("r1", "r2")) for r in rows)
    summary = dict(reference=ref_audit, cohort_count=len(rows),
                   comparative_count=len(specimens), focal_libraries=4,
                   race_candidate_biosamples={k: len(v) for k, v in race_stats.items()},
                   race_resolved_runs=Counter(r["group"] for r in rows
                                              if r["group"] in ("cythera", "vulcanus")),
                   ena_issues=ena_issues, fastq_pairs_present_unverified=available,
                   original142_metadata_matches=candidates_n,
                   requires_more_race_candidates=any(len(race_stats[r]) < 3
                                                     for r in ("cythera", "vulcanus")),
                   permitted_to_joint_call=False,
                   why_not_ready="No FASTQ digests, reference provenance for reuse, or tool versions verified")
    (out / "audit.json").write_text(json.dumps(summary, indent=2, default=dict) + "\n")
    print("PREFLIGHT_DIR\t" + str(out))
    print("COHORT\t" + str(len(rows)) + " libraries; " +
          str(len(specimens)) + " comparative specimens; 4 focal libraries")
    print("NEW_RACE_RUNS\t" + str(dict(summary["race_resolved_runs"])))
    print("ENA_UNRESOLVED\t" + str(len(ena_issues)))
    print("ANALYSIS142_CANDIDATE_LINES\t" + str(candidates_n))
    print("READY_FOR_SUBMISSION\tNO (input inventory and reference audit only)")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print("PREFLIGHT_ERROR\t" + str(error), file=sys.stderr)
        raise SystemExit(2)
