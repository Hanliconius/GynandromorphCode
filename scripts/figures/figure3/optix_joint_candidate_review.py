#!/usr/bin/env python3
"""Resolve historical cythera/vulcanus candidate accessions; read-only on Pegasus.

Uses the archived preflight metadata hits, validates ENA WGS paired FASTQ
records, and prints an abbreviated summary. No downloads / Slurm submissions.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
SEEDS = [
    ("cythera", "cyth.CJ2856", "ERS977691"),
    ("cythera", "cyth.CJ2857", "ERS977692"),
    ("vulcanus", "vul.CJ14632", "ERS977709"),
    ("vulcanus", "vul.CS519", "ERS977710"),
    ("vulcanus", "vul.CS10", "ERS1030540"),
    ("vulcanus", "vul.CS11", "ERS1030541"),
]
ACCESSIONS = re.compile(r"(?<![A-Z0-9])((?:ERR|SRR|DRR|ERS|SRS|DRS)\d+)(?!\d)", re.I)
RACES = ("cythera", "vulcanus")
RUN_FIELDS = ["run_accession", "sample_accession", "library_strategy",
              "library_layout", "fastq_ftp", "fastq_md5", "fastq_bytes"]
OUTPUT_FIELDS = [
    "subspecies", "biosample", "label", "run_accession",
    "library_strategy", "library_layout", "fastq_r1", "fastq_r2",
    "md5_r1", "md5_r2", "bytes_r1", "bytes_r2", "source_evidence",
]


def write_tsv(path, rows, columns):
    with path.open("w", newline="") as h:
        writer = csv.DictWriter(h, columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with Path(path).open(newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def find_evidence(preflight, existing, output, top):
    """Collect accession hints; no historical GT values are ever ingested."""
    source = preflight / "original142_race_metadata_hits.tsv"
    hits = read_tsv(source)
    evidence = {}
    file_counts = Counter()
    for r in hits:
        line = r.get("evidence", "")
        race = next((x for x in RACES if re.search(r"\b" + x + r"\b", line, re.I)), None)
        if not race:
            continue
        file_counts[r["file"]] += 1
        accessions = set(m.group(1).upper() for m in ACCESSIONS.finditer(line))
        for accession in accessions:
            if accession in existing:
                continue
            key = (race, accession)
            obj = evidence.setdefault(key, {
                "race": race, "accession": accession,
                "n_lines": 0, "wgs_lines": 0, "cs_lines": 0,
                "files": set(), "first_source": "",
                "first_evidence": "", "source_is_manifest": False,
            })
            obj["n_lines"] += 1
            obj["files"].add(r["file"])
            obj["wgs_lines"] += bool(re.search(r"\bWGS\b", line, re.I))
            obj["cs_lines"] += bool(re.search(r"\bCS\b|\bcapture\b", line, re.I))
            obj["source_is_manifest"] |= "manifest" in r["file"].lower()
            if not obj["first_source"]:
                obj["first_source"] = r["file"] + ":" + r["line"]
                obj["first_evidence"] = line[:200]

    rows = []
    for e in evidence.values():
        score = (
            25 * bool(e["wgs_lines"])
            + 12 * e["source_is_manifest"]
            + 6 * e["accession"].startswith(("ERS", "SRS", "DRS"))
            + min(12, math.log2(1 + e["n_lines"]) * 2)
        )
        rows.append(dict(
            subspecies=e["race"], accession=e["accession"],
            score=round(score, 2), lines=e["n_lines"],
            wgs_lines=e["wgs_lines"], cs_lines=e["cs_lines"],
            sources=len(e["files"]), first_source=e["first_source"],
            example=e["first_evidence"],
        ))
    rows.sort(key=lambda x: (-x["score"], x["subspecies"], x["accession"]))
    write_tsv(output / "original142_accession_hints.tsv", rows,
              ["subspecies", "accession", "score", "lines",
               "wgs_lines", "cs_lines", "sources", "first_source", "example"])
    files = [dict(path=k, matching_lines=v) for k, v in file_counts.most_common()]
    write_tsv(output / "original142_matching_files.tsv", files,
              ["path", "matching_lines"])
    priority = defaultdict(list)
    for row in rows:
        if len(priority[row["subspecies"]]) < top:
            priority[row["subspecies"]].append(row)
    return len(hits), rows, priority, file_counts


def ena_record(query):
    params = urlencode({
        "accession": query, "result": "read_run", "format": "tsv",
        "fields": ",".join(RUN_FIELDS),
    })
    uri = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + params
    last_error = None
    for attempt in range(2):
        try:
            with urlopen(uri, timeout=35) as response:
                content = response.read(5_000_000).decode()
            records = list(csv.DictReader(content.splitlines(), delimiter="\t"))
            return records, ""
        except (HTTPError, URLError, TimeoutError, ValueError) as error:
            last_error = str(error)
            if attempt == 0:
                time.sleep(1)
    return [], "ENA_ERROR: " + str(last_error)


def parse_pair(row):
    if row.get("library_strategy", "").upper() != "WGS":
        return None
    if row.get("library_layout", "").upper() != "PAIRED":
        return None
    fields = [row.get(key, "").split(";") for key in
              ("fastq_ftp", "fastq_md5", "fastq_bytes")]
    paths, md5s, sizes = fields
    if not all(len(z) == 2 for z in fields):
        return None
    pairs = {}
    for p, checksum, size in zip(paths, md5s, sizes):
        if not re.fullmatch(r"[a-fA-F0-9]{32}", checksum) or not size.isdigit():
            return None
        if p.endswith(("_1.fastq.gz", "_1.fq.gz")):
            mate = 1
        elif p.endswith(("_2.fastq.gz", "_2.fq.gz")):
            mate = 2
        else:
            return None
        if mate in pairs:
            return None
        pairs[mate] = (p, checksum, size)
    if set(pairs) != {1, 2}:
        return None
    accession = row.get("run_accession", "")
    sample = row.get("sample_accession", "")
    if not re.fullmatch(r"[SED]RR\d+", accession):
        return None
    if not re.fullmatch(r"[SED]RS\d+", sample):
        return None
    return accession, sample, pairs


def lookup_one(item):
    race, label, query, origin = item
    records, error = ena_record(query)
    candidates = []
    for row in records:
        parsed = parse_pair(row)
        if parsed is None:
            continue
        run, sample, pairs = parsed
        # A seed BioSample must never resolve to a different BioSample.
        if query.startswith(("ERS", "SRS", "DRS")) and sample != query:
            continue
        if query.startswith(("ERR", "SRR", "DRR")) and run != query:
            continue
        candidates.append(dict(
            subspecies=race, biosample=sample, label=label or sample,
            run_accession=run, library_strategy="WGS",
            library_layout="PAIRED",
            fastq_r1=pairs[1][0], fastq_r2=pairs[2][0],
            md5_r1=pairs[1][1], md5_r2=pairs[2][1],
            bytes_r1=pairs[1][2], bytes_r2=pairs[2][2],
            source_evidence=origin,
        ))
    return item, candidates, error


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--preflight", required=True, type=Path)
    p.add_argument("--top-per-race", type=int, default=16)
    p.add_argument("--ena-workers", type=int, default=4)
    p.add_argument("--ena", action="store_true", help="Resolve archival candidate runs via ENA")
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    base = args.preflight.resolve()
    required = ("audit.json", "manifest.preflight.tsv",
                "original142_race_metadata_hits.tsv")
    if any(not (base / name).is_file() for name in required):
        p.error("Missing required preflight outputs in " + str(base))
    if not 1 <= args.top_per_race <= 80 or not 1 <= args.ena_workers <= 8:
        p.error("top-per-race must be 1..80; ena-workers 1..8")
    output = args.output or (ROOT / "analysis/optix_joint_candidate_review" /
                            datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    output = output.resolve()
    if output.exists():
        p.error("Output already exists: " + str(output))
    output.mkdir(parents=True)
    samples = read_tsv(base / "manifest.preflight.tsv")
    previous = {s["sample"] for s in samples}
    n, hints, priority, files = find_evidence(
        base, previous, output, args.top_per_race)
    print("EVIDENCE_LINES\t" + str(n))
    print("ACCESSION_HINTS\t" + str(len(hints)))
    print("TOP_METADATA_FILES")
    for path, count in files.most_common(6):
        print("  " + str(count) + "\t" + path)
    for race in RACES:
        print("HINTS_" + race.upper() + "\t" + str(
            len([h for h in hints if h["subspecies"] == race])))

    tools = {name: shutil.which(name) or ""
             for name in ("gatk", "java", "samtools", "bwa", "bcftools",
                          "picard", "bgzip", "tabix", "sbatch")}
    (output / "tools_on_login.json").write_text(json.dumps(tools, indent=2) + "\n")
    print("TOOLS\t" + " ".join(f"{name}={'yes' if where else 'NO'}"
                               for name, where in tools.items()))

    requests = {}
    for race, label, query in SEEDS:
        requests[query] = (race, label, query, "Martin_2016_archival_seed")
    for race in RACES:
        for hint in priority[race]:
            key = hint["accession"]
            if key not in requests:
                requests[key] = (race, "", key, hint["first_source"])
    verified = []
    errors = []
    if args.ena:
        with ThreadPoolExecutor(max_workers=args.ena_workers) as pool:
            futures = {pool.submit(lookup_one, item): item for item in requests.values()}
            for f in as_completed(futures):
                item, candidates, error = f.result()
                if error:
                    errors.append(dict(race=item[0], accession=item[2], reason=error))
                elif not candidates:
                    errors.append(dict(race=item[0], accession=item[2],
                                       reason="No paired WGS run with two validated FASTQ links"))
                verified.extend(candidates)
    else:
        print("ENA\tNot queried; pass --ena for validation")

    # Repeated runs per BioSample are alternatives, not independent individuals.
    by_biosample = defaultdict(list)
    for row in verified:
        by_biosample[row["biosample"]].append(row)
    ambiguous = [(bio, sorted({r["subspecies"] for r in group}))
                 for bio, group in by_biosample.items()
                 if len({r["subspecies"] for r in group}) > 1]
    if ambiguous:
        raise RuntimeError("Same BioSample assigned multiple races: " + repr(ambiguous))
    selected = []
    for bio, group in by_biosample.items():
        # Prefer a previously identified archival specimen, then larger paired FASTQ volume.
        group.sort(key=lambda r: (
            -int(r["source_evidence"] == "Martin_2016_archival_seed"),
            -int(r["bytes_r1"]) - int(r["bytes_r2"]), r["run_accession"]))
        selected.append(group[0])
    selected.sort(key=lambda r: (r["subspecies"], r["biosample"]))
    write_tsv(output / "all_resolved_runs.tsv", verified, OUTPUT_FIELDS)
    write_tsv(output / "verified_unique_wgs_candidates.tsv", selected, OUTPUT_FIELDS)
    write_tsv(output / "ena_unresolved.tsv", errors, ["race", "accession", "reason"])
    write_tsv(output / "candidate_races_for_preflight.tsv", [
        dict(subspecies=x["subspecies"], biosample=x["biosample"],
             label=x["label"], run_accession=x["run_accession"]) for x in selected
    ], ["subspecies", "biosample", "label", "run_accession"])
    for race in RACES:
        group = [x for x in selected if x["subspecies"] == race]
        print(f"VERIFIED_{race.upper()}\t{len(group)} unique WGS BioSamples")
        for x in group[:8]:
            print("  " + x["label"] + "\t" + x["biosample"] +
                  "\t" + x["run_accession"])
    print("ENA_UNRESOLVED\t" + str(len(errors)))
    print("CANDIDATE_DIR\t" + str(output))
    print("READY_3_PER_RACE\t" + ("YES" if args.ena and all(
        sum(r["subspecies"] == race for r in selected) >= 3 for race in RACES
    ) else "NO"))
    print("SUBMITTED_JOBS\t0")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("CANDIDATE_REVIEW_ERROR\t" + str(exc), file=sys.stderr)
        raise SystemExit(2)
