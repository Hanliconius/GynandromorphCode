#!/usr/bin/env python3
"""Read-only ENA project/taxonomy sweep and taxonomic audit for Optix races.

Purpose: find additional independent H. m. cythera WGS specimens beyond the
historical two, and flag vulcanus labels that are not directly supported by
ENA metadata. Does not download reads, submit jobs, or edit previous results.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, quote
from urllib.request import urlopen

import optix_joint_candidate_review as prior

ROOT = prior.ROOT
ENA_ROOT = "https://www.ebi.ac.uk"
FIELDS = prior.RUN_FIELDS + [
    "scientific_name", "sample_title", "sample_alias",
    "sample_description", "sub_species", "tax_id", "study_accession",
]
PROJECTS = ["PRJNA485723"]
EXTRA_SAMPLES = ["SAMEA6447028"]  # Darragh et al. genome-sample table
ANNOTATED = [
    "source_query", "subspecies", "biosample", "secondary_sample_accession",
    "run_accession", "scientific_name", "sample_title", "sample_alias",
    "sample_description", "sub_species", "tax_id", "study_accession",
    "library_strategy", "library_layout", "wgs_pair", "evidence",
    "race_match", "already_in_previous",
    "fastq_r1", "fastq_r2", "md5_r1", "md5_r2",
    "bytes_r1", "bytes_r2", "reason",
]
AUDIT_COLS = [
    "previous_group", "previous_biosample", "previous_run",
    "ena_primary_sample", "ena_secondary_sample",
    "ena_scientific_name", "ena_sample_title",
    "ena_sub_species", "status", "notes",
]


def field(record, key):
    return str(record.get(key) or "").strip()


def fetch_tsv(endpoint, params):
    uri = ENA_ROOT + endpoint + "?" + urlencode(params)
    err = ""
    for attempt in range(3):
        try:
            with urlopen(uri, timeout=55) as response:
                data = response.read(10_000_000).decode()
            reader = csv.DictReader(data.splitlines(), delimiter="\t")
            if not reader.fieldnames:
                return [], "NO_HEADER"
            return [
                {k: (v if isinstance(v, str) else "")
                 for k, v in record.items() if k is not None}
                for record in reader
            ], ""
        except (HTTPError, URLError, TimeoutError, ValueError) as exc:
            err = f"{type(exc).__name__}: {exc}"
            if attempt < 2:
                time.sleep(1 + attempt)
    return [], err


def ena_filereport(acc):
    return fetch_tsv("/ena/portal/api/filereport", {
        "accession": acc, "result": "read_run",
        "fields": ",".join(FIELDS), "format": "tsv",
    })


def ena_taxonomy_query():
    """Use exact taxon name, not a substring matching unrelated records."""
    url = ENA_ROOT + "/ena/taxonomy/rest/scientific-name/" + quote(
        "Heliconius melpomene cythera")
    with urlopen(url, timeout=55) as response:
        rows = json.loads(response.read(2_000_000).decode())
    exact = [r for r in rows if
             str(r.get("scientificName", "")).lower() ==
             "heliconius melpomene cythera"]
    if len(exact) != 1:
        raise RuntimeError(f"Exact cythera taxonomy not uniquely resolved: {exact}")
    tax_id = str(exact[0]["taxId"])
    return tax_id, fetch_tsv("/ena/portal/api/search", {
        "result": "read_run", "query": f"tax_eq({tax_id})",
        "fields": ",".join(FIELDS), "format": "tsv", "limit": "0",
    })


def race_evidence(row):
    scientific = field(row, "scientific_name").lower()
    title = " ".join(field(row, key) for key in (
        "sample_title", "sample_alias", "sample_description", "sub_species")).lower()
    target = "cythera"
    explicit = re.search(r"\bheliconius\s+melpomene\s+cythera\b", scientific)
    text_mention = bool(re.search(r"\bcythera\b", title))
    other = bool(re.search(r"\b(?:vulcanus|rosina|plesseni|malleti|ecuadorensis|amaryllis)\b",
                           scientific + " " + title))
    if explicit:
        return "cythera", "TAXON_CONFIRMED", "scientific_name"
    if text_mention and not other:
        return "cythera", "MANUAL_REVIEW", "title/alias/description only"
    return "", "NOT_CONFIRMED", "no direct cythera evidence"


def annotate(row, source, previous_sample_ids):
    run = field(row, "run_accession")
    bio = field(row, "sample_accession")
    secondary = field(row, "secondary_sample_accession")
    race, evidence, detail = race_evidence(row)
    parsed, reason = prior.parse_pair(row)
    record = {k: "" for k in ANNOTATED}
    for k in ("run_accession", "scientific_name", "sample_title",
              "sample_alias", "sample_description", "sub_species",
              "tax_id", "study_accession", "library_strategy", "library_layout"):
        record[k] = field(row, k)
    record.update(
        source_query=source,
        subspecies=race, biosample=bio,
        secondary_sample_accession=secondary,
        evidence=evidence,
        race_match=detail,
        already_in_previous=("yes" if bio in previous_sample_ids else "no"),
        wgs_pair=("yes" if parsed is not None else "no"),
        reason=reason,
    )
    if parsed is not None:
        _, _, _, pairs = parsed
        record.update(
            fastq_r1=pairs[1][0], fastq_r2=pairs[2][0],
            md5_r1=pairs[1][1], md5_r2=pairs[2][1],
            bytes_r1=pairs[1][2], bytes_r2=pairs[2][2],
        )
    return record


def audit_previous(previous, annotated):
    # Keep candidate race inferred from prior historical metadata separate
    # from the organism's actual BioSample taxonomy and submission metadata.
    indices = defaultdict(list)
    for row in annotated:
        if row["source_query"].startswith("audit:"):
            indices[row["run_accession"]].append(row)
    report = []
    for p in previous:
        result = indices.get(p["run_accession"], [])
        row = result[0] if result else {}
        primary = row.get("biosample", "")
        secondary = row.get("secondary_sample_accession", "")
        scientific = row.get("scientific_name", "")
        title = row.get("sample_title", "")
        sub = row.get("sub_species", "")
        source = " ".join((scientific, title, sub)).lower()
        group = p["subspecies"]
        notes = []
        status = "UNVERIFIED_TAXONOMY"
        if not row:
            status = "ENA_LOOKUP_MISSING"
        elif p["biosample"] != primary:
            status = "BIOSAMPLE_CHANGED"
        elif re.search(r"\b" + re.escape(group) + r"\b", source):
            status = "RACE_METADATA_SUPPORTED"
        else:
            conflicts = [x for x in ("cythera", "vulcanus", "rosina",
                                      "plesseni", "malleti", "ecuadorensis")
                         if x != group and re.search(r"\b" + x + r"\b", source)]
            if conflicts:
                status = "TAXON_CONFLICT"
                notes.append(",".join(conflicts))
        report.append(dict(
            previous_group=group, previous_biosample=p["biosample"],
            previous_run=p["run_accession"],
            ena_primary_sample=primary, ena_secondary_sample=secondary,
            ena_scientific_name=scientific, ena_sample_title=title,
            ena_sub_species=sub, status=status,
            notes="; ".join(notes),
        ))
    return report


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidate-run", required=True, type=Path,
                    help="Successful 20261014T203210Z_fixed discovery run")
    ap.add_argument("--output", type=Path)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    if not 1 <= args.workers <= 8:
        ap.error("workers must be 1..8")
    previous_path = args.candidate_run / "verified_unique_wgs_candidates.tsv"
    if not previous_path.is_file():
        ap.error(f"Missing previous verified sample table: {previous_path}")
    previous = prior.read_tsv(previous_path)
    out = args.output or (ROOT / "analysis/optix_joint_project_sweep" /
                          datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"))
    if out.exists():
        ap.error(f"Output already exists: {out}")
    out.mkdir(parents=True)
    prev_ids = {r["biosample"] for r in previous}
    queries = [(f"project:{s}", s) for s in PROJECTS]
    queries += [(f"sample:{s}", s) for s in EXTRA_SAMPLES]
    queries += [(f"audit:{r['subspecies']}:{r['run_accession']}",
                 r["run_accession"]) for r in previous]
    rows = []
    errors = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(ena_filereport, acc): label for label, acc in queries}
        for future in as_completed(futures):
            label = futures[future]
            try:
                records, error = future.result()
                if error:
                    errors.append(dict(source=label, problem=error))
                elif not records:
                    errors.append(dict(source=label, problem="NO_ROWS"))
                rows.extend(annotate(r, label, prev_ids) for r in records)
            except Exception as exc:
                errors.append(dict(source=label, problem=repr(exc)))
    try:
        taxid, (taxrows, taxerror) = ena_taxonomy_query()
        if taxerror:
            errors.append(dict(source="cythera_taxon:" + taxid,
                               problem=taxerror))
        rows.extend(annotate(r, "taxon:" + taxid, prev_ids) for r in taxrows)
    except Exception as exc:
        errors.append(dict(source="cythera_taxon", problem=repr(exc)))
    # Multiple retrieval routes return the same run: keep the strongest evidence.
    unique = {}
    for row in rows:
        key = (row["run_accession"], row["biosample"])
        if not key[0] or not key[1]:
            continue
        if key not in unique or (
            row["source_query"].startswith("project:") and
            not unique[key]["source_query"].startswith("project:")
        ):
            unique[key] = row
    dedup = list(unique.values())
    dedup.sort(key=lambda x: (x["source_query"], x["scientific_name"],
                              x["biosample"], x["run_accession"]))
    prior.write_tsv(out / "ena_all_project_taxon_runs.tsv", dedup, ANNOTATED)
    review = [
        r for r in dedup if r["subspecies"] == "cythera" and
        r["wgs_pair"] == "yes" and r["already_in_previous"] == "no"
    ]
    by_bio = defaultdict(list)
    for r in review:
        by_bio[r["biosample"]].append(r)
    candidate = []
    for bio, group in by_bio.items():
        group.sort(key=lambda r: (
            0 if r["evidence"] == "TAXON_CONFIRMED" else 1,
            -int(r["bytes_r1"]) - int(r["bytes_r2"]), r["run_accession"]))
        candidate.append(group[0])
    candidate.sort(key=lambda r: (r["evidence"], r["biosample"]))
    prior.write_tsv(out / "additional_cythera_wgs_candidates.tsv", candidate, ANNOTATED)

    audit = audit_previous(previous, rows)
    prior.write_tsv(out / "previous_race_taxonomy_audit.tsv", audit, AUDIT_COLS)
    prior.write_tsv(out / "retrieval_issues.tsv", errors,
                    ["source", "problem"])
    print("PREVIOUS\t" + str(dict(Counter(
        p["subspecies"] for p in previous))))
    print("PROJECT_TAXON_RECORDS\t" + str(len(dedup)))
    print("NEW_CYTHERA_PAIRED_WGS_BIOSAMPLES\t" + str(len(candidate)))
    for row in candidate[:12]:
        print("  " + row["biosample"] + "\t" + row["run_accession"] +
              "\t" + row["evidence"] + "\t" + row["source_query"])
    print("RACE_AUDIT\t" + str(dict(Counter(x["status"] for x in audit))))
    print("ENA_RETRIEVAL_ISSUES\t" + str(len(errors)))
    print("OUTPUT\t" + str(out))
    print("JOBS_SUBMITTED\t0")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("PROJECT_SWEEP_ERROR\t" + str(exc), file=sys.stderr)
        raise SystemExit(2)
