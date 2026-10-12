#!/usr/bin/env python3
"""Inventory existing local inputs for a matched Hmel2.5/HmRos Yb rebuild.

Reads the existing genomewide manifest; never downloads or submits jobs.
Writes complete evidence to disk and prints a bounded focal-input summary.
"""
import argparse
import csv
import json
from pathlib import Path

def describe(value):
    p = Path(value).expanduser()
    absolute = str(p.resolve())
    return {"path": absolute, "exists": p.is_file(),
            "bytes": p.stat().st_size if p.is_file() else 0}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(
        "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"))
    args = parser.parse_args()
    root = args.project_root.resolve()
    meta = root / "analysis/genomewide_mixing_input/metadata"
    out = root / "analysis/yb_raw_from_scratch/local_input_inventory"
    out.mkdir(parents=True, exist_ok=True)
    manifest = meta / "resolved_download_manifest.tsv"
    rows = []
    errors = []
    if manifest.is_file():
        with manifest.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("role", "").lower() != "focal":
                    continue
                entry = {"library_id": row.get("library_id", ""),
                         "side": row.get("side", ""), "leg": row.get("leg", ""),
                         "source_accession": row.get("source_accession", ""),
                         "resolved_run": row.get("resolved_run", "")}
                for field in ("fastq1", "fastq2"):
                    value = row.get(field, "")
                    if value and Path(value).is_absolute():
                        entry[field] = describe(value)
                    else:
                        entry[field] = {"path": value, "exists": False, "bytes": 0}
                        errors.append(f"{entry['library_id']} {field}: missing absolute path")
                rows.append(entry)
    else:
        errors.append(f"Missing existing manifest: {manifest}")
    if len(rows) != 4:
        errors.append(f"Expected four focal libraries; found {len(rows)}")
    refs = {"HmRos": [], "Hmel2.5": []}
    pathfile = meta / "reference_fasta.path"
    if pathfile.is_file():
        value = pathfile.read_text().strip().splitlines()[0]
        if value and Path(value).is_absolute():
            refs["HmRos"].append(describe(value))
        else:
            errors.append("HmRos reference path is not absolute")
    else:
        errors.append(f"Missing HmRos reference path file: {pathfile}")
    candidates = [
        root / "analysis/yb_exact_coordinate_transfer/reference/hmel2_5.chromosomes.fa",
        root.parent / "Hmel2.5_PLUS_MITO_GENOME/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa.gz"]
    for p in candidates:
        refs["Hmel2.5"].append(describe(p))
    report = {"manifest": str(manifest), "focal_libraries": rows,
              "references": refs, "errors": errors,
              "required_references": ["Hmel2.5", "HmRos"],
              "processing": "matched raw-read alignment and coverage recipe on both references"}
    destination = out / "local_inputs.json"
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print("YB_LOCAL_INPUT_INVENTORY")
    for row in rows[:4]:
        print("LIBRARY\t" + row["library_id"] + "\t" + row["side"] + "\t" + row["leg"])
        for field in ("fastq1", "fastq2"):
            d = row[field]
            print(f"{field.upper()}\t{'OK' if d['exists'] and d['bytes'] else 'MISSING'}\t{d['path']}")
    for name, entries in refs.items():
        for d in entries:
            print(f"REFERENCE\t{name}\t{'OK' if d['exists'] and d['bytes'] else 'MISSING'}\t{d['path']}")
    for error in errors[:5]:
        print("CHECK\t" + error)
    print("REPORT\t" + str(destination))

if __name__ == "__main__":
    main()
