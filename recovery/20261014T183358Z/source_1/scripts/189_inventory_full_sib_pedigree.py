#!/usr/bin/env python3

import csv
import io
import json
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path


OUTPUT_ROOT = Path(
    "analysis/genomewide_mixing_corrected_v1/"
    "full_sib_validation"
)
METADATA_DIR = OUTPUT_ROOT / "metadata"

PROJECTS = (
    "PRJEB11288",
    "PRJEB7581",
)

REFERENCE_BASES = 298_068_923
TARGET_COVERAGE = 15.0
OFFSPRING_COUNT = 12

EXCLUDED_NEXTERA_PROGENY = {
    11,
    16,
    17,
    18,
}

FIELDS = (
    "run_accession,"
    "experiment_accession,"
    "sample_accession,"
    "scientific_name,"
    "sample_title,"
    "experiment_title,"
    "library_strategy,"
    "library_layout,"
    "library_source,"
    "library_selection,"
    "instrument_platform,"
    "instrument_model,"
    "base_count,"
    "read_count,"
    "fastq_ftp,"
    "fastq_bytes"
)


def fetch_project(project):
    parameters = urllib.parse.urlencode(
        {
            "accession": project,
            "result": "read_run",
            "format": "tsv",
            "download": "false",
            "fields": FIELDS,
        }
    )

    url = (
        "https://www.ebi.ac.uk/ena/portal/api/"
        "filereport?"
        + parameters
    )

    last_error = None

    for attempt in range(1, 6):
        try:
            request = urllib.request.Request(
                url,
                headers={
                    "User-Agent":
                    "GynandroCalling-full-sib-audit/1.0"
                },
            )

            with urllib.request.urlopen(
                request,
                timeout=120,
            ) as handle:
                text = handle.read().decode("utf-8")

            rows = list(
                csv.DictReader(
                    io.StringIO(text),
                    delimiter="\t",
                )
            )

            if not rows:
                raise RuntimeError(
                    f"ENA returned no rows for {project}"
                )

            return text, rows

        except Exception as error:
            last_error = error

            if attempt < 5:
                time.sleep(5 * attempt)

    raise RuntimeError(
        f"ENA query failed for {project}: {last_error}"
    )


def integer(value):
    if value in {"", None}:
        return 0

    return int(value)


def split_integers(value):
    return [
        integer(piece)
        for piece in (value or "").split(";")
        if piece
    ]


def https_url(value):
    value = value.strip()

    if value.startswith(
        (
            "http://",
            "https://",
        )
    ):
        return value

    return "https://" + value


def progeny_number(title):
    match = re.search(
        r"\bprogeny\s+(\d+)\b",
        title,
        flags=re.I,
    )

    if match:
        return int(match.group(1))

    return None


def role_for(title):
    lowered = title.lower()

    if "f1 mother" in lowered:
        return "F1_mother"

    if "f1 father" in lowered:
        return "F1_father"

    if "f0 grandmother" in lowered:
        return "F0_grandmother"

    if progeny_number(title) is not None:
        return "F2_offspring"

    return "other"


def summarize_samples(rows):
    samples = defaultdict(
        lambda: {
            "titles": set(),
            "runs": [],
            "base_count": 0,
            "read_count": 0,
            "fastq_bytes": 0,
            "layouts": set(),
            "strategies": set(),
            "selections": set(),
            "models": set(),
            "mate1": [],
            "mate2": [],
        }
    )

    for row in rows:
        title = row["sample_title"]

        if (
            "melpomene x rosina cross"
            not in title.lower()
        ):
            continue

        sample = samples[
            row["sample_accession"]
        ]

        sample["titles"].add(title)
        sample["runs"].append(
            row["run_accession"]
        )
        sample["base_count"] += integer(
            row["base_count"]
        )
        sample["read_count"] += integer(
            row["read_count"]
        )
        sample["fastq_bytes"] += sum(
            split_integers(
                row["fastq_bytes"]
            )
        )
        sample["layouts"].add(
            row["library_layout"]
        )
        sample["strategies"].add(
            row["library_strategy"]
        )
        sample["selections"].add(
            row["library_selection"]
        )
        sample["models"].add(
            row["instrument_model"]
        )

        urls = [
            https_url(value)
            for value
            in row["fastq_ftp"].split(";")
            if value
        ]

        mate1 = [
            value
            for value in urls
            if re.search(
                r"_1\.fastq\.gz$",
                value,
            )
        ]

        mate2 = [
            value
            for value in urls
            if re.search(
                r"_2\.fastq\.gz$",
                value,
            )
        ]

        if (
            row["library_layout"] == "PAIRED"
            and (
                len(mate1) != 1
                or len(mate2) != 1
            )
        ):
            raise RuntimeError(
                "Expected one mate-1 and one "
                "mate-2 URL for "
                f"{row['run_accession']}; found "
                f"{len(mate1)} and {len(mate2)}"
            )

        sample["mate1"].extend(mate1)
        sample["mate2"].extend(mate2)

    result = []

    for accession, values in samples.items():
        if len(values["titles"]) != 1:
            raise RuntimeError(
                "Conflicting titles for "
                f"{accession}: "
                f"{values['titles']}"
            )

        title = next(iter(values["titles"]))
        number = progeny_number(title)

        result.append(
            {
                "sample_accession": accession,
                "sample_title": title,
                "role": role_for(title),
                "progeny_number": (
                    ""
                    if number is None
                    else number
                ),
                "run_accessions": ",".join(
                    sorted(values["runs"])
                ),
                "run_count": len(
                    values["runs"]
                ),
                "library_layouts": ",".join(
                    sorted(values["layouts"])
                ),
                "library_strategies": ",".join(
                    sorted(values["strategies"])
                ),
                "library_selections": ",".join(
                    sorted(values["selections"])
                ),
                "instrument_models": ",".join(
                    sorted(values["models"])
                ),
                "base_count":
                    values["base_count"],
                "read_count":
                    values["read_count"],
                "estimated_coverage": (
                    values["base_count"]
                    / REFERENCE_BASES
                ),
                "fastq_bytes":
                    values["fastq_bytes"],
                "mate1_urls": ",".join(
                    values["mate1"]
                ),
                "mate2_urls": ",".join(
                    values["mate2"]
                ),
            }
        )

    return sorted(
        result,
        key=lambda row: (
            row["role"],
            row["sample_accession"],
        ),
    )


def write_tsv(path, rows, fields):
    temporary = Path(
        str(path) + ".new"
    )

    with temporary.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fields,
        )
        writer.writeheader()
        writer.writerows(rows)

    temporary.replace(path)


def main():
    METADATA_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    project_rows = {}

    for project in PROJECTS:
        text, rows = fetch_project(project)

        raw_path = (
            METADATA_DIR
            / f"{project}.read_run.tsv"
        )
        temporary = Path(
            str(raw_path) + ".new"
        )
        temporary.write_text(text)
        temporary.replace(raw_path)

        project_rows[project] = rows

    samples = summarize_samples(
        project_rows["PRJEB11288"]
    )

    inventory_fields = list(
        samples[0]
    )

    write_tsv(
        METADATA_DIR
        / (
            "PRJEB11288."
            "melpomene_rosina.samples.tsv"
        ),
        samples,
        inventory_fields,
    )

    parents = [
        row
        for row in samples
        if row["role"]
        in {
            "F1_mother",
            "F1_father",
        }
    ]

    if {
        row["role"]
        for row in parents
    } != {
        "F1_mother",
        "F1_father",
    }:
        raise RuntimeError(
            "Could not identify exactly the "
            "F1 mother and F1 father"
        )

    eligible = []

    for row in samples:
        number = row["progeny_number"]

        if row["role"] != "F2_offspring":
            continue

        if integer(
            number
        ) in EXCLUDED_NEXTERA_PROGENY:
            continue

        if (
            row["library_layouts"]
            != "PAIRED"
        ):
            continue

        if (
            row["library_strategies"]
            != "WGS"
        ):
            continue

        if (
            not row["mate1_urls"]
            or not row["mate2_urls"]
        ):
            continue

        eligible.append(row)

    selected_children = sorted(
        eligible,
        key=lambda row: (
            abs(
                float(
                    row["estimated_coverage"]
                )
                - TARGET_COVERAGE
            ),
            integer(
                row["progeny_number"]
            ),
        ),
    )[:OFFSPRING_COUNT]

    if (
        len(selected_children)
        != OFFSPRING_COUNT
    ):
        raise RuntimeError(
            "Expected "
            f"{OFFSPRING_COUNT} eligible "
            "offspring; found "
            f"{len(selected_children)}"
        )

    role_order = {
        "F1_mother": 0,
        "F1_father": 1,
        "F2_offspring": 2,
    }

    selected = sorted(
        parents + selected_children,
        key=lambda row: (
            role_order[row["role"]],
            integer(
                row["progeny_number"]
            ),
        ),
    )

    manifest = []

    for array_id, row in enumerate(
        selected
    ):
        manifest.append(
            {
                "array_id": array_id,
                **dict(row),
            }
        )

    manifest_fields = (
        ["array_id"]
        + inventory_fields
    )

    manifest_path = (
        METADATA_DIR
        / (
            "initial_full_sib_"
            "download_manifest.tsv"
        )
    )

    write_tsv(
        manifest_path,
        manifest,
        manifest_fields,
    )

    first_samples = {
        row["sample_accession"]
        for row
        in project_rows["PRJEB11288"]
    }

    second_samples = {
        row["sample_accession"]
        for row
        in project_rows["PRJEB7581"]
    }

    first_runs = {
        row["run_accession"]
        for row
        in project_rows["PRJEB11288"]
    }

    second_runs = {
        row["run_accession"]
        for row
        in project_rows["PRJEB7581"]
    }

    audit = {
        "primary_project":
            "PRJEB11288",
        "secondary_project":
            "PRJEB7581",
        "PRJEB11288_runs":
            len(
                project_rows[
                    "PRJEB11288"
                ]
            ),
        "PRJEB11288_samples":
            len(first_samples),
        "PRJEB7581_runs":
            len(
                project_rows[
                    "PRJEB7581"
                ]
            ),
        "PRJEB7581_samples":
            len(second_samples),
        "overlapping_sample_accessions":
            len(
                first_samples
                & second_samples
            ),
        "overlapping_run_accessions":
            len(
                first_runs
                & second_runs
            ),
        "PRJEB7581_interpretation": (
            "additional runs from "
            "overlapping biological "
            "samples; not an independent "
            "pedigree"
        ),
        "selected_F1_parents":
            len(parents),
        "selected_F2_offspring":
            len(selected_children),
        "available_sibling_pairs": (
            len(selected_children)
            * (
                len(selected_children)
                - 1
            )
            // 2
        ),
        "estimated_fastq_GB": (
            sum(
                integer(
                    row["fastq_bytes"]
                )
                for row in selected
            )
            / 1e9
        ),
        "target_offspring_coverage":
            TARGET_COVERAGE,
        "excluded_Nextera_progeny":
            sorted(
                EXCLUDED_NEXTERA_PROGENY
            ),
        "reference_for_coverage":
            "GCF_056956725.1",
        "reference_bases":
            REFERENCE_BASES,
        "manifest":
            str(manifest_path),
    }

    audit_path = (
        METADATA_DIR
        / "full_sib_inventory.audit.json"
    )

    temporary = Path(
        str(audit_path) + ".new"
    )
    temporary.write_text(
        json.dumps(
            audit,
            indent=2,
        )
        + "\n"
    )
    temporary.replace(audit_path)

    print(
        "FULL-SIB PEDIGREE "
        "INVENTORY COMPLETE"
    )
    print(
        json.dumps(
            audit,
            indent=2,
        )
    )

    print("\nSELECTED SAMPLES")
    print(
        "array_id\t"
        "sample_accession\t"
        "role\t"
        "progeny_number\t"
        "runs\t"
        "coverage\t"
        "fastq_GB"
    )

    for row in manifest:
        print(
            f"{row['array_id']}\t"
            f"{row['sample_accession']}\t"
            f"{row['role']}\t"
            f"{row['progeny_number']}\t"
            f"{row['run_accessions']}\t"
            f"{float(row['estimated_coverage']):.2f}\t"
            f"{integer(row['fastq_bytes']) / 1e9:.3f}"
        )


if __name__ == "__main__":
    main()
