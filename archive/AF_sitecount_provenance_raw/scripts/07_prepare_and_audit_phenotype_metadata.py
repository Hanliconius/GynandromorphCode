#!/usr/bin/env python3

import csv
import re
from collections import Counter
from pathlib import Path

from openpyxl import load_workbook
from openpyxl.utils import get_column_letter


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
WORKBOOK = PROJECT / "metadata" / "Wing_Pattern_MetaData.xlsx"

SAMPLE_LIST = (
    PROJECT
    / "analysis/checkpoints/vcf_sample_lists"
    / "data_raw_region_panel_optix_filtered.vcf.gz.samples.txt"
)

DERIVED_DIR = PROJECT / "data/derived"
CHECKPOINTS = PROJECT / "analysis/checkpoints"

OUTPUT_METADATA = DERIVED_DIR / "region_panel_metadata.tsv"
OUTPUT_COUNTS = CHECKPOINTS / "metadata_candidate_column_counts.tsv"

FOCAL = {
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
}


def clean(value):
    if value is None:
        return ""
    return str(value).strip()


def normalized_header(value):
    value = clean(value).lower()
    value = re.sub(r"[^a-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def display_value(value):
    value = clean(value)
    return value if value else "<missing>"


def candidate_column(header):
    normalized = normalized_header(header)

    exact = {
        "b",
        "d",
        "r",
        "taxon",
        "species",
        "subspecies",
        "race",
        "country",
        "data type",
        "figure",
        "optix",
        "ivory",
        "yb",
        "yb dorsal",
        "yb ventral",
    }

    return (
        normalized in exact
        or normalized.startswith("yb ")
        or "phenotype" in normalized
        or "figure" in normalized
        or normalized.startswith("optix ")
        or normalized.startswith("ivory ")
    )


def main():
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)

    if not WORKBOOK.is_file():
        print(f"ERROR: missing workbook: {WORKBOOK}")
        return 1

    if not SAMPLE_LIST.is_file():
        print(f"ERROR: missing VCF sample list: {SAMPLE_LIST}")
        return 1

    with SAMPLE_LIST.open() as handle:
        vcf_samples = [
            line.strip()
            for line in handle
            if line.strip()
        ]

    vcf_sample_set = set(vcf_samples)

    workbook = load_workbook(
        WORKBOOK,
        read_only=True,
        data_only=True,
    )

    if "Metadata" not in workbook.sheetnames:
        print("ERROR: workbook does not contain a Metadata worksheet.")
        workbook.close()
        return 1

    worksheet = workbook["Metadata"]
    row_iterator = worksheet.iter_rows(values_only=True)

    try:
        raw_headers = next(row_iterator)
    except StopIteration:
        print("ERROR: Metadata worksheet is empty.")
        workbook.close()
        return 1

    headers = [
        clean(value) if clean(value) else f"unnamed_{index}"
        for index, value in enumerate(raw_headers, start=1)
    ]

    normalized_headers = [
        normalized_header(header)
        for header in headers
    ]

    accession_indexes = [
        index
        for index, header in enumerate(normalized_headers)
        if header == "accession"
    ]

    if len(accession_indexes) != 1:
        print(
            "ERROR: expected exactly one Accession column; "
            f"found {len(accession_indexes)}."
        )
        workbook.close()
        return 1

    accession_index = accession_indexes[0]

    metadata_by_accession = {}
    duplicate_accessions = set()
    rows_scanned = 0
    rows_with_accessions = 0

    for row in row_iterator:
        rows_scanned += 1

        padded = list(row) + [""] * max(0, len(headers) - len(row))
        padded = padded[: len(headers)]

        accession = clean(padded[accession_index])

        if not accession:
            continue

        rows_with_accessions += 1

        if accession in metadata_by_accession:
            duplicate_accessions.add(accession)
        else:
            metadata_by_accession[accession] = [
                clean(value)
                for value in padded
            ]

    workbook.close()

    if duplicate_accessions:
        print(
            "ERROR: duplicate accessions detected in metadata: "
            f"{len(duplicate_accessions)}"
        )
        return 1

    missing_metadata = [
        sample
        for sample in vcf_samples
        if sample not in metadata_by_accession
    ]

    if missing_metadata:
        print(
            "ERROR: VCF samples without metadata: "
            f"{len(missing_metadata)}"
        )
        return 1

    with OUTPUT_METADATA.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)

        for sample in vcf_samples:
            writer.writerow(metadata_by_accession[sample])

    selected_columns = [
        index
        for index, header in enumerate(headers)
        if candidate_column(header)
    ]

    count_rows = []

    for column_index in selected_columns:
        header = headers[column_index]

        all_counts = Counter()
        comparative_counts = Counter()

        for sample in vcf_samples:
            value = display_value(
                metadata_by_accession[sample][column_index]
            )
            all_counts[value] += 1

            if sample not in FOCAL:
                comparative_counts[value] += 1

        all_unique = len(all_counts)
        comparative_unique = len(comparative_counts)

        for value, count in sorted(
            comparative_counts.items(),
            key=lambda item: (-item[1], item[0]),
        ):
            count_rows.append(
                [
                    get_column_letter(column_index + 1),
                    header,
                    value,
                    count,
                    all_counts.get(value, 0),
                    comparative_unique,
                    all_unique,
                ]
            )

    with OUTPUT_COUNTS.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "column",
                "header",
                "value",
                "comparative_n",
                "all_502_n",
                "comparative_unique_values",
                "all_502_unique_values",
            ]
        )
        writer.writerows(count_rows)

    print("METADATA STRUCTURE")
    print(f"Worksheet: Metadata")
    print(f"Rows scanned after header: {rows_scanned}")
    print(f"Rows with accessions: {rows_with_accessions}")
    print(f"VCF samples matched: {len(vcf_samples)}")
    print(f"Comparative samples: {len(vcf_sample_set - FOCAL)}")
    print(f"Focal samples: {len(vcf_sample_set & FOCAL)}")
    print(f"Columns: {len(headers)}")

    print("\nCOLUMN HEADERS")
    print("column\theader")
    for index, header in enumerate(headers, start=1):
        print(f"{get_column_letter(index)}\t{header}")

    print("\nCANDIDATE PHENOTYPE/METADATA COUNTS")
    print(
        "column\theader\tvalue\tcomparative_n\tall_502_n"
    )

    for row in count_rows:
        print(
            f"{row[0]}\t{row[1]}\t{row[2]}\t"
            f"{row[3]}\t{row[4]}"
        )

    print(f"\nJoined metadata table: {OUTPUT_METADATA}")
    print(f"Candidate-column counts: {OUTPUT_COUNTS}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
