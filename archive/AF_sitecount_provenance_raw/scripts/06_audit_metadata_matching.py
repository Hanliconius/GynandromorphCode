#!/usr/bin/env python3

import csv
import hashlib
from collections import defaultdict
from pathlib import Path

from openpyxl import load_workbook
from openpyxl.utils import get_column_letter


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
CHECKPOINTS = PROJECT / "analysis" / "checkpoints"
SAMPLE_DIR = CHECKPOINTS / "vcf_sample_lists"
WORKBOOK = PROJECT / "metadata" / "Wing_Pattern_MetaData.xlsx"

FOCAL = {
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
}


def read_samples(path):
    with path.open() as handle:
        return [line.strip() for line in handle if line.strip()]


def sample_hash(samples):
    content = "\n".join(samples).encode("utf-8")
    return hashlib.sha256(content).hexdigest()


def clean_value(value):
    if value is None:
        return ""
    return str(value).strip()


def main():
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)

    sample_files = sorted(
        SAMPLE_DIR.glob("data_raw_region_panel_*.samples.txt")
    )

    if not sample_files:
        print("ERROR: region-panel sample lists were not found.")
        return 1

    sample_sets_output = CHECKPOINTS / "region_panel_sample_sets.tsv"

    canonical_samples = read_samples(sample_files[0])
    canonical_set = set(canonical_samples)
    canonical_hash = sample_hash(canonical_samples)

    with sample_sets_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sample_file",
                "samples",
                "unique_samples",
                "order_sha256",
                "identical_to_first",
            ]
        )

        for path in sample_files:
            samples = read_samples(path)
            writer.writerow(
                [
                    path.name,
                    len(samples),
                    len(set(samples)),
                    sample_hash(samples),
                    int(
                        samples == canonical_samples
                        and sample_hash(samples) == canonical_hash
                    ),
                ]
            )

    workbook = load_workbook(
        WORKBOOK,
        read_only=True,
        data_only=True,
    )

    sheet_output = CHECKPOINTS / "metadata_workbook_sheets.tsv"
    match_output = CHECKPOINTS / "metadata_sample_column_matches.tsv"

    sheet_rows = []
    match_rows = []

    for worksheet in workbook.worksheets:
        headers = {}
        matched_by_column = defaultdict(set)
        cells_matching_by_column = defaultdict(int)

        for row_number, row in enumerate(
            worksheet.iter_rows(values_only=True),
            start=1,
        ):
            for column_number, value in enumerate(row, start=1):
                cleaned = clean_value(value)

                if (
                    row_number <= 10
                    and column_number not in headers
                    and cleaned
                ):
                    headers[column_number] = cleaned

                if cleaned in canonical_set:
                    matched_by_column[column_number].add(cleaned)
                    cells_matching_by_column[column_number] += 1

        sheet_rows.append(
            [
                worksheet.title,
                worksheet.max_row,
                worksheet.max_column,
                sum(len(values) for values in matched_by_column.values()),
                len(matched_by_column),
            ]
        )

        for column_number, matched_samples in matched_by_column.items():
            match_rows.append(
                [
                    worksheet.title,
                    get_column_letter(column_number),
                    headers.get(column_number, ""),
                    len(matched_samples),
                    cells_matching_by_column[column_number],
                    len(canonical_set - matched_samples),
                ]
            )

    workbook.close()

    with sheet_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sheet",
                "rows",
                "columns",
                "unique_sample_matches_across_columns",
                "columns_with_matches",
            ]
        )
        writer.writerows(sheet_rows)

    match_rows.sort(key=lambda row: (-row[3], row[0], row[1]))

    with match_output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "sheet",
                "column",
                "header_candidate",
                "unique_vcf_samples_matched",
                "matching_cells",
                "vcf_samples_not_matched",
            ]
        )
        writer.writerows(match_rows)

    focal_present = sorted(canonical_set & FOCAL)
    focal_missing = sorted(FOCAL - canonical_set)

    identical_files = 0
    for path in sample_files:
        if read_samples(path) == canonical_samples:
            identical_files += 1

    print("REGION-PANEL SAMPLE CONSISTENCY")
    print(f"VCF sample lists: {len(sample_files)}")
    print(f"Lists identical and in the same order: {identical_files}")
    print(f"Samples per panel: {len(canonical_samples)}")
    print(f"Unique samples: {len(canonical_set)}")
    print(f"Focal samples present: {','.join(focal_present)}")
    print(
        "Focal samples missing: "
        + (",".join(focal_missing) if focal_missing else "none")
    )
    print(f"Non-focal comparative samples: {len(canonical_set - FOCAL)}")

    print("\nWORKBOOK SHEETS")
    print(
        "sheet\trows\tcolumns\t"
        "sample_matches_across_columns\tcolumns_with_matches"
    )
    for row in sheet_rows:
        print("\t".join(str(value) for value in row))

    print("\nBEST SAMPLE-COLUMN MATCHES")
    print(
        "sheet\tcolumn\theader_candidate\t"
        "unique_matches\tmatching_cells\tunmatched_vcf_samples"
    )
    for row in match_rows[:15]:
        print("\t".join(str(value) for value in row))

    print(f"\nFull sample-set audit: {sample_sets_output}")
    print(f"Full workbook summary: {sheet_output}")
    print(f"Full column-match audit: {match_output}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
