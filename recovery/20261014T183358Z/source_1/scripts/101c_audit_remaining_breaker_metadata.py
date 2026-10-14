#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


project = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

metadata_file = (
    project
    / "data"
    / "derived"
    / "region_panel_metadata.tsv"
)

manifest_file = (
    project
    / "analysis"
    / "drop_problematic_reanalysis"
    / "checkpoints"
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)

samples = [
    "ERR3650481",
    "ERR3652402",
]


def normalized(value):

    return (
        str(value)
        .strip()
        .lower()
        .replace(" ", "")
        .replace("_", "")
        .replace("-", "")
    )


def find_column(columns, target):

    target = normalized(target)

    for column in columns:
        if normalized(column) == target:
            return column

    raise RuntimeError(
        f"Column not found: {target}; "
        f"available columns are {list(columns)}"
    )


metadata = pd.read_csv(
    metadata_file,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

accession_column = find_column(
    metadata.columns,
    "Accession",
)

wanted_columns = [
    accession_column,
    find_column(metadata.columns, "Taxon Name"),
    find_column(metadata.columns, "Country"),
    find_column(metadata.columns, "Data Type"),
    find_column(metadata.columns, "B"),
    find_column(metadata.columns, "D"),
    find_column(metadata.columns, "R"),
    find_column(metadata.columns, "Optix Figure"),
]

print("ALL LOGICAL METADATA RECORDS")

for sample in samples:

    subset = metadata.loc[
        metadata[accession_column].astype(str).str.strip()
        == sample,
        wanted_columns,
    ].copy()

    print()
    print(sample)
    print(f"Records: {len(subset)}")

    if len(subset):
        print(
            subset.to_string(
                index=False
            )
        )
    else:
        print("No metadata record found.")

manifest = pd.read_csv(
    manifest_file,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

print()
print("CORRECTED FILTER-MANIFEST RECORDS")

subset = manifest.loc[
    manifest["sample"].isin(samples),
    [
        "cohort",
        "module",
        "sample",
        "taxon_name",
        "data_type",
        "callable_fraction",
        "excluded_taxon",
        "passes_callability",
        "keep",
    ],
].copy()

print(
    subset.to_string(
        index=False
    )
)

print()
print("DUPLICATE-ACCESSION AUDIT")

counts = (
    metadata[accession_column]
    .astype(str)
    .str.strip()
    .value_counts()
)

for sample in samples:
    print(
        f"{sample}\tmetadata_records="
        f"{int(counts.get(sample, 0))}"
    )
