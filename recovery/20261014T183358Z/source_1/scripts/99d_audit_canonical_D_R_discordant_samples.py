#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

METADATA = (
    PROJECT
    / "data"
    / "derived"
    / "region_panel_metadata.tsv"
)


def clean_header(value):
    return " ".join(
        str(value).split()
    )


def binary(value):
    value = pd.to_numeric(
        value,
        errors="coerce",
    )

    if pd.isna(value):
        return None

    if float(value) == 0:
        return 0

    if float(value) == 1:
        return 1

    return None


def main():
    data = pd.read_csv(
        METADATA,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    data.columns = [
        clean_header(column)
        for column in data.columns
    ]

    required = {
        "Accession",
        "Taxon Name",
        "Data Type",
        "B",
        "R",
        "D",
    }

    missing = required.difference(data.columns)

    if missing:
        raise RuntimeError(
            "Missing columns: "
            + ", ".join(sorted(missing))
        )

    data["R_binary"] = data["R"].map(binary)
    data["D_binary"] = data["D"].map(binary)

    data["discordant_group"] = ""

    data.loc[
        (data["D_binary"] == 1)
        & (data["R_binary"] == 0),
        "discordant_group",
    ] = "D1R0"

    data.loc[
        (data["D_binary"] == 0)
        & (data["R_binary"] == 1),
        "discordant_group",
    ] = "D0R1"

    discordant = data[
        data["discordant_group"] != ""
    ].copy()

    discordant["piano_key"] = (
        discordant["Taxon Name"]
        .str.lower()
        .str.contains(
            "piano keys",
            regex=False,
        )
    )

    print("CANONICAL D/R-DISCORDANT METADATA")
    print()

    summary = (
        discordant
        .groupby(
            [
                "discordant_group",
                "Data Type",
                "piano_key",
            ],
            dropna=False,
        )
        .size()
        .reset_index(name="samples")
    )

    print(
        summary.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print()
    print("NON-PIANO-KEY DISCORDANT SAMPLES")
    print(
        "group\taccession\tdata_type\ttaxon"
    )

    retained = discordant[
        ~discordant["piano_key"]
    ].sort_values(
        [
            "discordant_group",
            "Taxon Name",
            "Accession",
        ]
    )

    for _, row in retained.iterrows():
        taxon = clean_header(
            row["Taxon Name"]
        )

        print(
            f'{row["discordant_group"]}\t'
            f'{row["Accession"]}\t'
            f'{row["Data Type"]}\t'
            f'{taxon}'
        )

    print()
    print("COUNTS USED FOR RECIPROCAL D/R TEST")

    print(
        retained[
            "discordant_group"
        ].value_counts().to_string()
    )


if __name__ == "__main__":
    main()
