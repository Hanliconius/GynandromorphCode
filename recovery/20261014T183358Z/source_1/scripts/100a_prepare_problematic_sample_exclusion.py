#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

INPUT_METADATA = (
    PROJECT / "data/derived/region_panel_metadata.tsv"
)

OUTPUT_METADATA = (
    PROJECT
    / "data/derived"
    / "region_panel_metadata.drop_problematic.tsv"
)

EXCLUSION_MANIFEST = (
    PROJECT
    / "analysis/drop_problematic_reanalysis/checkpoints"
    / "excluded_comparative_samples.tsv"
)


def clean(value):
    return " ".join(str(value).replace("\n", " ").split())


def main():
    EXCLUSION_MANIFEST.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    metadata = pd.read_csv(
        INPUT_METADATA,
        sep="\t",
        low_memory=False,
    )

    original_columns = list(metadata.columns)
    normalized = {
        clean(column).lower(): column
        for column in metadata.columns
    }

    accession_column = normalized["accession"]
    taxon_column = normalized["taxon name"]
    B_column = normalized["b"]
    R_column = normalized["r"]
    D_column = normalized["d"]

    B = pd.to_numeric(metadata[B_column], errors="coerce")
    R = pd.to_numeric(metadata[R_column], errors="coerce")
    D = pd.to_numeric(metadata[D_column], errors="coerce")

    contigua_D0R0 = metadata[
        metadata[taxon_column]
        .astype(str)
        .str.contains("contigua", case=False, na=False)
        & B.eq(0)
        & D.eq(0)
        & R.eq(0)
    ][accession_column].dropna().astype(str).tolist()

    if len(contigua_D0R0) != 1:
        raise RuntimeError(
            "Expected exactly one contigua B0D0R0 specimen; "
            f"found {len(contigua_D0R0)}: {contigua_D0R0}"
        )

    explicit = [
        "ERR3653118",
        "ERR3653120",
        "ERR3653298",
        "ERR3653306",
        "ERR3653309",
    ]

    excluded = explicit + contigua_D0R0

    if len(set(excluded)) != 6:
        raise RuntimeError(
            f"Expected six unique exclusions; recovered: {excluded}"
        )

    metadata_accessions = set(
        metadata[accession_column].dropna().astype(str)
    )

    missing = [
        sample for sample in excluded
        if sample not in metadata_accessions
    ]

    if missing:
        raise RuntimeError(
            "Requested exclusions absent from metadata: "
            + ",".join(missing)
        )

    manifest_rows = []

    for sample in explicit:
        row = metadata[
            metadata[accession_column].astype(str).eq(sample)
        ].iloc[0]

        manifest_rows.append({
            "sample": sample,
            "reason": "user_requested_comparative_outlier",
            "taxon": clean(row[taxon_column]),
            "B": row[B_column],
            "R": row[R_column],
            "D": row[D_column],
        })

    for sample in contigua_D0R0:
        row = metadata[
            metadata[accession_column].astype(str).eq(sample)
        ].iloc[0]

        manifest_rows.append({
            "sample": sample,
            "reason": "contigua_not_D0R1_B0D0R0",
            "taxon": clean(row[taxon_column]),
            "B": row[B_column],
            "R": row[R_column],
            "D": row[D_column],
        })

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(
        EXCLUSION_MANIFEST,
        sep="\t",
        index=False,
    )

    filtered = metadata[
        ~metadata[accession_column].astype(str).isin(excluded)
    ].copy()

    filtered = filtered[original_columns]

    filtered.to_csv(
        OUTPUT_METADATA,
        sep="\t",
        index=False,
    )

    print("PROBLEMATIC COMPARATIVE SAMPLE EXCLUSION")
    print()
    print(manifest.to_csv(sep="\t", index=False).rstrip())
    print()
    print(f"Original metadata rows: {len(metadata)}")
    print(f"Excluded rows: {len(manifest)}")
    print(f"Filtered metadata rows: {len(filtered)}")
    print()
    print(f"Filtered metadata: {OUTPUT_METADATA}")
    print(f"Exclusion manifest: {EXCLUSION_MANIFEST}")


if __name__ == "__main__":
    main()
