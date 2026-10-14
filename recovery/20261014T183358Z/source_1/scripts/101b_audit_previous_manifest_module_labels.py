#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

MANIFEST = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)


def truthy(value):
    return str(value).strip().lower() in {
        "true", "1", "yes", "t"
    }


def main():
    if not MANIFEST.is_file():
        raise FileNotFoundError(
            f"Manifest is missing: {MANIFEST}"
        )

    table = pd.read_csv(
        MANIFEST,
        sep="\t",
        low_memory=False,
    )

    required = {
        "cohort",
        "module",
        "sample",
        "keep",
    }

    missing = sorted(required - set(table.columns))

    if missing:
        raise RuntimeError(
            "Manifest is missing required columns: "
            + ",".join(missing)
        )

    table["keep_bool"] = table["keep"].map(truthy)

    print("PREVIOUS FILTERED-COHORT MANIFEST")
    print(f"File: {MANIFEST}")
    print(f"Rows: {len(table)}")

    print()
    print("ALL MODULE LABELS")
    module_counts = (
        table.groupby("module", dropna=False)
        .size()
        .reset_index(name="rows")
        .sort_values("module")
    )
    print(
        module_counts.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print()
    print("UNIQUE SAMPLES BY COHORT AND MODULE")
    cohort_counts = (
        table.groupby(
            ["cohort", "module"],
            dropna=False,
        )["sample"]
        .nunique()
        .reset_index(name="samples")
        .sort_values(["cohort", "module"])
    )
    print(
        cohort_counts.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print()
    print("KEPT SAMPLES BY COHORT AND MODULE")
    kept_counts = (
        table[table["keep_bool"]]
        .groupby(
            ["cohort", "module"],
            dropna=False,
        )["sample"]
        .nunique()
        .reset_index(name="kept_samples")
        .sort_values(["cohort", "module"])
    )
    print(
        kept_counts.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print()
    print("MODULE LABELS CONTAINING B, D, OR R")

    mask = (
        table["module"]
        .astype(str)
        .str.contains(
            "band|dennis|ray|^b$|^d$|^r$",
            case=False,
            regex=True,
            na=False,
        )
    )

    matching = (
        table.loc[mask, ["module"]]
        .drop_duplicates()
        .sort_values("module")
    )

    if matching.empty:
        print("None")
    else:
        print(
            matching.to_csv(
                sep="\t",
                index=False,
            ).rstrip()
        )


if __name__ == "__main__":
    main()
