#!/usr/bin/env python3

from pathlib import Path
import re
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
CHECKPOINTS = PROJECT / "analysis/original_502/checkpoints"
SCRIPTS = PROJECT / "scripts"


def truthy(series):
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin(["true", "1", "yes", "t"])
    )


def clean(value):
    return " ".join(str(value).replace("\n", " ").split())


def main():
    print("PREVIOUS COHORT MANIFESTS")

    manifests = sorted(
        set(CHECKPOINTS.glob("*cohort*manifest*.tsv"))
        | set(CHECKPOINTS.glob("*filtered*manifest*.tsv"))
        | set(CHECKPOINTS.glob("*score_manifest*.tsv"))
    )

    for path in manifests:
        print(path)

    print()
    print("FILTERED COHORT MANIFEST AUDIT")

    target = (
        CHECKPOINTS
        / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
    )

    if not target.is_file():
        print(f"MISSING: {target}")
    else:
        table = pd.read_csv(
            target,
            sep="\t",
            low_memory=False,
        )

        print(f"Manifest: {target}")
        print(f"Rows: {len(table)}")
        print("Columns: " + ",".join(table.columns))

        if "keep" in table.columns:
            table["keep_bool"] = truthy(table["keep"])
        else:
            table["keep_bool"] = True

        if "excluded_taxon" in table.columns:
            table["excluded_taxon_bool"] = truthy(
                table["excluded_taxon"]
            )
        else:
            table["excluded_taxon_bool"] = False

        if "passes_callability" in table.columns:
            table["passes_callability_bool"] = truthy(
                table["passes_callability"]
            )
        else:
            table["passes_callability_bool"] = True

        sample_column = (
            "sample"
            if "sample" in table.columns
            else None
        )

        taxon_column = next(
            (
                column for column in table.columns
                if column.lower() in {
                    "taxon",
                    "taxon_name",
                    "taxon name",
                }
            ),
            None,
        )

        if sample_column and taxon_column:
            deduplicated = table.drop_duplicates(
                [sample_column, taxon_column]
            ).copy()

            print()
            print("UNIQUE-SAMPLE COUNTS BY FILTER STATUS")
            summary = (
                deduplicated.groupby(
                    [
                        "keep_bool",
                        "excluded_taxon_bool",
                        "passes_callability_bool",
                    ],
                    dropna=False,
                )[sample_column]
                .nunique()
                .reset_index(name="samples")
            )
            print(summary.to_csv(sep="\t", index=False).rstrip())

            print()
            print("TAXA MARKED AS EXCLUDED_TAXON")
            excluded = deduplicated[
                deduplicated["excluded_taxon_bool"]
            ]

            excluded_summary = (
                excluded.groupby(taxon_column, dropna=False)
                [sample_column]
                .nunique()
                .reset_index(name="samples")
                .sort_values(
                    ["samples", taxon_column],
                    ascending=[False, True],
                )
            )

            if excluded_summary.empty:
                print("None")
            else:
                print(
                    excluded_summary.to_csv(
                        sep="\t",
                        index=False,
                    ).rstrip()
                )

            print()
            print("NUMATA RECORDS")
            numata = deduplicated[
                deduplicated[taxon_column]
                .astype(str)
                .str.contains("numata", case=False, na=False)
            ]

            requested = [
                column for column in [
                    sample_column,
                    taxon_column,
                    "cohort",
                    "module",
                    "data_type",
                    "callable_fraction",
                    "excluded_taxon",
                    "passes_callability",
                    "keep",
                ]
                if column in numata.columns
            ]

            if numata.empty:
                print("None")
            else:
                print(
                    numata[requested]
                    .sort_values(sample_column)
                    .to_csv(
                        sep="\t",
                        index=False,
                    )
                    .rstrip()
                )

            print()
            print("ALL NON-KEPT TAXA")
            not_kept = deduplicated[
                ~deduplicated["keep_bool"]
            ]

            not_kept_summary = (
                not_kept.groupby(
                    [
                        taxon_column,
                        "excluded_taxon_bool",
                        "passes_callability_bool",
                    ],
                    dropna=False,
                )[sample_column]
                .nunique()
                .reset_index(name="samples")
                .sort_values(
                    ["samples", taxon_column],
                    ascending=[False, True],
                )
            )

            print(
                not_kept_summary.to_csv(
                    sep="\t",
                    index=False,
                ).rstrip()
            )

    print()
    print("FILTER LOGIC IN EXISTING SCANNERS")

    source_candidates = sorted(
        set(SCRIPTS.glob("94*.py"))
        | set(SCRIPTS.glob("95*.py"))
        | set(SCRIPTS.glob("96*.py"))
    )

    pattern = re.compile(
        r"numata|cydno|pachinus|elevatus|besckei|"
        r"excluded_taxon|passes_callability|"
        r"callable_fraction|piano|focal",
        re.IGNORECASE,
    )

    for path in source_candidates:
        matching = []

        for number, line in enumerate(
            path.read_text(errors="replace").splitlines(),
            start=1,
        ):
            if pattern.search(line):
                matching.append(
                    f"{number}\t{line.strip()}"
                )

        if matching:
            print()
            print(f"### {path.name}")
            for line in matching:
                print(line)

    print()
    print("AVAILABLE 96-SERIES SCRIPTS")

    for path in sorted(SCRIPTS.glob("96*")):
        print(path.name)


if __name__ == "__main__":
    main()
