#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

CK = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
)

INDIVIDUAL = (
    CK
    / "band_51pass_vs_502fail_individual_summary.tsv"
)

ERROR_LONG = (
    CK
    / "band_51pass_vs_502fail_misclassified.tsv"
)

COHORT = (
    CK
    / "optix_original_502_cohort_manifest.tsv"
)

META_FILES = {
    "Optix_metadata": (
        PROJECT
        / "external"
        / "krista_genotype_plots_original"
        / "Genotype_Plots"
        / "Optix_metadata.csv"
    ),
    "Krista_metadata": (
        PROJECT
        / "external"
        / "krista_analysis"
        / "metadata.csv"
    ),
    "region_panel_metadata": (
        PROJECT
        / "data"
        / "derived"
        / "region_panel_metadata.tsv"
    ),
    "figure2_metadata": (
        PROJECT
        / "analysis"
        / "checkpoints"
        / "nonfocal_figure2_samples.tsv"
    ),
}

OUT_MATCHES = (
    CK
    / "band_error_sample_metadata_source_rows.tsv"
)

OUT_COMPARE = (
    CK
    / "band_error_sample_metadata_comparison.tsv"
)


def read_any(path):
    sep = "," if path.suffix.lower() == ".csv" else "\t"

    return pd.read_csv(
        path,
        sep=sep,
        dtype=str,
        keep_default_na=False,
        engine="python",
    )


individual = read_any(
    INDIVIDUAL
)

errors = read_any(
    ERROR_LONG
)

cohort = read_any(
    COHORT
)


# ------------------------------------------------------------
# Worst offenders first.
# Require >= 3 errors so we focus on systematic discrepancies,
# rather than one-off site failures.
# ------------------------------------------------------------

individual["sites_misclassified_num"] = pd.to_numeric(
    individual["sites_misclassified"],
    errors="coerce",
)

individual["error_fraction_num"] = pd.to_numeric(
    individual["error_fraction"],
    errors="coerce",
)

targets = (
    individual.loc[
        individual["sites_misclassified_num"] >= 3
    ]
    .sort_values(
        [
            "sites_misclassified_num",
            "error_fraction_num",
        ],
        ascending=False,
    )
    .copy()
)

target_samples = list(
    targets["sample"]
    .astype(str)
    .unique()
)

print(
    f"\nSystematic-error individuals: "
    f"{len(target_samples)}"
)


# ------------------------------------------------------------
# Search every metadata row across every column for exact sample
# identifiers. This avoids assuming which identifier column is
# actually carrying ERR accessions.
# ------------------------------------------------------------

match_rows = []

for source_name, path in META_FILES.items():

    if not path.exists():
        print(
            f"WARNING missing metadata source: {path}"
        )
        continue

    df = read_any(
        path
    )

    print(
        f"\n{source_name}: "
        f"{len(df)} rows x {len(df.columns)} columns"
    )

    df_strings = df.astype(str)

    for sample in target_samples:

        mask = df_strings.apply(
            lambda col:
                col.str.strip().eq(sample)
        ).any(axis=1)

        hits = df.loc[
            mask
        ].copy()

        if hits.empty:
            continue

        for row_index, row in hits.iterrows():

            rec = {
                "sample":
                    sample,
                "metadata_source":
                    source_name,
                "metadata_file":
                    str(path),
                "source_row_index":
                    row_index,
            }

            for column in df.columns:
                rec[
                    f"source__{column}"
                ] = row[column]

            match_rows.append(
                rec
            )


matches = pd.DataFrame(
    match_rows
)

matches.to_csv(
    OUT_MATCHES,
    sep="\t",
    index=False,
    na_rep="NA",
)


# ------------------------------------------------------------
# Cohort values actually used in scoring.
# ------------------------------------------------------------

band_cohort = cohort.loc[
    (
        cohort["panel"].str.lower()
        == "band"
    )
    |
    (
        cohort["trait"].str.upper()
        == "B"
    )
].copy()

band_lookup = (
    band_cohort
    .drop_duplicates("sample")
    .set_index("sample")
)


# ------------------------------------------------------------
# One compact row per target sample, including all matched
# metadata rows collapsed as text so conflicts are visible.
# ------------------------------------------------------------

comparison_rows = []

for sample in target_samples:

    target = targets.loc[
        targets["sample"]
        == sample
    ].iloc[0]

    rec = {
        "sample":
            sample,

        "sites_callable":
            target.get(
                "sites_callable",
                "",
            ),

        "sites_misclassified":
            target.get(
                "sites_misclassified",
                "",
            ),

        "error_fraction":
            target.get(
                "error_fraction",
                "",
            ),

        "in_51sample_VCF":
            target.get(
                "in_51sample_VCF",
                "",
            ),
    }

    if sample in band_lookup.index:

        crow = band_lookup.loc[
            sample
        ]

        for column in [
            "phenotype",
            "reason",
            "taxon",
            "data_type",
        ]:
            if column in crow.index:
                rec[
                    f"scoring_manifest_{column}"
                ] = crow[column]

    sample_matches = matches.loc[
        matches["sample"]
        == sample
    ]

    rec[
        "metadata_sources_found"
    ] = ",".join(
        sorted(
            sample_matches[
                "metadata_source"
            ].unique()
        )
    )

    rec[
        "metadata_rows_found"
    ] = len(
        sample_matches
    )

    # Preserve every non-empty metadata value, column-labelled,
    # rather than choosing one silently.
    values = []

    for _, mrow in sample_matches.iterrows():

        source = mrow[
            "metadata_source"
        ]

        for column in sample_matches.columns:

            if not column.startswith(
                "source__"
            ):
                continue

            value = str(
                mrow[column]
            ).strip()

            if (
                value
                and value
                not in {
                    "NA",
                    "nan",
                }
            ):
                values.append(
                    f"{source}:"
                    f"{column.removeprefix('source__')}"
                    f"={value}"
                )

    rec[
        "all_matching_metadata"
    ] = " | ".join(
        dict.fromkeys(values)
    )

    comparison_rows.append(
        rec
    )


comparison = pd.DataFrame(
    comparison_rows
)

comparison.to_csv(
    OUT_COMPARE,
    sep="\t",
    index=False,
    na_rep="NA",
)


print(
    "\nMETADATA SOURCE RECOVERY"
)

if matches.empty:
    print(
        "No systematic-error samples found "
        "in secondary metadata files."
    )
else:
    print(
        matches.groupby(
            "metadata_source"
        )["sample"]
        .nunique()
        .to_string()
    )


print(
    "\nTOP 50 SYSTEMATIC ERRORS + METADATA"
)

show = comparison.head(
    50
)

columns = [
    "sample",
    "scoring_manifest_taxon",
    "scoring_manifest_phenotype",
    "sites_callable",
    "sites_misclassified",
    "in_51sample_VCF",
    "metadata_sources_found",
    "all_matching_metadata",
]

columns = [
    c
    for c in columns
    if c in show.columns
]

print(
    show[
        columns
    ].to_string(
        index=False
    )
)


print(
    "\nOutputs:"
)

print(
    f"  {OUT_MATCHES}"
)

print(
    f"  {OUT_COMPARE}"
)
