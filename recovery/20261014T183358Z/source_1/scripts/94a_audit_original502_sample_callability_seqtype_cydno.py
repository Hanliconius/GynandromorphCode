#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

CK = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
)

LONG = (
    CK
    / "optix_original502_selected_sites_genotype_plot.long.tsv"
)

SITES = (
    CK
    / "optix_original502_selected_sites_for_plotting.tsv"
)

OUT_SAMPLE_MODULE = (
    CK
    / "original502_sample_callability_by_module.tsv"
)

OUT_SAMPLE = (
    CK
    / "original502_sample_callability_summary.tsv"
)

METADATA_CANDIDATES = [
    PROJECT
    / "data"
    / "derived"
    / "region_panel_metadata.tsv",

    PROJECT
    / "external"
    / "krista_analysis"
    / "metadata.csv",

    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
    / "Optix_metadata.csv",
]


def read_table(path):

    if path.suffix.lower() == ".csv":
        return pd.read_csv(
            path,
            dtype=str,
            keep_default_na=False,
        )

    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )


def norm(x):

    if x is None:
        return ""

    x = str(x).strip()

    if x.lower() in {
        "",
        "na",
        "nan",
        "none",
        ".",
    }:
        return ""

    return x


def boolish(x):

    return str(x).strip().lower() in {
        "true",
        "t",
        "1",
        "yes",
        "y",
    }


def find_sample_column(df):

    preferred = [
        "sample",
        "Sample",
        "Accession",
        "accession",
        "Run",
        "run",
        "SRA",
        "sra",
        "ERR",
        "SRR",
    ]

    for col in preferred:
        if col in df.columns:
            return col

    return None


def candidate_sequence_columns(df):

    keys = (
        "seq",
        "strategy",
        "library",
        "data",
        "type",
        "source",
        "assay",
        "method",
        "experiment",
    )

    return [
        col
        for col in df.columns
        if any(
            key in col.lower()
            for key in keys
        )
    ]


def classify_sequence_text(text):

    text = str(text).lower()

    # Give explicit selective-sequencing labels precedence.
    if re.search(
        r"\bselective\b|targeted|target[ -]?capture|capture sequencing",
        text,
    ):
        return "selective_sequencing"

    if re.search(
        r"\bwgs\b|whole[ -]?genome|whole genome sequencing",
        text,
    ):
        return "WGS"

    return "unknown"


def recover_sequence_metadata(samples):

    result = pd.DataFrame(
        {
            "sample": sorted(samples),
            "sequencing_type": "unknown",
            "sequencing_metadata_source": "",
            "sequencing_metadata_evidence": "",
        }
    )

    result = result.set_index(
        "sample",
        drop=False,
    )

    print(
        "\n=================================================="
    )
    print(
        "SEQUENCING-TYPE METADATA SEARCH"
    )
    print(
        "=================================================="
    )

    for path in METADATA_CANDIDATES:

        if not path.exists():
            print(
                f"\nMetadata file absent: {path}"
            )
            continue

        df = read_table(
            path
        )

        sample_col = find_sample_column(
            df
        )

        print(
            f"\nMetadata file: {path}"
        )

        print(
            f"rows: {len(df)}"
        )

        print(
            f"sample column: {sample_col}"
        )

        seq_cols = candidate_sequence_columns(
            df
        )

        print(
            "candidate sequencing columns:"
        )

        if seq_cols:
            for col in seq_cols:
                print(
                    f"  {col}"
                )
        else:
            print(
                "  <none detected>"
            )

        if (
            sample_col is None
            or not seq_cols
        ):
            continue

        df["_sample"] = (
            df[sample_col]
            .map(norm)
        )

        df = df.loc[
            df["_sample"].isin(
                samples
            )
        ].copy()

        for _, row in df.iterrows():

            sample = row[
                "_sample"
            ]

            evidence_parts = []

            for col in seq_cols:

                value = norm(
                    row[col]
                )

                if value:
                    evidence_parts.append(
                        f"{col}={value}"
                    )

            if not evidence_parts:
                continue

            evidence = " | ".join(
                evidence_parts
            )

            seq_type = classify_sequence_text(
                evidence
            )

            # Only overwrite unknown with an actual classification.
            if (
                seq_type != "unknown"
                and result.loc[
                    sample,
                    "sequencing_type"
                ] == "unknown"
            ):

                result.loc[
                    sample,
                    "sequencing_type"
                ] = seq_type

                result.loc[
                    sample,
                    "sequencing_metadata_source"
                ] = str(
                    path
                )

                result.loc[
                    sample,
                    "sequencing_metadata_evidence"
                ] = evidence

    return result.reset_index(
        drop=True
    )


def main():

    long = pd.read_csv(
        LONG,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    sites = pd.read_csv(
        SITES,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    long["sample"] = long[
        "sample"
    ].map(norm)

    long["taxon_name"] = long[
        "taxon_name"
    ].map(norm)

    long["QC_pass"] = long[
        "QC_pass_DP2_GQ20"
    ].map(boolish)

    # -------------------------------------------------------
    # Only sites that really exist with matching REF/ALT in
    # the old regional VCFs contribute to sample missingness.
    #
    # A sample must not be penalised for a column that is
    # technically unavailable for every individual.
    # -------------------------------------------------------

    usable_sites = sites.loc[
        sites["source_status"]
        ==
        "exact_REF_ALT_match"
    ].copy()

    usable_keys = set(
        usable_sites[
            "site_key"
        ]
    )

    usable_long = long.loc[
        long["site_key"].isin(
            usable_keys
        )
    ].copy()

    expected_sites = (
        usable_sites.groupby(
            "module"
        )
        .size()
        .to_dict()
    )

    print(
        "\n=================================================="
    )
    print(
        "TECHNICALLY REPRESENTABLE FROZEN SITES"
    )
    print(
        "=================================================="
    )

    total_usable = 0

    for module in [
        "Band",
        "Ray2",
        "Dennis2",
    ]:

        n = int(
            expected_sites.get(
                module,
                0,
            )
        )

        total_usable += n

        print(
            f"{module:8s} {n}"
        )

    print(
        f"TOTAL    {total_usable}"
    )

    # -------------------------------------------------------
    # Per-sample, per-module callability.
    # -------------------------------------------------------

    records = []

    grouped = usable_long.groupby(
        [
            "module",
            "sample",
        ],
        sort=False,
    )

    for (
        module,
        sample,
    ), d in grouped:

        n_available = int(
            expected_sites[
                module
            ]
        )

        n_callable = int(
            d["QC_pass"].sum()
        )

        taxon_values = [
            x
            for x in d[
                "taxon_name"
            ].map(norm).unique()
            if x
        ]

        taxon = (
            taxon_values[0]
            if taxon_values
            else ""
        )

        records.append(
            {
                "sample":
                    sample,

                "taxon_name":
                    taxon,

                "module":
                    module,

                "available_selected_sites":
                    n_available,

                "callable_selected_sites":
                    n_callable,

                "missing_selected_sites":
                    n_available
                    - n_callable,

                "callable_fraction":
                    (
                        n_callable
                        /
                        n_available
                        if n_available
                        else np.nan
                    ),
            }
        )

    by_module = pd.DataFrame(
        records
    )

    by_module.to_csv(
        OUT_SAMPLE_MODULE,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # -------------------------------------------------------
    # Global callability across Band + Ray2 + Dennis2.
    # -------------------------------------------------------

    global_rows = []

    for sample, d in (
        by_module.groupby(
            "sample",
            sort=False,
        )
    ):

        n_available = int(
            d[
                "available_selected_sites"
            ].sum()
        )

        n_callable = int(
            d[
                "callable_selected_sites"
            ].sum()
        )

        taxon_values = [
            x
            for x in d[
                "taxon_name"
            ].map(norm).unique()
            if x
        ]

        taxon = (
            taxon_values[0]
            if taxon_values
            else ""
        )

        global_rows.append(
            {
                "sample":
                    sample,

                "taxon_name":
                    taxon,

                "available_selected_sites":
                    n_available,

                "callable_selected_sites":
                    n_callable,

                "missing_selected_sites":
                    n_available
                    - n_callable,

                "callable_fraction":
                    (
                        n_callable
                        /
                        n_available
                        if n_available
                        else np.nan
                    ),

                "all_selected_sites_missing":
                    n_callable == 0,

                "contains_cydno":
                    (
                        "cydno"
                        in taxon.lower()
                    ),
            }
        )

    sample_summary = pd.DataFrame(
        global_rows
    )

    seq_meta = recover_sequence_metadata(
        set(
            sample_summary[
                "sample"
            ]
        )
    )

    sample_summary = (
        sample_summary.merge(
            seq_meta,
            on="sample",
            how="left",
        )
    )

    sample_summary = (
        sample_summary.sort_values(
            [
                "callable_fraction",
                "taxon_name",
                "sample",
            ]
        )
        .reset_index(
            drop=True
        )
    )

    sample_summary.to_csv(
        OUT_SAMPLE,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # -------------------------------------------------------
    # Report all-zero individuals.
    # -------------------------------------------------------

    zero = sample_summary.loc[
        sample_summary[
            "all_selected_sites_missing"
        ]
    ].copy()

    print(
        "\n=================================================="
    )
    print(
        "INDIVIDUALS WITH 0 CALLABLE FROZEN SITES"
    )
    print(
        "=================================================="
    )

    print(
        f"N = {len(zero)}"
    )

    if len(zero):

        print(
            zero[
                [
                    "sample",
                    "taxon_name",
                    "sequencing_type",
                    "contains_cydno",
                ]
            ].to_string(
                index=False
            )
        )

    # -------------------------------------------------------
    # Missingness distribution.
    # -------------------------------------------------------

    print(
        "\n=================================================="
    )
    print(
        "GLOBAL CALLABILITY DISTRIBUTION"
    )
    print(
        "=================================================="
    )

    x = sample_summary[
        "callable_fraction"
    ].astype(float)

    print(
        x.describe(
            percentiles=[
                0.01,
                0.05,
                0.10,
                0.25,
                0.50,
                0.75,
                0.90,
                0.95,
                0.99,
            ]
        ).to_string()
    )

    print(
        "\nSample retention at candidate callability cutoffs:"
    )

    for cutoff in [
        0.00,
        0.10,
        0.25,
        0.50,
        0.75,
        0.90,
    ]:

        retained = int(
            (
                sample_summary[
                    "callable_fraction"
                ]
                >= cutoff
            ).sum()
        )

        retained_no_cydno = int(
            (
                (
                    sample_summary[
                        "callable_fraction"
                    ]
                    >= cutoff
                )
                &
                (
                    ~sample_summary[
                        "contains_cydno"
                    ]
                )
            ).sum()
        )

        print(
            f"  >= {cutoff:0.2f}: "
            f"{retained:3d} retained; "
            f"{retained_no_cydno:3d} "
            f"after cydno removal"
        )

    # -------------------------------------------------------
    # Cydno audit.
    # -------------------------------------------------------

    cydno = sample_summary.loc[
        sample_summary[
            "contains_cydno"
        ]
    ].copy()

    print(
        "\n=================================================="
    )
    print(
        "CYDNO TAXON EXCLUSION"
    )
    print(
        "=================================================="
    )

    print(
        f"N containing 'cydno' = {len(cydno)}"
    )

    if len(cydno):

        print(
            "\nCydno taxa:"
        )

        print(
            cydno[
                "taxon_name"
            ]
            .value_counts()
            .to_string()
        )

    # -------------------------------------------------------
    # Sequencing type.
    # -------------------------------------------------------

    print(
        "\n=================================================="
    )
    print(
        "SEQUENCING TYPE — FULL 481 COHORT"
    )
    print(
        "=================================================="
    )

    print(
        sample_summary[
            "sequencing_type"
        ]
        .value_counts(
            dropna=False
        )
        .to_string()
    )

    print(
        "\nSEQUENCING TYPE — AFTER CYDNO REMOVAL"
    )

    print(
        sample_summary.loc[
            ~sample_summary[
                "contains_cydno"
            ],
            "sequencing_type",
        ]
        .value_counts(
            dropna=False
        )
        .to_string()
    )

    print(
        "\nCALLABILITY BY SEQUENCING TYPE"
    )

    seq_stats = (
        sample_summary.groupby(
            "sequencing_type"
        )[
            "callable_fraction"
        ]
        .agg(
            [
                "count",
                "min",
                "median",
                "mean",
                "max",
            ]
        )
    )

    print(
        seq_stats.to_string()
    )

    print(
        "\n=================================================="
    )
    print(
        "20 LOWEST-CALLABILITY INDIVIDUALS"
    )
    print(
        "=================================================="
    )

    print(
        sample_summary[
            [
                "sample",
                "taxon_name",
                "sequencing_type",
                "callable_selected_sites",
                "available_selected_sites",
                "callable_fraction",
                "contains_cydno",
            ]
        ]
        .head(
            20
        )
        .to_string(
            index=False
        )
    )

    print(
        "\nOutputs:"
    )

    print(
        f"  {OUT_SAMPLE_MODULE}"
    )

    print(
        f"  {OUT_SAMPLE}"
    )


if __name__ == "__main__":
    main()
