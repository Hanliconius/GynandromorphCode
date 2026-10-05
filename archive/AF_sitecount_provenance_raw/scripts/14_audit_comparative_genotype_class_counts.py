#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
GENOTYPE_DIR = PROJECT / "data/derived/region_genotypes"
METADATA_FILE = PROJECT / "data/derived/region_panel_metadata.tsv"
CHECKPOINT_DIR = PROJECT / "analysis/checkpoints"

FOCAL = {
    "TN2L_S1", "TN2R_S2", "TN3L_S3",
    "TN3R_S4", "TN_Left", "TN_Right",
}

PANELS = {
    "optix": {
        "figure_column": "Optix Figure",
    },
    "ivory": {
        "figure_column": "Ivory Figure",
    },
}


def normalize_name(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def resolve_column(dataframe, requested):
    target = normalize_name(requested)
    matches = [
        column for column in dataframe.columns
        if normalize_name(column) == target
    ]

    if len(matches) != 1:
        raise ValueError(
            f"Could not resolve {requested!r}; matches={matches}"
        )

    return matches[0]


def numeric(series):
    return pd.to_numeric(
        series.replace({"": np.nan, ".": np.nan}),
        errors="coerce",
    )


def dosage(value):
    value = str(value)

    if value in {"", ".", "./.", ".|.", "nan"}:
        return np.nan

    fields = re.split(r"[/|]", value)

    if len(fields) != 2 or any(field == "." for field in fields):
        return np.nan

    try:
        alleles = [int(field) for field in fields]
    except ValueError:
        return np.nan

    if any(allele not in {0, 1} for allele in alleles):
        return np.nan

    return float(sum(alleles))


def quantile(series, probability):
    series = pd.to_numeric(series, errors="coerce").dropna()

    if series.empty:
        return np.nan

    return float(series.quantile(probability))


def main():
    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    accession_column = resolve_column(metadata, "Accession")
    taxon_column = resolve_column(metadata, "Taxon Name")

    metadata["sample"] = metadata[accession_column].astype(str)
    metadata["taxon_qc"] = metadata[taxon_column].astype(str)
    metadata["role_qc"] = "comparative"

    metadata.loc[
        metadata["taxon_qc"].str.contains(
            "piano keys",
            case=False,
            na=False,
        ),
        "role_qc",
    ] = "piano_keys"

    metadata.loc[
        metadata["sample"].isin(FOCAL),
        "role_qc",
    ] = "focal"

    metadata_by_sample = metadata.set_index("sample")

    site_rows = []

    for panel, configuration in PANELS.items():
        figure_column = resolve_column(
            metadata,
            configuration["figure_column"],
        )

        figure_by_sample = numeric(
            metadata.set_index("sample")[figure_column]
        )

        calls = pd.read_csv(
            GENOTYPE_DIR / f"{panel}.genotypes.tsv.gz",
            sep="\t",
            dtype=str,
            compression="gzip",
            keep_default_na=False,
        )

        calls["role"] = calls["sample"].map(
            metadata_by_sample["role_qc"]
        )
        calls["figure"] = calls["sample"].map(
            figure_by_sample
        )
        calls["dosage"] = calls["GT"].map(dosage)
        calls["DP_numeric"] = numeric(calls["DP"])
        calls["GQ_numeric"] = numeric(calls["GQ"])

        calls["passes_qc"] = (
            calls["dosage"].notna()
            & calls["DP_numeric"].ge(2)
            & calls["GQ_numeric"].ge(20)
        )

        calls["biallelic_SNP"] = (
            ~calls["alt"].str.contains(",", regex=False)
            & calls["ref"].str.len().eq(1)
            & calls["alt"].str.len().eq(1)
        )

        cohort = calls.loc[
            calls["role"].eq("comparative")
            & calls["figure"].eq(1)
            & calls["biallelic_SNP"]
        ].copy()

        for site_key, site in cohort.groupby(
            ["chrom", "pos", "ref", "alt"],
            sort=False,
        ):
            passing = site.loc[site["passes_qc"]]
            counts = passing["dosage"].value_counts()

            ref_hom = int(counts.get(0.0, 0))
            het = int(counts.get(1.0, 0))
            alt_hom = int(counts.get(2.0, 0))
            called = ref_hom + het + alt_hom

            alt_alleles = het + (2 * alt_hom)
            allele_number = 2 * called
            alt_frequency = (
                alt_alleles / allele_number
                if allele_number else np.nan
            )
            minor_allele_frequency = (
                min(alt_frequency, 1 - alt_frequency)
                if np.isfinite(alt_frequency)
                else np.nan
            )

            alt_dominant_min_class = min(
                het + alt_hom,
                ref_hom,
            )
            ref_dominant_min_class = min(
                ref_hom + het,
                alt_hom,
            )
            alt_recessive_min_class = min(
                alt_hom,
                ref_hom + het,
            )
            ref_recessive_min_class = min(
                ref_hom,
                het + alt_hom,
            )

            site_rows.append({
                "panel": panel,
                "chrom": site_key[0],
                "pos": int(site_key[1]),
                "ref": site_key[2],
                "alt": site_key[3],
                "cohort_n": int(site["sample"].nunique()),
                "called_n": called,
                "call_rate": (
                    called / len(site) if len(site) else np.nan
                ),
                "ref_hom_n": ref_hom,
                "het_n": het,
                "alt_hom_n": alt_hom,
                "alt_allele_frequency": alt_frequency,
                "minor_allele_frequency": minor_allele_frequency,
                "best_dominant_min_class": max(
                    alt_dominant_min_class,
                    ref_dominant_min_class,
                ),
                "best_recessive_min_class": max(
                    alt_recessive_min_class,
                    ref_recessive_min_class,
                ),
            })

    sites = pd.DataFrame(site_rows)

    output_file = (
        CHECKPOINT_DIR
        / "primary_cohort_genotype_class_counts.tsv"
    )
    sites.to_csv(
        output_file,
        sep="\t",
        index=False,
        float_format="%.8g",
    )

    summaries = []

    for panel, subset in sites.groupby("panel", sort=False):
        summaries.append({
            "panel": panel,
            "sites": len(subset),
            "polymorphic_sites": int(
                subset["minor_allele_frequency"].gt(0).sum()
            ),
            "MAF_p10": quantile(
                subset["minor_allele_frequency"], 0.10
            ),
            "MAF_median": quantile(
                subset["minor_allele_frequency"], 0.50
            ),
            "MAF_p90": quantile(
                subset["minor_allele_frequency"], 0.90
            ),
            "MAF_max": subset["minor_allele_frequency"].max(),
            "sites_dominant_min_class_ge_1": int(
                subset["best_dominant_min_class"].ge(1).sum()
            ),
            "sites_dominant_min_class_ge_5": int(
                subset["best_dominant_min_class"].ge(5).sum()
            ),
            "sites_dominant_min_class_ge_10": int(
                subset["best_dominant_min_class"].ge(10).sum()
            ),
            "sites_recessive_min_class_ge_1": int(
                subset["best_recessive_min_class"].ge(1).sum()
            ),
            "sites_recessive_min_class_ge_5": int(
                subset["best_recessive_min_class"].ge(5).sum()
            ),
            "sites_recessive_min_class_ge_10": int(
                subset["best_recessive_min_class"].ge(10).sum()
            ),
        })

    summary = pd.DataFrame(summaries)

    most_variable = (
        sites.sort_values(
            [
                "panel",
                "minor_allele_frequency",
                "called_n",
            ],
            ascending=[True, False, False],
        )
        .groupby("panel", sort=False)
        .head(12)
    )

    display_columns = [
        "panel", "chrom", "pos", "ref", "alt",
        "called_n", "ref_hom_n", "het_n", "alt_hom_n",
        "minor_allele_frequency",
        "best_dominant_min_class",
        "best_recessive_min_class",
    ]

    print("GENOTYPE-CLASS SUMMARY")
    print(
        summary.to_csv(
            sep="\t",
            index=False,
            float_format="%.5g",
        ).rstrip()
    )

    print("\nMOST VARIABLE SITES IN EACH PRIMARY COHORT")
    print(
        most_variable[display_columns].to_csv(
            sep="\t",
            index=False,
            float_format="%.5g",
        ).rstrip()
    )

    print(f"\nFull site table: {output_file}")


if __name__ == "__main__":
    main()
