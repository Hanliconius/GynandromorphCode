#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
GENOTYPE_FILE = (
    PROJECT / "data/derived/region_genotypes/optix.genotypes.tsv.gz"
)
METADATA_FILE = (
    PROJECT / "data/derived/region_panel_metadata.tsv"
)
CHECKPOINT_DIR = PROJECT / "analysis/checkpoints"

FOCAL = {
    "TN2L_S1", "TN2R_S2", "TN3L_S3",
    "TN3R_S4", "TN_Left", "TN_Right",
}

BAND_START = 1_112_000
BAND_END = 1_116_000


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


def summarize_site(subset):
    passing = subset.loc[subset["passes_qc"]]
    counts = passing["dosage"].value_counts()

    ref_hom = int(counts.get(0.0, 0))
    het = int(counts.get(1.0, 0))
    alt_hom = int(counts.get(2.0, 0))
    called = ref_hom + het + alt_hom

    if called:
        alt_frequency = (het + 2 * alt_hom) / (2 * called)
        maf = min(alt_frequency, 1 - alt_frequency)
    else:
        maf = np.nan

    return {
        "samples": int(subset["sample"].nunique()),
        "called_n": called,
        "call_rate": (
            called / len(subset) if len(subset) else np.nan
        ),
        "ref_hom_n": ref_hom,
        "het_n": het,
        "alt_hom_n": alt_hom,
        "MAF": maf,
    }


def main():
    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    accession = resolve_column(metadata, "Accession")
    taxon = resolve_column(metadata, "Taxon Name")
    figure = resolve_column(metadata, "Optix Figure")

    metadata["sample"] = metadata[accession].astype(str)
    metadata["taxon_qc"] = metadata[taxon].astype(str)
    metadata["figure_qc"] = numeric(metadata[figure])
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

    calls = pd.read_csv(
        GENOTYPE_FILE,
        sep="\t",
        dtype=str,
        compression="gzip",
        keep_default_na=False,
    )

    calls["role"] = calls["sample"].map(
        metadata_by_sample["role_qc"]
    )
    calls["figure"] = calls["sample"].map(
        metadata_by_sample["figure_qc"]
    )
    calls["dosage"] = calls["GT"].map(dosage)
    calls["DP_numeric"] = numeric(calls["DP"])
    calls["GQ_numeric"] = numeric(calls["GQ"])
    calls["pos_numeric"] = numeric(calls["pos"])

    calls["biallelic_SNP"] = (
        ~calls["alt"].str.contains(",", regex=False)
        & calls["ref"].str.len().eq(1)
        & calls["alt"].str.len().eq(1)
    )

    calls["passes_qc"] = (
        calls["dosage"].notna()
        & calls["DP_numeric"].ge(2)
        & calls["GQ_numeric"].ge(20)
    )

    calls = calls.loc[calls["biallelic_SNP"]].copy()

    strata = {
        "all_comparative": calls["role"].eq("comparative"),
        "Optix_Figure_0": (
            calls["role"].eq("comparative")
            & calls["figure"].eq(0)
        ),
        "Optix_Figure_1": (
            calls["role"].eq("comparative")
            & calls["figure"].eq(1)
        ),
        "piano_keys_excluded": calls["role"].eq("piano_keys"),
        "focal_pooled": calls["sample"].isin(
            ["TN_Left", "TN_Right"]
        ),
    }

    site_rows = []

    for site_key, site in calls.groupby(
        ["chrom", "pos", "ref", "alt"],
        sort=False,
    ):
        for stratum, mask in strata.items():
            subset = site.loc[mask.loc[site.index]]
            summary = summarize_site(subset)

            site_rows.append({
                "chrom": site_key[0],
                "pos": int(site_key[1]),
                "ref": site_key[2],
                "alt": site_key[3],
                "stratum": stratum,
                **summary,
            })

    site_table = pd.DataFrame(site_rows)

    summary_rows = []

    for stratum, subset in site_table.groupby(
        "stratum",
        sort=False,
    ):
        summary_rows.append({
            "stratum": stratum,
            "samples": int(subset["samples"].max()),
            "sites": len(subset),
            "polymorphic_sites": int(subset["MAF"].gt(0).sum()),
            "sites_MAF_ge_0.01": int(
                subset["MAF"].ge(0.01).sum()
            ),
            "sites_MAF_ge_0.05": int(
                subset["MAF"].ge(0.05).sum()
            ),
            "sites_MAF_ge_0.10": int(
                subset["MAF"].ge(0.10).sum()
            ),
            "maximum_MAF": subset["MAF"].max(),
        })

    summary = pd.DataFrame(summary_rows)

    interval = site_table.loc[
        site_table["pos"].between(BAND_START, BAND_END)
    ].copy()

    interval_wide = interval.pivot_table(
        index=["chrom", "pos", "ref", "alt"],
        columns="stratum",
        values=[
            "called_n",
            "ref_hom_n",
            "het_n",
            "alt_hom_n",
            "MAF",
        ],
        aggfunc="first",
    )

    interval_wide.columns = [
        f"{metric}__{stratum}"
        for metric, stratum in interval_wide.columns
    ]
    interval_wide = interval_wide.reset_index()

    site_file = (
        CHECKPOINT_DIR
        / "optix_genotype_diversity_by_figure_flag.tsv"
    )
    interval_file = (
        CHECKPOINT_DIR
        / "optix_band_1.112_1.116Mb_genotype_counts.tsv"
    )

    site_table.to_csv(
        site_file,
        sep="\t",
        index=False,
        float_format="%.8g",
    )
    interval_wide.to_csv(
        interval_file,
        sep="\t",
        index=False,
        float_format="%.8g",
    )

    print("OPTIX DIVERSITY BY FIGURE FLAG")
    print(
        summary.to_csv(
            sep="\t",
            index=False,
            float_format="%.5g",
        ).rstrip()
    )

    print(
        f"\nBAND INTERVAL {BAND_START}-{BAND_END}"
    )

    if interval_wide.empty:
        print("No biallelic SNP records occur in this interval.")
    else:
        print(
            interval_wide.to_csv(
                sep="\t",
                index=False,
                float_format="%.5g",
            ).rstrip()
        )

    print("\nDetailed outputs:")
    print(site_file)
    print(interval_file)


if __name__ == "__main__":
    main()
