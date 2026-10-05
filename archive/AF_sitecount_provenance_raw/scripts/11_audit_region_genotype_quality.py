#!/usr/bin/env python3

from pathlib import Path
import re
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=RuntimeWarning)

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
GENOTYPE_DIR = PROJECT / "data/derived/region_genotypes"
METADATA_FILE = PROJECT / "data/derived/region_panel_metadata.tsv"
CHECKPOINT_DIR = PROJECT / "analysis/checkpoints"

PANELS = [
    "optix",
    "ivory",
    "WntA",
    "neutral1",
    "neutral2",
    "neutral3",
    "neutral4",
]

FOCAL_INDIVIDUALS = {
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
}

FOCAL_POOLED = {
    "TN_Left",
    "TN_Right",
}


def normalize_name(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def resolve_column(dataframe, requested):
    requested_norm = normalize_name(requested)
    matches = [
        column for column in dataframe.columns
        if normalize_name(column) == requested_norm
    ]

    if len(matches) != 1:
        raise ValueError(
            f"Could not uniquely resolve metadata column {requested!r}; "
            f"matches={matches}"
        )

    return matches[0]


def numeric(series):
    return pd.to_numeric(
        series.replace({".": np.nan, "": np.nan}),
        errors="coerce",
    )


def quantile(series, probability):
    values = pd.to_numeric(series, errors="coerce").dropna()

    if values.empty:
        return np.nan

    return float(values.quantile(probability))


def genotype_alleles(value):
    value = str(value)

    if value in {".", "./.", ".|.", "nan", "None"}:
        return None

    alleles = re.split(r"[/|]", value)

    if len(alleles) != 2 or any(allele == "." for allele in alleles):
        return None

    try:
        return tuple(int(allele) for allele in alleles)
    except ValueError:
        return None


def maximum_probability(value):
    if pd.isna(value):
        return np.nan

    text = str(value)

    if text in {"", "."}:
        return np.nan

    try:
        values = [
            float(item)
            for item in text.split(",")
            if item not in {"", "."}
        ]
    except ValueError:
        return np.nan

    return max(values) if values else np.nan


def heterozygote_alt_balance(ad_value, alleles, is_biallelic):
    if alleles is None or not is_biallelic:
        return np.nan

    if set(alleles) != {0, 1}:
        return np.nan

    if pd.isna(ad_value):
        return np.nan

    fields = str(ad_value).split(",")

    if len(fields) < 2 or fields[0] == "." or fields[1] == ".":
        return np.nan

    try:
        ref_depth = float(fields[0])
        alt_depth = float(fields[1])
    except ValueError:
        return np.nan

    total = ref_depth + alt_depth

    if total <= 0:
        return np.nan

    return alt_depth / total


def summarize_genotypes(subset, panel, stratum):
    called = subset["called"]
    called_subset = subset.loc[called]

    return {
        "panel": panel,
        "stratum": stratum,
        "samples": int(subset["sample"].nunique()),
        "genotypes": int(len(subset)),
        "called_genotypes": int(called.sum()),
        "call_rate": float(called.mean()) if len(subset) else np.nan,
        "DP_available_fraction": float(
            subset["DP_numeric"].notna().mean()
        ) if len(subset) else np.nan,
        "DP_p10": quantile(called_subset["DP_numeric"], 0.10),
        "DP_median": quantile(called_subset["DP_numeric"], 0.50),
        "DP_p90": quantile(called_subset["DP_numeric"], 0.90),
        "GQ_available_fraction": float(
            subset["GQ_numeric"].notna().mean()
        ) if len(subset) else np.nan,
        "GQ_p10": quantile(called_subset["GQ_numeric"], 0.10),
        "GQ_median": quantile(called_subset["GQ_numeric"], 0.50),
        "GQ_p90": quantile(called_subset["GQ_numeric"], 0.90),
        "GP_available_fraction": float(
            subset["max_GP"].notna().mean()
        ) if len(subset) else np.nan,
        "max_GP_p10": quantile(called_subset["max_GP"], 0.10),
        "max_GP_median": quantile(called_subset["max_GP"], 0.50),
        "heterozygous_calls": int(subset["is_heterozygous"].sum()),
        "het_AB_available": int(subset["het_alt_balance"].notna().sum()),
        "het_AB_p10": quantile(subset["het_alt_balance"], 0.10),
        "het_AB_median": quantile(subset["het_alt_balance"], 0.50),
        "het_AB_p90": quantile(subset["het_alt_balance"], 0.90),
    }


def summarize_sample(subset, panel, sample):
    called_subset = subset.loc[subset["called"]]

    return {
        "panel": panel,
        "sample": sample,
        "data_type": subset["data_type"].iloc[0],
        "role": subset["role"].iloc[0],
        "taxon": subset["taxon"].iloc[0],
        "sites": int(len(subset)),
        "called_sites": int(subset["called"].sum()),
        "call_rate": float(subset["called"].mean()),
        "DP_median": quantile(called_subset["DP_numeric"], 0.50),
        "GQ_median": quantile(called_subset["GQ_numeric"], 0.50),
        "max_GP_median": quantile(called_subset["max_GP"], 0.50),
        "heterozygous_calls": int(subset["is_heterozygous"].sum()),
        "het_AB_median": quantile(subset["het_alt_balance"], 0.50),
    }


def main():
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    accession_column = resolve_column(metadata, "Accession")
    data_type_column = resolve_column(metadata, "Data Type")
    taxon_column = resolve_column(metadata, "Taxon Name")
    optix_figure_column = resolve_column(metadata, "Optix Figure")
    ivory_figure_column = resolve_column(metadata, "Ivory Figure")

    metadata = metadata.copy()
    metadata["sample"] = metadata[accession_column].astype(str)
    metadata["data_type_qc"] = metadata[data_type_column].astype(str)
    metadata["taxon_qc"] = metadata[taxon_column].astype(str)
    metadata["optix_figure_qc"] = numeric(metadata[optix_figure_column])
    metadata["ivory_figure_qc"] = numeric(metadata[ivory_figure_column])

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
        metadata["sample"].isin(FOCAL_INDIVIDUALS),
        "role_qc",
    ] = "focal_individual"

    metadata.loc[
        metadata["sample"].isin(FOCAL_POOLED),
        "role_qc",
    ] = "focal_pooled"

    metadata_by_sample = metadata.set_index("sample")

    variant_rows = []
    genotype_summary_rows = []
    sample_summary_rows = []
    site_summary_rows = []

    for panel in PANELS:
        genotype_file = GENOTYPE_DIR / f"{panel}.genotypes.tsv.gz"

        if not genotype_file.is_file():
            raise FileNotFoundError(genotype_file)

        calls = pd.read_csv(
            genotype_file,
            sep="\t",
            dtype=str,
            compression="gzip",
            keep_default_na=False,
        )

        missing_metadata = sorted(
            set(calls["sample"]) - set(metadata_by_sample.index)
        )

        if missing_metadata:
            raise ValueError(
                f"{panel}: {len(missing_metadata)} samples lack metadata"
            )

        calls["data_type"] = calls["sample"].map(
            metadata_by_sample["data_type_qc"]
        )
        calls["taxon"] = calls["sample"].map(
            metadata_by_sample["taxon_qc"]
        )
        calls["role"] = calls["sample"].map(
            metadata_by_sample["role_qc"]
        )
        calls["optix_figure"] = calls["sample"].map(
            metadata_by_sample["optix_figure_qc"]
        )
        calls["ivory_figure"] = calls["sample"].map(
            metadata_by_sample["ivory_figure_qc"]
        )

        calls["alleles"] = calls["GT"].map(genotype_alleles)
        calls["called"] = calls["alleles"].notna()
        calls["is_heterozygous"] = calls["alleles"].map(
            lambda value: (
                value is not None and value[0] != value[1]
            )
        )

        calls["DP_numeric"] = numeric(calls["DP"])
        calls["GQ_numeric"] = numeric(calls["GQ"])
        calls["max_GP"] = calls["GP"].map(maximum_probability)

        calls["is_biallelic"] = ~calls["alt"].str.contains(
            ",",
            regex=False,
        )

        calls["het_alt_balance"] = [
            heterozygote_alt_balance(ad, alleles, biallelic)
            for ad, alleles, biallelic in zip(
                calls["AD"],
                calls["alleles"],
                calls["is_biallelic"],
            )
        ]

        site_columns = [
            "chrom",
            "pos",
            "ref",
            "alt",
            "qual",
            "filter",
        ]

        sites = calls[site_columns].drop_duplicates().copy()
        sites["is_biallelic"] = ~sites["alt"].str.contains(
            ",",
            regex=False,
        )
        sites["is_snp"] = (
            sites["is_biallelic"]
            & sites["ref"].str.len().eq(1)
            & sites["alt"].str.len().eq(1)
        )
        sites["is_indel_or_complex"] = ~sites["is_snp"]

        variant_rows.append({
            "panel": panel,
            "variants": int(len(sites)),
            "PASS_or_dot": int(
                sites["filter"].isin(["PASS", "."]).sum()
            ),
            "filtered": int(
                (~sites["filter"].isin(["PASS", "."])).sum()
            ),
            "biallelic": int(sites["is_biallelic"].sum()),
            "multiallelic": int((~sites["is_biallelic"]).sum()),
            "biallelic_SNPs": int(sites["is_snp"].sum()),
            "indel_or_complex": int(
                sites["is_indel_or_complex"].sum()
            ),
        })

        comparative_no_pk = (
            calls["role"].eq("comparative")
        )

        strata = {
            "comparative_no_piano_keys": comparative_no_pk,
            "comparative_CS": (
                comparative_no_pk & calls["data_type"].eq("CS")
            ),
            "comparative_WGS": (
                comparative_no_pk & calls["data_type"].eq("WGS")
            ),
            "piano_keys_excluded": calls["role"].eq("piano_keys"),
            "focal_individuals_QC_only": calls["role"].eq(
                "focal_individual"
            ),
            "focal_pooled_main": calls["role"].eq("focal_pooled"),
        }

        if panel == "optix":
            strata["primary_Optix_Figure1"] = (
                comparative_no_pk
                & calls["optix_figure"].eq(1)
            )

        if panel == "ivory":
            strata["primary_Ivory_Figure1"] = (
                comparative_no_pk
                & calls["ivory_figure"].eq(1)
            )

        for stratum_name, mask in strata.items():
            subset = calls.loc[mask].copy()

            if not subset.empty:
                genotype_summary_rows.append(
                    summarize_genotypes(
                        subset,
                        panel,
                        stratum_name,
                    )
                )

        for sample, subset in calls.groupby("sample", sort=False):
            sample_summary_rows.append(
                summarize_sample(subset, panel, sample)
            )

        site_key = ["chrom", "pos", "ref", "alt"]

        for site_values, subset in calls.groupby(
            site_key,
            sort=False,
            dropna=False,
        ):
            row = {
                "panel": panel,
                "chrom": site_values[0],
                "pos": site_values[1],
                "ref": site_values[2],
                "alt": site_values[3],
            }

            for label, mask in {
                "comparative": subset["role"].eq("comparative"),
                "CS": (
                    subset["role"].eq("comparative")
                    & subset["data_type"].eq("CS")
                ),
                "WGS": (
                    subset["role"].eq("comparative")
                    & subset["data_type"].eq("WGS")
                ),
            }.items():
                site_subset = subset.loc[mask]
                row[f"{label}_n"] = int(len(site_subset))
                row[f"{label}_called"] = int(
                    site_subset["called"].sum()
                )
                row[f"{label}_call_rate"] = (
                    float(site_subset["called"].mean())
                    if len(site_subset)
                    else np.nan
                )

            site_summary_rows.append(row)

    variants = pd.DataFrame(variant_rows)
    genotype_summary = pd.DataFrame(genotype_summary_rows)
    sample_summary = pd.DataFrame(sample_summary_rows)
    site_summary = pd.DataFrame(site_summary_rows)

    variants_file = (
        CHECKPOINT_DIR / "region_variant_structure.tsv"
    )
    genotype_file = (
        CHECKPOINT_DIR / "region_genotype_qc_by_stratum.tsv"
    )
    sample_file = (
        CHECKPOINT_DIR / "region_genotype_qc_by_sample.tsv"
    )
    site_file = (
        CHECKPOINT_DIR / "region_genotype_qc_by_site.tsv"
    )

    variants.to_csv(variants_file, sep="\t", index=False)
    genotype_summary.to_csv(
        genotype_file,
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    sample_summary.to_csv(
        sample_file,
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    site_summary.to_csv(
        site_file,
        sep="\t",
        index=False,
        float_format="%.6g",
    )

    display_columns = [
        "panel",
        "stratum",
        "samples",
        "call_rate",
        "DP_p10",
        "DP_median",
        "DP_p90",
        "GQ_p10",
        "GQ_median",
        "GP_available_fraction",
        "max_GP_p10",
        "max_GP_median",
        "heterozygous_calls",
        "het_AB_median",
    ]

    core_strata = {
        "comparative_no_piano_keys",
        "comparative_CS",
        "comparative_WGS",
        "primary_Optix_Figure1",
        "primary_Ivory_Figure1",
        "piano_keys_excluded",
        "focal_individuals_QC_only",
        "focal_pooled_main",
    }

    display = genotype_summary.loc[
        genotype_summary["stratum"].isin(core_strata),
        display_columns,
    ].copy()

    low_call_summary = (
        sample_summary
        .assign(
            below_90=sample_summary["call_rate"] < 0.90,
            below_95=sample_summary["call_rate"] < 0.95,
            below_99=sample_summary["call_rate"] < 0.99,
        )
        .groupby(["panel", "data_type", "role"], dropna=False)
        .agg(
            samples=("sample", "size"),
            call_rate_median=("call_rate", "median"),
            samples_below_90=("below_90", "sum"),
            samples_below_95=("below_95", "sum"),
            samples_below_99=("below_99", "sum"),
        )
        .reset_index()
    )

    print("VARIANT STRUCTURE")
    print(variants.to_csv(sep="\t", index=False).rstrip())

    print("\nGENOTYPE QC BY CORE STRATUM")
    print(
        display.to_csv(
            sep="\t",
            index=False,
            float_format="%.4g",
        ).rstrip()
    )

    print("\nLOW-CALL-RATE SAMPLE COUNTS")
    print(
        low_call_summary.to_csv(
            sep="\t",
            index=False,
            float_format="%.4g",
        ).rstrip()
    )

    print("\nDetailed outputs:")
    print(variants_file)
    print(genotype_file)
    print(sample_file)
    print(site_file)


if __name__ == "__main__":
    main()
