#!/usr/bin/env python3

from pathlib import Path
import re
import subprocess

import numpy as np
import pandas as pd

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
GENOTYPE_DIR = PROJECT / "data/derived/region_genotypes"
VCF_DIR = PROJECT / "data/raw/region_panel"
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

VCF_NAMES = {
    panel: f"{panel}_filtered.vcf.gz"
    for panel in PANELS
}

FOCAL = {
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
}

QC_SCHEMES = {
    "called_GT": {
        "minimum_DP": None,
        "minimum_GQ": None,
    },
    "DP2_GQ20": {
        "minimum_DP": 2,
        "minimum_GQ": 20,
    },
    "DP5_GQ20": {
        "minimum_DP": 5,
        "minimum_GQ": 20,
    },
    "DP10_GQ30": {
        "minimum_DP": 10,
        "minimum_GQ": 30,
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


def genotype_called(value):
    value = str(value)

    if value in {"", ".", "./.", ".|.", "nan"}:
        return False

    alleles = re.split(r"[/|]", value)

    return (
        len(alleles) == 2
        and all(allele not in {"", "."} for allele in alleles)
    )


def qc_mask(calls, scheme):
    mask = calls["GT_called"].copy()

    minimum_dp = scheme["minimum_DP"]
    minimum_gq = scheme["minimum_GQ"]

    if minimum_dp is not None:
        mask &= (
            calls["DP_numeric"].notna()
            & calls["DP_numeric"].ge(minimum_dp)
        )

    if minimum_gq is not None:
        mask &= (
            calls["GQ_numeric"].notna()
            & calls["GQ_numeric"].ge(minimum_gq)
        )

    return mask


def quantile(values, probability):
    values = pd.to_numeric(values, errors="coerce").dropna()

    if values.empty:
        return np.nan

    return float(values.quantile(probability))


def extract_format_headers(panel, vcf):
    completed = subprocess.run(
        ["bcftools", "view", "-h", str(vcf)],
        check=True,
        capture_output=True,
        text=True,
    )

    rows = []

    for line in completed.stdout.splitlines():
        if not line.startswith("##FORMAT=<"):
            continue

        identifier_match = re.search(r"ID=([^,>]+)", line)
        number_match = re.search(r"Number=([^,>]+)", line)
        type_match = re.search(r"Type=([^,>]+)", line)
        description_match = re.search(
            r'Description="([^"]*)"',
            line,
        )

        rows.append({
            "panel": panel,
            "ID": (
                identifier_match.group(1)
                if identifier_match else ""
            ),
            "Number": (
                number_match.group(1)
                if number_match else ""
            ),
            "Type": (
                type_match.group(1)
                if type_match else ""
            ),
            "Description": (
                description_match.group(1)
                if description_match else ""
            ),
        })

    return rows


def main():
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    accession = resolve_column(metadata, "Accession")
    data_type = resolve_column(metadata, "Data Type")
    taxon = resolve_column(metadata, "Taxon Name")
    optix_figure = resolve_column(metadata, "Optix Figure")
    ivory_figure = resolve_column(metadata, "Ivory Figure")
    band = resolve_column(metadata, "B")
    ray = resolve_column(metadata, "R")
    dennis = resolve_column(metadata, "D")
    yb_dorsal = resolve_column(metadata, "Yb-dorsal")
    yb_ventral = resolve_column(metadata, "Yb-ventral")

    metadata = metadata.copy()
    metadata["sample"] = metadata[accession].astype(str)
    metadata["data_type_qc"] = metadata[data_type].astype(str)
    metadata["taxon_qc"] = metadata[taxon].astype(str)
    metadata["optix_figure_qc"] = numeric(
        metadata[optix_figure]
    )
    metadata["ivory_figure_qc"] = numeric(
        metadata[ivory_figure]
    )
    metadata["B_qc"] = numeric(metadata[band])
    metadata["R_qc"] = numeric(metadata[ray])
    metadata["D_qc"] = numeric(metadata[dennis])
    metadata["Ybd_qc"] = numeric(metadata[yb_dorsal])
    metadata["Ybv_qc"] = numeric(metadata[yb_ventral])

    piano_key_mask = metadata["taxon_qc"].str.contains(
        "piano keys",
        case=False,
        na=False,
    )

    metadata["role_qc"] = "comparative"
    metadata.loc[piano_key_mask, "role_qc"] = "piano_keys"
    metadata.loc[
        metadata["sample"].isin(FOCAL),
        "role_qc",
    ] = "focal"

    metadata_by_sample = metadata.set_index("sample")

    format_rows = []
    callability_rows = []
    trait_site_rows = []
    phenotype_count_rows = []

    for panel in PANELS:
        vcf = VCF_DIR / VCF_NAMES[panel]
        genotype_file = (
            GENOTYPE_DIR / f"{panel}.genotypes.tsv.gz"
        )

        format_rows.extend(
            extract_format_headers(panel, vcf)
        )

        calls = pd.read_csv(
            genotype_file,
            sep="\t",
            dtype=str,
            compression="gzip",
            keep_default_na=False,
        )

        calls["data_type"] = calls["sample"].map(
            metadata_by_sample["data_type_qc"]
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

        for trait in ["B", "R", "D", "Ybd", "Ybv"]:
            calls[trait] = calls["sample"].map(
                metadata_by_sample[f"{trait}_qc"]
            )

        calls["DP_numeric"] = numeric(calls["DP"])
        calls["GQ_numeric"] = numeric(calls["GQ"])
        calls["GT_called"] = calls["GT"].map(genotype_called)

        calls["biallelic_SNP"] = (
            ~calls["alt"].str.contains(",", regex=False)
            & calls["ref"].str.len().eq(1)
            & calls["alt"].str.len().eq(1)
        )

        comparative = calls["role"].eq("comparative")

        if panel == "optix":
            primary = (
                comparative
                & calls["optix_figure"].eq(1)
            )
            traits = ["B", "R", "D"]
        elif panel == "ivory":
            primary = (
                comparative
                & calls["ivory_figure"].eq(1)
            )
            traits = ["Ybd", "Ybv"]
        else:
            primary = comparative
            traits = []

        cohort_masks = {
            "comparative_all_assays": comparative,
            "comparative_CS": (
                comparative & calls["data_type"].eq("CS")
            ),
            "comparative_WGS": (
                comparative & calls["data_type"].eq("WGS")
            ),
        }

        if panel == "optix":
            cohort_masks["primary_Figure1"] = primary
            cohort_masks["primary_Figure1_CS"] = (
                primary & calls["data_type"].eq("CS")
            )
            cohort_masks["primary_Figure1_WGS"] = (
                primary & calls["data_type"].eq("WGS")
            )

        if panel == "ivory":
            cohort_masks["primary_Figure1"] = primary
            cohort_masks["primary_Figure1_CS"] = (
                primary & calls["data_type"].eq("CS")
            )
            cohort_masks["primary_Figure1_WGS"] = (
                primary & calls["data_type"].eq("WGS")
            )

        site_columns = ["chrom", "pos", "ref", "alt"]

        for scheme_name, scheme in QC_SCHEMES.items():
            calls["passes_current_qc"] = qc_mask(
                calls,
                scheme,
            )

            for cohort_name, cohort_mask in cohort_masks.items():
                cohort = calls.loc[
                    cohort_mask & calls["biallelic_SNP"]
                ].copy()

                site_rates = (
                    cohort
                    .groupby(site_columns, sort=False)
                    ["passes_current_qc"]
                    .agg(["size", "sum", "mean"])
                    .reset_index()
                )

                callability_rows.append({
                    "panel": panel,
                    "scheme": scheme_name,
                    "cohort": cohort_name,
                    "sites": int(len(site_rates)),
                    "samples": int(
                        cohort["sample"].nunique()
                    ),
                    "call_rate_p10": quantile(
                        site_rates["mean"], 0.10
                    ),
                    "call_rate_p25": quantile(
                        site_rates["mean"], 0.25
                    ),
                    "call_rate_median": quantile(
                        site_rates["mean"], 0.50
                    ),
                    "call_rate_p75": quantile(
                        site_rates["mean"], 0.75
                    ),
                    "call_rate_p90": quantile(
                        site_rates["mean"], 0.90
                    ),
                    "sites_ge_0.50": int(
                        site_rates["mean"].ge(0.50).sum()
                    ),
                    "sites_ge_0.70": int(
                        site_rates["mean"].ge(0.70).sum()
                    ),
                    "sites_ge_0.80": int(
                        site_rates["mean"].ge(0.80).sum()
                    ),
                    "sites_ge_0.90": int(
                        site_rates["mean"].ge(0.90).sum()
                    ),
                })

            for trait in traits:
                phenotype_metadata = metadata.loc[
                    (
                        metadata["role_qc"].eq("comparative")
                        & (
                            metadata["optix_figure_qc"].eq(1)
                            if panel == "optix"
                            else metadata["ivory_figure_qc"].eq(1)
                        )
                        & metadata[f"{trait}_qc"].isin([0, 1])
                    ),
                    ["sample", f"{trait}_qc"],
                ].copy()

                phenotype_counts = (
                    phenotype_metadata[f"{trait}_qc"]
                    .value_counts()
                    .to_dict()
                )

                if scheme_name == "called_GT":
                    phenotype_count_rows.append({
                        "panel": panel,
                        "trait": trait,
                        "phenotype_0_n": int(
                            phenotype_counts.get(0, 0)
                        ),
                        "phenotype_1_n": int(
                            phenotype_counts.get(1, 0)
                        ),
                    })

                trait_calls = calls.loc[
                    primary
                    & calls["biallelic_SNP"]
                    & calls[trait].isin([0, 1])
                ].copy()

                grouped = (
                    trait_calls
                    .groupby(
                        site_columns + [trait],
                        sort=False,
                    )["passes_current_qc"]
                    .agg(["size", "sum", "mean"])
                    .reset_index()
                )

                wide_n = grouped.pivot_table(
                    index=site_columns,
                    columns=trait,
                    values="sum",
                    fill_value=0,
                )

                wide_rate = grouped.pivot_table(
                    index=site_columns,
                    columns=trait,
                    values="mean",
                    fill_value=0,
                )

                for state in [0.0, 1.0]:
                    if state not in wide_n.columns:
                        wide_n[state] = 0
                    if state not in wide_rate.columns:
                        wide_rate[state] = 0

                summary = pd.DataFrame({
                    "state0_called": wide_n[0.0],
                    "state1_called": wide_n[1.0],
                    "state0_call_rate": wide_rate[0.0],
                    "state1_call_rate": wide_rate[1.0],
                }).reset_index()

                summary["min_state_called"] = summary[
                    ["state0_called", "state1_called"]
                ].min(axis=1)

                summary["min_state_call_rate"] = summary[
                    ["state0_call_rate", "state1_call_rate"]
                ].min(axis=1)

                trait_site_rows.append(pd.DataFrame({
                    "panel": panel,
                    "trait": trait,
                    "scheme": scheme_name,
                    "sites": int(len(summary)),
                    "min_state_called_median": quantile(
                        summary["min_state_called"], 0.50
                    ),
                    "sites_both_states_ge_10": int(
                        summary["min_state_called"].ge(10).sum()
                    ),
                    "sites_both_states_ge_20": int(
                        summary["min_state_called"].ge(20).sum()
                    ),
                    "sites_both_states_ge_30": int(
                        summary["min_state_called"].ge(30).sum()
                    ),
                    "sites_both_state_call_rates_ge_0.50": int(
                        summary["min_state_call_rate"]
                        .ge(0.50)
                        .sum()
                    ),
                    "sites_both_state_call_rates_ge_0.70": int(
                        summary["min_state_call_rate"]
                        .ge(0.70)
                        .sum()
                    ),
                }, index=[0]))

    format_table = pd.DataFrame(format_rows)
    callability_table = pd.DataFrame(callability_rows)
    trait_site_table = pd.concat(
        trait_site_rows,
        ignore_index=True,
    )
    phenotype_counts = pd.DataFrame(phenotype_count_rows)

    format_file = CHECKPOINT_DIR / "region_FORMAT_definitions.tsv"
    callability_file = (
        CHECKPOINT_DIR / "site_callability_qc_sensitivity.tsv"
    )
    trait_file = (
        CHECKPOINT_DIR / "phenotype_state_site_usability.tsv"
    )
    phenotype_file = (
        CHECKPOINT_DIR / "primary_phenotype_state_counts.tsv"
    )

    format_table.to_csv(format_file, sep="\t", index=False)
    callability_table.to_csv(
        callability_file,
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    trait_site_table.to_csv(
        trait_file,
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    phenotype_counts.to_csv(
        phenotype_file,
        sep="\t",
        index=False,
    )

    selected_format = format_table.loc[
        format_table["ID"].isin(
            ["GT", "DP", "GQ", "GP", "PL", "AD"]
        )
    ].drop_duplicates(
        subset=["ID", "Number", "Type", "Description"]
    )

    primary_callability = callability_table.loc[
        callability_table["cohort"].isin(
            [
                "primary_Figure1",
                "primary_Figure1_CS",
                "primary_Figure1_WGS",
            ]
        )
    ]

    print("RELEVANT FORMAT DEFINITIONS")
    print(
        selected_format.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print("\nPRIMARY PHENOTYPE COUNTS")
    print(
        phenotype_counts.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print("\nPRIMARY-COHORT SITE CALLABILITY")
    print(
        primary_callability.to_csv(
            sep="\t",
            index=False,
            float_format="%.4g",
        ).rstrip()
    )

    print("\nPHENOTYPE-STATE SITE USABILITY")
    print(
        trait_site_table.to_csv(
            sep="\t",
            index=False,
            float_format="%.4g",
        ).rstrip()
    )

    print("\nDetailed outputs:")
    print(format_file)
    print(callability_file)
    print(trait_file)
    print(phenotype_file)


if __name__ == "__main__":
    main()
