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
    "optix": "Optix Figure",
    "ivory": "Ivory Figure",
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


def gt_dosage(value):
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


def parse_gp(value):
    value = str(value)

    if value in {"", ".", "nan"}:
        return (np.nan, np.nan, np.nan, np.nan)

    fields = value.split(",")

    if len(fields) != 3:
        return (np.nan, np.nan, np.nan, np.nan)

    try:
        probabilities = np.asarray(
            [float(field) for field in fields],
            dtype=float,
        )
    except ValueError:
        return (np.nan, np.nan, np.nan, np.nan)

    total = probabilities.sum()

    if (
        not np.all(np.isfinite(probabilities))
        or total <= 0
        or np.any(probabilities < 0)
    ):
        return (np.nan, np.nan, np.nan, np.nan)

    probabilities = probabilities / total

    expected_alt_dosage = (
        probabilities[1] + 2 * probabilities[2]
    )
    maximum_probability = float(probabilities.max())
    most_likely_dosage = float(np.argmax(probabilities))

    return (
        expected_alt_dosage,
        maximum_probability,
        most_likely_dosage,
        total,
    )


def quantile(series, probability):
    values = pd.to_numeric(series, errors="coerce").dropna()

    if values.empty:
        return np.nan

    return float(values.quantile(probability))


def main():
    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    accession = resolve_column(metadata, "Accession")
    taxon = resolve_column(metadata, "Taxon Name")

    metadata["sample"] = metadata[accession].astype(str)
    metadata["taxon_qc"] = metadata[taxon].astype(str)
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

    availability_rows = []
    site_rows = []

    for panel, figure_name in PANELS.items():
        figure_column = resolve_column(metadata, figure_name)
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

        calls["biallelic_SNP"] = (
            ~calls["alt"].str.contains(",", regex=False)
            & calls["ref"].str.len().eq(1)
            & calls["alt"].str.len().eq(1)
        )

        calls = calls.loc[
            calls["role"].eq("comparative")
            & calls["figure"].eq(1)
            & calls["biallelic_SNP"]
        ].copy()

        calls["GT_dosage"] = calls["GT"].map(gt_dosage)

        gp_values = calls["GP"].map(parse_gp)
        calls["GP_expected_dosage"] = gp_values.map(
            lambda value: value[0]
        )
        calls["GP_max_probability"] = gp_values.map(
            lambda value: value[1]
        )
        calls["GP_most_likely_dosage"] = gp_values.map(
            lambda value: value[2]
        )
        calls["GP_original_sum"] = gp_values.map(
            lambda value: value[3]
        )

        calls["GT_status"] = np.where(
            calls["GT_dosage"].notna(),
            "hard_GT_present",
            "hard_GT_missing",
        )

        for status, subset in calls.groupby(
            "GT_status",
            sort=False,
        ):
            valid_gp = subset["GP_expected_dosage"].notna()
            comparable = (
                subset["GT_dosage"].notna()
                & subset["GP_most_likely_dosage"].notna()
            )

            availability_rows.append({
                "panel": panel,
                "GT_status": status,
                "rows": len(subset),
                "valid_GP_fraction": float(valid_gp.mean()),
                "GP_max_probability_p10": quantile(
                    subset["GP_max_probability"], 0.10
                ),
                "GP_max_probability_median": quantile(
                    subset["GP_max_probability"], 0.50
                ),
                "GP_max_probability_p90": quantile(
                    subset["GP_max_probability"], 0.90
                ),
                "GP_max_probability_ge_0.90_fraction": float(
                    subset["GP_max_probability"].ge(0.90).mean()
                ),
                "GP_max_probability_ge_0.99_fraction": float(
                    subset["GP_max_probability"].ge(0.99).mean()
                ),
                "GT_GP_argmax_concordance": (
                    float(
                        (
                            subset.loc[
                                comparable,
                                "GT_dosage",
                            ]
                            == subset.loc[
                                comparable,
                                "GP_most_likely_dosage",
                            ]
                        ).mean()
                    )
                    if comparable.any() else np.nan
                ),
            })

        for site_key, site in calls.groupby(
            ["chrom", "pos", "ref", "alt"],
            sort=False,
        ):
            valid_gp = site.loc[
                site["GP_expected_dosage"].notna()
            ]
            high_confidence_gp = site.loc[
                site["GP_max_probability"].ge(0.90)
            ]

            gp_alt_frequency = (
                valid_gp["GP_expected_dosage"].sum()
                / (2 * len(valid_gp))
                if len(valid_gp) else np.nan
            )

            high_confidence_alt_frequency = (
                high_confidence_gp[
                    "GP_expected_dosage"
                ].sum()
                / (2 * len(high_confidence_gp))
                if len(high_confidence_gp) else np.nan
            )

            gt = site.loc[site["GT_dosage"].notna()]
            gt_alt_frequency = (
                gt["GT_dosage"].sum() / (2 * len(gt))
                if len(gt) else np.nan
            )

            def maf(frequency):
                if not np.isfinite(frequency):
                    return np.nan
                return min(frequency, 1 - frequency)

            site_rows.append({
                "panel": panel,
                "chrom": site_key[0],
                "pos": int(site_key[1]),
                "ref": site_key[2],
                "alt": site_key[3],
                "cohort_n": len(site),
                "hard_GT_n": len(gt),
                "valid_GP_n": len(valid_gp),
                "high_confidence_GP_n": len(
                    high_confidence_gp
                ),
                "GT_MAF": maf(gt_alt_frequency),
                "GP_MAF": maf(gp_alt_frequency),
                "high_confidence_GP_MAF": maf(
                    high_confidence_alt_frequency
                ),
                "GP_max_probability_median": quantile(
                    site["GP_max_probability"], 0.50
                ),
            })

    availability = pd.DataFrame(availability_rows)
    sites = pd.DataFrame(site_rows)

    availability_file = (
        CHECKPOINT_DIR / "GP_availability_by_GT_status.tsv"
    )
    site_file = (
        CHECKPOINT_DIR / "GP_derived_site_allele_frequencies.tsv"
    )

    availability.to_csv(
        availability_file,
        sep="\t",
        index=False,
        float_format="%.8g",
    )
    sites.to_csv(
        site_file,
        sep="\t",
        index=False,
        float_format="%.8g",
    )

    site_summary_rows = []

    for panel, subset in sites.groupby("panel", sort=False):
        site_summary_rows.append({
            "panel": panel,
            "sites": len(subset),
            "sites_GT_MAF_ge_0.01": int(
                subset["GT_MAF"].ge(0.01).sum()
            ),
            "sites_GT_MAF_ge_0.05": int(
                subset["GT_MAF"].ge(0.05).sum()
            ),
            "sites_GP_MAF_ge_0.01": int(
                subset["GP_MAF"].ge(0.01).sum()
            ),
            "sites_GP_MAF_ge_0.05": int(
                subset["GP_MAF"].ge(0.05).sum()
            ),
            "sites_high_conf_GP_MAF_ge_0.01": int(
                subset["high_confidence_GP_MAF"]
                .ge(0.01)
                .sum()
            ),
            "sites_high_conf_GP_MAF_ge_0.05": int(
                subset["high_confidence_GP_MAF"]
                .ge(0.05)
                .sum()
            ),
            "valid_GP_n_median": quantile(
                subset["valid_GP_n"], 0.50
            ),
            "high_confidence_GP_n_median": quantile(
                subset["high_confidence_GP_n"], 0.50
            ),
        })

    site_summary = pd.DataFrame(site_summary_rows)

    most_variable = (
        sites.sort_values(
            ["panel", "GP_MAF", "valid_GP_n"],
            ascending=[True, False, False],
        )
        .groupby("panel", sort=False)
        .head(12)
    )

    display_columns = [
        "panel", "chrom", "pos", "ref", "alt",
        "cohort_n", "hard_GT_n", "valid_GP_n",
        "high_confidence_GP_n",
        "GT_MAF", "GP_MAF",
        "high_confidence_GP_MAF",
        "GP_max_probability_median",
    ]

    print("GP AVAILABILITY BY HARD-GT STATUS")
    print(
        availability.to_csv(
            sep="\t",
            index=False,
            float_format="%.5g",
        ).rstrip()
    )

    print("\nSITE-LEVEL GP INFORMATION SUMMARY")
    print(
        site_summary.to_csv(
            sep="\t",
            index=False,
            float_format="%.5g",
        ).rstrip()
    )

    print("\nMOST VARIABLE SITES BY GP EXPECTED DOSAGE")
    print(
        most_variable[display_columns].to_csv(
            sep="\t",
            index=False,
            float_format="%.5g",
        ).rstrip()
    )

    print("\nDetailed outputs:")
    print(availability_file)
    print(site_file)


if __name__ == "__main__":
    main()
