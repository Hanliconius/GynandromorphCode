#!/usr/bin/env python3

from pathlib import Path
import math
import re

import numpy as np
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

CHECKPOINT_DIR = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
)

SELECTED_SITES = (
    CHECKPOINT_DIR
    / "broad_optix_strict_comparative_SNPs.tsv"
)

GENOTYPES = (
    CHECKPOINT_DIR
    / "core_explanatory_optix_genotypes.long.tsv"
)

METADATA = (
    PROJECT
    / "data"
    / "derived"
    / "region_panel_metadata.tsv"
)

OUT_ALL = (
    CHECKPOINT_DIR
    / "B_vs_DR_canonical_discordant_evaluation.tsv"
)

OUT_EXACT = (
    CHECKPOINT_DIR
    / "B_vs_DR_exact_focal_zygosity.canonical.tsv"
)

OUT_EXACT_QC = (
    CHECKPOINT_DIR
    / "B_vs_DR_exact_focal_zygosity.canonical.DP2_GQ20.tsv"
)

OUT_EXACT_SPECIFIC = (
    CHECKPOINT_DIR
    / "B_vs_DR_exact_focal_zygosity."
      "canonical.discordantBA80.tsv"
)

OUT_SUMMARY = (
    CHECKPOINT_DIR
    / "B_vs_DR_exact_focal_zygosity."
      "canonical.summary.tsv"
)


def clean_header(value):
    return " ".join(
        str(value).split()
    )


def find_column(
    dataframe,
    candidates,
    required=True,
):
    lookup = {
        clean_header(column).lower(): column
        for column in dataframe.columns
    }

    for candidate in candidates:
        key = clean_header(candidate).lower()

        if key in lookup:
            return lookup[key]

    if required:
        raise RuntimeError(
            "Required column was not identified.\n"
            f"Tried: {', '.join(candidates)}\n"
            "Available columns:\n"
            + "\n".join(
                str(column)
                for column in dataframe.columns
            )
        )

    return None


def boolean(value):
    if isinstance(value, bool):
        return value

    if pd.isna(value):
        return False

    return str(value).strip().lower() in {
        "true",
        "t",
        "1",
        "yes",
    }


def numeric(value):
    result = pd.to_numeric(
        value,
        errors="coerce",
    )

    if pd.isna(result):
        return np.nan

    return float(result)


def binary(value):
    result = numeric(value)

    if not np.isfinite(result):
        return np.nan

    if abs(result) < 1e-9:
        return 0.0

    if abs(result - 1.0) < 1e-9:
        return 1.0

    return np.nan


def canonical_gt(value):
    if pd.isna(value):
        return None

    text = str(value).strip()

    if text in {
        "",
        ".",
        "./.",
        ".|.",
        "missing",
        "NA",
        "nan",
    }:
        return None

    alleles = re.split(r"[|/]", text)

    if (
        len(alleles) != 2
        or any(
            allele in {"", "."}
            for allele in alleles
        )
    ):
        return None

    try:
        numeric_alleles = sorted(
            int(allele)
            for allele in alleles
        )
    except ValueError:
        return None

    if any(
        allele not in {0, 1}
        for allele in numeric_alleles
    ):
        return None

    return "/".join(
        str(allele)
        for allele in numeric_alleles
    )


def normalize_biological_gt(value):
    if pd.isna(value):
        return None

    text = (
        str(value)
        .strip()
        .upper()
        .replace("/", "")
        .replace("|", "")
        .replace("-", "")
        .replace("_", "")
    )

    if text == "AA":
        return "A/A"

    if text in {"PA", "AP"}:
        return "P/A"

    if text == "PP":
        return "P/P"

    return None


def focal_class(module):
    if str(module).strip().lower() == "band":
        return "B"

    return "DR"


def comparative_trait(module):
    text = str(module).strip().lower()

    if text == "dennis":
        return "D"

    if text == "ray":
        return "R"

    if text == "band":
        return "B"

    raise RuntimeError(
        f"Unexpected module: {module}"
    )


def biological_gt_from_raw(gt, module):
    gt = canonical_gt(gt)

    if gt is None:
        return None

    if str(module).strip().lower() == "band":
        mapping = {
            "0/0": "P/P",
            "0/1": "P/A",
            "1/1": "A/A",
        }
    else:
        mapping = {
            "0/0": "A/A",
            "0/1": "P/A",
            "1/1": "P/P",
        }

    return mapping.get(gt)


def dosage(biological_gt):
    return {
        "A/A": 0.0,
        "P/A": 0.5,
        "P/P": 1.0,
    }.get(biological_gt, np.nan)


def expected_focal_states(module):
    if focal_class(module) == "B":
        return {
            "Left": {
                "GT": "A/A",
                "dosage": 0.0,
            },
            "Right": {
                "GT": "P/A",
                "dosage": 0.5,
            },
        }

    return {
        "Left": {
            "GT": "P/P",
            "dosage": 1.0,
        },
        "Right": {
            "GT": "P/A",
            "dosage": 0.5,
        },
    }


def safe_fraction(numerator, denominator):
    if denominator <= 0:
        return np.nan

    return numerator / denominator


def main():
    for path in [
        SELECTED_SITES,
        GENOTYPES,
        METADATA,
    ]:
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing required input: {path}"
            )

    sites = pd.read_csv(
        SELECTED_SITES,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    genotypes = pd.read_csv(
        GENOTYPES,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    metadata = pd.read_csv(
        METADATA,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    metadata.columns = [
        clean_header(column)
        for column in metadata.columns
    ]

    # ------------------------------------------------------------
    # Resolve selected-site columns
    # ------------------------------------------------------------

    site_module = find_column(
        sites,
        ["module", "source_module"],
    )

    site_chrom = find_column(
        sites,
        ["chrom", "CHROM", "chromosome"],
    )

    site_pos = find_column(
        sites,
        ["pos", "POS", "position"],
    )

    site_ref = find_column(
        sites,
        ["ref", "REF"],
    )

    site_alt = find_column(
        sites,
        ["alt", "ALT"],
    )

    site_ba = find_column(
        sites,
        [
            "minimum_BA",
            "minimum_comparative_BA",
            "minimum_BA_all_cohorts_schemes",
            "minimum_BA_across_cohorts_and_schemes",
        ],
        required=False,
    )

    selected = pd.DataFrame(
        {
            "module": sites[site_module].str.strip(),
            "chrom": sites[site_chrom].str.strip(),
            "pos": pd.to_numeric(
                sites[site_pos],
                errors="coerce",
            ),
            "ref": sites[site_ref].str.strip(),
            "alt": sites[site_alt].str.strip(),
            "minimum_comparative_BA": (
                pd.to_numeric(
                    sites[site_ba],
                    errors="coerce",
                )
                if site_ba is not None
                else np.nan
            ),
        }
    )

    selected = (
        selected
        .dropna(subset=["pos"])
        .drop_duplicates(
            [
                "module",
                "chrom",
                "pos",
                "ref",
                "alt",
            ]
        )
    )

    selected["pos"] = selected["pos"].astype(int)

    # ------------------------------------------------------------
    # Resolve genotype columns
    # ------------------------------------------------------------

    geno_module = find_column(
        genotypes,
        ["module", "source_module"],
    )

    geno_chrom = find_column(
        genotypes,
        ["chrom", "CHROM", "chromosome"],
    )

    geno_pos = find_column(
        genotypes,
        ["pos", "POS", "position"],
    )

    geno_ref = find_column(
        genotypes,
        ["ref", "REF"],
    )

    geno_alt = find_column(
        genotypes,
        ["alt", "ALT"],
    )

    geno_sample = find_column(
        genotypes,
        [
            "sample",
            "Sample",
            "accession",
            "Accession",
        ],
    )

    geno_gt = find_column(
        genotypes,
        [
            "GT",
            "GT_raw",
            "canonical_GT",
            "genotype",
        ],
        required=False,
    )

    geno_biological = find_column(
        genotypes,
        [
            "biological_GT",
            "biological_genotype",
            "biological_gt",
        ],
        required=False,
    )

    geno_ad = find_column(
        genotypes,
        ["AD", "AD_raw"],
        required=False,
    )

    geno_dp = find_column(
        genotypes,
        ["DP", "depth"],
        required=False,
    )

    geno_gq = find_column(
        genotypes,
        ["GQ", "genotype_quality"],
        required=False,
    )

    geno_qc = find_column(
        genotypes,
        [
            "DP2_GQ20_pass",
            "QC_pass_DP2_GQ20",
            "qc_pass",
        ],
        required=False,
    )

    calls = pd.DataFrame(
        {
            "module": genotypes[geno_module].str.strip(),
            "chrom": genotypes[geno_chrom].str.strip(),
            "pos": pd.to_numeric(
                genotypes[geno_pos],
                errors="coerce",
            ),
            "ref": genotypes[geno_ref].str.strip(),
            "alt": genotypes[geno_alt].str.strip(),
            "sample": genotypes[geno_sample].str.strip(),
            "GT": (
                genotypes[geno_gt]
                if geno_gt is not None
                else ""
            ),
            "biological_GT_input": (
                genotypes[geno_biological]
                if geno_biological is not None
                else ""
            ),
            "AD": (
                genotypes[geno_ad]
                if geno_ad is not None
                else ""
            ),
            "DP": (
                pd.to_numeric(
                    genotypes[geno_dp],
                    errors="coerce",
                )
                if geno_dp is not None
                else np.nan
            ),
            "GQ": (
                pd.to_numeric(
                    genotypes[geno_gq],
                    errors="coerce",
                )
                if geno_gq is not None
                else np.nan
            ),
            "QC_input": (
                genotypes[geno_qc]
                if geno_qc is not None
                else ""
            ),
        }
    )

    calls = calls.dropna(
        subset=["pos"]
    )

    calls["pos"] = calls["pos"].astype(int)

    calls = calls.merge(
        selected,
        on=[
            "module",
            "chrom",
            "pos",
            "ref",
            "alt",
        ],
        how="inner",
        validate="many_to_one",
    )

    calls["GT"] = calls["GT"].map(
        canonical_gt
    )

    calls["biological_GT"] = calls[
        "biological_GT_input"
    ].map(
        normalize_biological_gt
    )

    missing_bio = calls[
        "biological_GT"
    ].isna()

    calls.loc[
        missing_bio,
        "biological_GT",
    ] = calls.loc[
        missing_bio,
    ].apply(
        lambda row: biological_gt_from_raw(
            row["GT"],
            row["module"],
        ),
        axis=1,
    )

    calls["pattern_dosage"] = calls[
        "biological_GT"
    ].map(dosage)

    if geno_qc is not None:
        calls["DP2_GQ20_pass"] = calls[
            "QC_input"
        ].map(boolean)
    else:
        calls["DP2_GQ20_pass"] = (
            calls["GT"].notna()
            & calls["DP"].ge(2)
            & calls["GQ"].ge(20)
        )

    # ------------------------------------------------------------
    # Canonical metadata
    # ------------------------------------------------------------

    metadata_accession = find_column(
        metadata,
        ["Accession"],
    )

    metadata_taxon = find_column(
        metadata,
        ["Taxon Name"],
    )

    metadata_b = find_column(
        metadata,
        ["B"],
    )

    metadata_r = find_column(
        metadata,
        ["R"],
    )

    metadata_d = find_column(
        metadata,
        ["D"],
    )

    meta = pd.DataFrame(
        {
            "sample": metadata[
                metadata_accession
            ].str.strip(),
            "taxon": metadata[
                metadata_taxon
            ].map(clean_header),
            "B": metadata[metadata_b].map(binary),
            "R": metadata[metadata_r].map(binary),
            "D": metadata[metadata_d].map(binary),
        }
    )

    meta = meta.drop_duplicates(
        "sample"
    )

    calls = calls.merge(
        meta,
        on="sample",
        how="left",
        validate="many_to_one",
    )

    calls["piano_key"] = (
        calls["taxon"]
        .fillna("")
        .str.lower()
        .str.contains(
            "piano keys",
            regex=False,
        )
    )

    calls["discordant_group"] = ""

    calls.loc[
        (calls["D"] == 1)
        & (calls["R"] == 0),
        "discordant_group",
    ] = "D1R0"

    calls.loc[
        (calls["D"] == 0)
        & (calls["R"] == 1),
        "discordant_group",
    ] = "D0R1"

    # ------------------------------------------------------------
    # Evaluate each independently selected SNP
    # ------------------------------------------------------------

    site_keys = [
        "module",
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    results = []

    for key, site_calls in calls.groupby(
        site_keys,
        sort=False,
        dropna=False,
    ):
        module, chrom, pos, ref, alt = key

        site_calls = site_calls.copy()

        expectation = expected_focal_states(
            module
        )

        focal_lookup = {}

        for side, sample in {
            "Left": "TN_Left",
            "Right": "TN_Right",
        }.items():
            subset = site_calls[
                site_calls["sample"] == sample
            ]

            focal_lookup[side] = (
                subset.iloc[0]
                if len(subset) > 0
                else None
            )

        def focal_value(side):
            row = focal_lookup[side]

            if row is None:
                return {
                    "GT": "missing",
                    "biological_GT": "missing",
                    "dosage": np.nan,
                    "AD": "",
                    "DP": np.nan,
                    "GQ": np.nan,
                    "QC": False,
                }

            return {
                "GT": row["GT"] or "missing",
                "biological_GT": (
                    row["biological_GT"]
                    if pd.notna(
                        row["biological_GT"]
                    )
                    else "missing"
                ),
                "dosage": row["pattern_dosage"],
                "AD": row["AD"],
                "DP": row["DP"],
                "GQ": row["GQ"],
                "QC": bool(
                    row["DP2_GQ20_pass"]
                ),
            }

        left = focal_value("Left")
        right = focal_value("Right")

        focal_exact = (
            np.isfinite(left["dosage"])
            and np.isfinite(right["dosage"])
            and left["dosage"]
            == expectation["Left"]["dosage"]
            and right["dosage"]
            == expectation["Right"]["dosage"]
        )

        focal_exact_qc = (
            focal_exact
            and left["QC"]
            and right["QC"]
        )

        trait = comparative_trait(module)

        if trait in {"D", "R"}:
            discordant = site_calls[
                site_calls[
                    "discordant_group"
                ].isin(["D1R0", "D0R1"])
                & ~site_calls["piano_key"]
                & site_calls["DP2_GQ20_pass"]
                & site_calls[
                    "pattern_dosage"
                ].notna()
            ].copy()

            discordant[
                "predicted_presence"
            ] = (
                discordant[
                    "pattern_dosage"
                ] > 0
            )

            discordant[
                "expected_presence"
            ] = (
                discordant[trait] == 1
            )

            group_results = {}

            for group in [
                "D1R0",
                "D0R1",
            ]:
                group_rows = discordant[
                    discordant[
                        "discordant_group"
                    ] == group
                ]

                called = len(group_rows)

                correct = int(
                    (
                        group_rows[
                            "predicted_presence"
                        ]
                        == group_rows[
                            "expected_presence"
                        ]
                    ).sum()
                )

                errors = called - correct

                accuracy = safe_fraction(
                    correct,
                    called,
                )

                group_results[group] = {
                    "called": called,
                    "correct": correct,
                    "errors": errors,
                    "accuracy": accuracy,
                }

            both_groups_callable = (
                group_results["D1R0"][
                    "called"
                ] > 0
                and group_results["D0R1"][
                    "called"
                ] > 0
            )

            discordant_ba = (
                (
                    group_results["D1R0"][
                        "accuracy"
                    ]
                    + group_results["D0R1"][
                        "accuracy"
                    ]
                ) / 2
                if both_groups_callable
                else np.nan
            )

        else:
            group_results = {
                "D1R0": {
                    "called": 0,
                    "correct": 0,
                    "errors": 0,
                    "accuracy": np.nan,
                },
                "D0R1": {
                    "called": 0,
                    "correct": 0,
                    "errors": 0,
                    "accuracy": np.nan,
                },
            }

            both_groups_callable = True
            discordant_ba = np.nan

        specificity_ba80 = (
            True
            if trait == "B"
            else (
                both_groups_callable
                and np.isfinite(
                    discordant_ba
                )
                and discordant_ba >= 0.80
            )
        )

        minimum_ba = pd.to_numeric(
            site_calls[
                "minimum_comparative_BA"
            ],
            errors="coerce",
        ).max()

        results.append(
            {
                "focal_class": focal_class(module),
                "source_module": module,
                "comparative_trait": trait,
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "minimum_comparative_BA": minimum_ba,
                "selection_used_focal_data": False,
                "D1R0_called": group_results[
                    "D1R0"
                ]["called"],
                "D1R0_correct": group_results[
                    "D1R0"
                ]["correct"],
                "D1R0_errors": group_results[
                    "D1R0"
                ]["errors"],
                "D1R0_accuracy": group_results[
                    "D1R0"
                ]["accuracy"],
                "D0R1_called": group_results[
                    "D0R1"
                ]["called"],
                "D0R1_correct": group_results[
                    "D0R1"
                ]["correct"],
                "D0R1_errors": group_results[
                    "D0R1"
                ]["errors"],
                "D0R1_accuracy": group_results[
                    "D0R1"
                ]["accuracy"],
                "both_discordant_groups_callable": (
                    both_groups_callable
                ),
                "discordant_balanced_accuracy": (
                    discordant_ba
                ),
                "discordant_BA80": specificity_ba80,
                "left_expected_GT": expectation[
                    "Left"
                ]["GT"],
                "left_GT": left["GT"],
                "left_biological_GT": left[
                    "biological_GT"
                ],
                "left_pattern_dosage": left[
                    "dosage"
                ],
                "left_AD": left["AD"],
                "left_DP": left["DP"],
                "left_GQ": left["GQ"],
                "left_DP2_GQ20_pass": left[
                    "QC"
                ],
                "right_expected_GT": expectation[
                    "Right"
                ]["GT"],
                "right_GT": right["GT"],
                "right_biological_GT": right[
                    "biological_GT"
                ],
                "right_pattern_dosage": right[
                    "dosage"
                ],
                "right_AD": right["AD"],
                "right_DP": right["DP"],
                "right_GQ": right["GQ"],
                "right_DP2_GQ20_pass": right[
                    "QC"
                ],
                "exact_focal_zygosity_match": (
                    focal_exact
                ),
                "exact_focal_zygosity_DP2_GQ20": (
                    focal_exact_qc
                ),
                "exact_and_discordant_BA80": (
                    focal_exact_qc
                    and specificity_ba80
                ),
            }
        )

    result = pd.DataFrame(results)

    result = result.sort_values(
        [
            "focal_class",
            "source_module",
            "pos",
        ]
    ).reset_index(drop=True)

    exact = result[
        result[
            "exact_focal_zygosity_match"
        ]
    ].copy()

    exact_qc = result[
        result[
            "exact_focal_zygosity_DP2_GQ20"
        ]
    ].copy()

    exact_specific = result[
        result[
            "exact_and_discordant_BA80"
        ]
    ].copy()

    result.to_csv(
        OUT_ALL,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    exact.to_csv(
        OUT_EXACT,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    exact_qc.to_csv(
        OUT_EXACT_QC,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    exact_specific.to_csv(
        OUT_EXACT_SPECIFIC,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary_rows = []

    for module in [
        "Band",
        "Dennis",
        "Ray",
    ]:
        subset = result[
            result["source_module"] == module
        ]

        summary_rows.append(
            {
                "focal_class": (
                    "B"
                    if module == "Band"
                    else "DR"
                ),
                "source_module": module,
                "comparative_selected_SNPs": len(
                    subset
                ),
                "exact_focal_zygosity": int(
                    subset[
                        "exact_focal_zygosity_match"
                    ].sum()
                ),
                "exact_focal_zygosity_DP2_GQ20": int(
                    subset[
                        "exact_focal_zygosity_DP2_GQ20"
                    ].sum()
                ),
                "discordant_BA80": int(
                    subset["discordant_BA80"].sum()
                ),
                "exact_focal_and_discordant_BA80": int(
                    subset[
                        "exact_and_discordant_BA80"
                    ].sum()
                ),
            }
        )

    summary = pd.DataFrame(
        summary_rows
    )

    summary.to_csv(
        OUT_SUMMARY,
        sep="\t",
        index=False,
    )

    print(
        "CANONICAL B VERSUS DR ZYGOSITY EVALUATION"
    )
    print()
    print(
        summary.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print()
    print(
        "EXACT FOCAL SITES PASSING RECIPROCAL "
        "D1R0/D0R1 BA >= 0.80"
    )
    print(
        "class\tmodule\tposition\tcomparative_BA\t"
        "discordant_BA\tD1R0\tD0R1\tleft\tright"
    )

    for _, row in exact_specific.iterrows():
        comparative_ba = row[
            "minimum_comparative_BA"
        ]

        discordant_score = row[
            "discordant_balanced_accuracy"
        ]

        comparative_text = (
            f"{comparative_ba:.4f}"
            if np.isfinite(comparative_ba)
            else "NA"
        )

        discordant_text = (
            f"{discordant_score:.4f}"
            if np.isfinite(discordant_score)
            else "NA"
        )

        print(
            f'{row["focal_class"]}\t'
            f'{row["source_module"]}\t'
            f'{row["chrom"]}:{row["pos"]}\t'
            f'{comparative_text}\t'
            f'{discordant_text}\t'
            f'{row["D1R0_correct"]}/'
            f'{row["D1R0_called"]}\t'
            f'{row["D0R1_correct"]}/'
            f'{row["D0R1_called"]}\t'
            f'{row["left_biological_GT"]}\t'
            f'{row["right_biological_GT"]}'
        )

    print()
    print(f"All evaluations: {OUT_ALL}")
    print(f"All exact focal matches: {OUT_EXACT}")
    print(f"QC exact focal matches: {OUT_EXACT_QC}")
    print(
        "QC exact plus reciprocal specificity: "
        f"{OUT_EXACT_SPECIFIC}"
    )
    print(f"Summary: {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
