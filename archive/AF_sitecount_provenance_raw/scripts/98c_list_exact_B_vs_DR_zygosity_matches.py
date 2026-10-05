#!/usr/bin/env python3

from pathlib import Path
import math
import re

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

GENOTYPE_TABLE = (
    CHECKPOINT_DIR
    / "core_explanatory_optix_genotypes.long.tsv"
)

OUTPUT_ALL_EVALUATED = (
    CHECKPOINT_DIR
    / "B_vs_DR_focal_zygosity_evaluation.tsv"
)

OUTPUT_EXACT_CALLED = (
    CHECKPOINT_DIR
    / "B_vs_DR_exact_zygosity_matches.all_called.tsv"
)

OUTPUT_EXACT_QC = (
    CHECKPOINT_DIR
    / "B_vs_DR_exact_zygosity_matches.DP2_GQ20.tsv"
)

OUTPUT_SUMMARY = (
    CHECKPOINT_DIR
    / "B_vs_DR_exact_zygosity_matches.summary.tsv"
)


def find_column(dataframe, candidates, required=True):
    lookup = {
        str(column).strip().lower(): column
        for column in dataframe.columns
    }

    for candidate in candidates:
        key = candidate.strip().lower()

        if key in lookup:
            return lookup[key]

    if required:
        raise RuntimeError(
            "Could not identify a required column.\n"
            f"Tried: {', '.join(candidates)}\n"
            f"Available: {', '.join(map(str, dataframe.columns))}"
        )

    return None


def as_bool(value):
    if pd.isna(value):
        return False

    return str(value).strip().lower() in {
        "true",
        "t",
        "1",
        "yes",
        "y",
    }


def as_number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return math.nan

    if not math.isfinite(result):
        return math.nan

    return result


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
        or any(allele in {"", "."} for allele in alleles)
    ):
        return None

    try:
        numeric = sorted(int(allele) for allele in alleles)
    except ValueError:
        return None

    return "/".join(str(allele) for allele in numeric)


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

    if text in {"AA", "A0A0"}:
        return "A/A"

    if text in {"PA", "AP"}:
        return "P/A"

    if text in {"PP"}:
        return "P/P"

    return None


def biological_gt_from_vcf(gt, focal_class):
    """
    Biological polarity:

      B:
        HMR reference is P/P.
        0/0 -> P/P
        0/1 -> P/A
        1/1 -> A/A

      DR:
        HMR reference is A/A.
        0/0 -> A/A
        0/1 -> P/A
        1/1 -> P/P
    """
    gt = canonical_gt(gt)

    if gt is None:
        return None

    if focal_class == "B":
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


def genotype_dosage(biological_gt):
    return {
        "A/A": 0.0,
        "P/A": 0.5,
        "P/P": 1.0,
    }.get(biological_gt, math.nan)


def classify_module(module):
    text = str(module).strip().lower()

    if "band" in text:
        return "B"

    if "ray" in text or "dennis" in text:
        return "DR"

    return None


def expected_states(focal_class):
    if focal_class == "B":
        return {
            "left_gt": "A/A",
            "left_dosage": 0.0,
            "right_gt": "P/A",
            "right_dosage": 0.5,
        }

    if focal_class == "DR":
        return {
            "left_gt": "P/P",
            "left_dosage": 1.0,
            "right_gt": "P/A",
            "right_dosage": 0.5,
        }

    raise RuntimeError(
        f"Unexpected focal class: {focal_class}"
    )


def main():
    if not GENOTYPE_TABLE.is_file():
        raise FileNotFoundError(
            f"Missing existing genotype table: {GENOTYPE_TABLE}"
        )

    data = pd.read_csv(
        GENOTYPE_TABLE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    module_col = find_column(
        data,
        ["module", "trait_module", "source_module"],
    )

    chrom_col = find_column(
        data,
        ["chrom", "CHROM", "chromosome"],
    )

    pos_col = find_column(
        data,
        ["pos", "POS", "position"],
    )

    ref_col = find_column(
        data,
        ["ref", "REF"],
    )

    alt_col = find_column(
        data,
        ["alt", "ALT"],
    )

    sample_col = find_column(
        data,
        ["sample", "Sample", "accession", "Accession"],
    )

    gt_col = find_column(
        data,
        ["GT", "GT_raw", "canonical_GT", "genotype"],
        required=False,
    )

    biological_gt_col = find_column(
        data,
        [
            "biological_GT",
            "biological_genotype",
            "biological_gt",
        ],
        required=False,
    )

    ad_col = find_column(
        data,
        ["AD", "AD_raw", "allelic_depth"],
        required=False,
    )

    dp_col = find_column(
        data,
        ["DP", "depth"],
        required=False,
    )

    gq_col = find_column(
        data,
        ["GQ", "genotype_quality"],
        required=False,
    )

    qc_col = find_column(
        data,
        [
            "DP2_GQ20_pass",
            "QC_pass_DP2_GQ20",
            "qc_pass",
        ],
        required=False,
    )

    pattern_fraction_col = find_column(
        data,
        [
            "pattern_allele_fraction",
            "presence_allele_fraction",
            "P_allele_fraction",
        ],
        required=False,
    )

    alt_fraction_col = find_column(
        data,
        ["ALT_fraction", "ALT_fraction_raw"],
        required=False,
    )

    comparative_ba_col = find_column(
        data,
        [
            "minimum_BA",
            "minimum_comparative_BA",
            "minimum_BA_all_cohorts_schemes",
        ],
        required=False,
    )

    focal = data[
        data[sample_col].isin(
            ["TN_Left", "TN_Right"]
        )
    ].copy()

    focal["focal_class"] = focal[module_col].map(
        classify_module
    )

    focal = focal[
        focal["focal_class"].isin(["B", "DR"])
    ].copy()

    focal["side"] = focal[sample_col].map(
        {
            "TN_Left": "Left",
            "TN_Right": "Right",
        }
    )

    focal["_pos_numeric"] = pd.to_numeric(
        focal[pos_col],
        errors="coerce",
    )

    site_columns = [
        module_col,
        chrom_col,
        pos_col,
        ref_col,
        alt_col,
    ]

    rows = []

    grouped = focal.groupby(
        site_columns,
        sort=False,
        dropna=False,
    )

    for site_key, site_data in grouped:
        source_module, chrom, pos, ref, alt = site_key
        focal_class = classify_module(source_module)

        if focal_class is None:
            continue

        side_rows = {}

        for side in ["Left", "Right"]:
            subset = site_data[
                site_data["side"] == side
            ]

            if subset.empty:
                side_rows[side] = None
            else:
                side_rows[side] = subset.iloc[0]

        left = side_rows["Left"]
        right = side_rows["Right"]

        expectation = expected_states(focal_class)

        def extract_side(row, side):
            if row is None:
                return {
                    "GT": "missing",
                    "biological_gt": None,
                    "dosage": math.nan,
                    "AD": "",
                    "DP": math.nan,
                    "GQ": math.nan,
                    "QC": False,
                    "pattern_fraction": math.nan,
                }

            raw_gt = (
                canonical_gt(row[gt_col])
                if gt_col is not None
                else None
            )

            biological_gt = None

            if biological_gt_col is not None:
                biological_gt = normalize_biological_gt(
                    row[biological_gt_col]
                )

            if biological_gt is None:
                biological_gt = biological_gt_from_vcf(
                    raw_gt,
                    focal_class,
                )

            dosage = genotype_dosage(
                biological_gt
            )

            dp = (
                as_number(row[dp_col])
                if dp_col is not None
                else math.nan
            )

            gq = (
                as_number(row[gq_col])
                if gq_col is not None
                else math.nan
            )

            if qc_col is not None:
                qc_pass = as_bool(row[qc_col])
            else:
                qc_pass = (
                    raw_gt is not None
                    and math.isfinite(dp)
                    and math.isfinite(gq)
                    and dp >= 2
                    and gq >= 20
                )

            pattern_fraction = math.nan

            if pattern_fraction_col is not None:
                pattern_fraction = as_number(
                    row[pattern_fraction_col]
                )

            elif alt_fraction_col is not None:
                alt_fraction = as_number(
                    row[alt_fraction_col]
                )

                if math.isfinite(alt_fraction):
                    if focal_class == "B":
                        pattern_fraction = 1.0 - alt_fraction
                    else:
                        pattern_fraction = alt_fraction

            return {
                "GT": raw_gt or "missing",
                "biological_gt": biological_gt or "missing",
                "dosage": dosage,
                "AD": (
                    row[ad_col]
                    if ad_col is not None
                    else ""
                ),
                "DP": dp,
                "GQ": gq,
                "QC": qc_pass,
                "pattern_fraction": pattern_fraction,
            }

        left_call = extract_side(left, "Left")
        right_call = extract_side(right, "Right")

        both_called = (
            math.isfinite(left_call["dosage"])
            and math.isfinite(right_call["dosage"])
        )

        if both_called:
            left_distance = abs(
                left_call["dosage"]
                - expectation["left_dosage"]
            )

            right_distance = abs(
                right_call["dosage"]
                - expectation["right_dosage"]
            )

            pair_distance = (
                left_distance + right_distance
            )

            # Transparent normalized score:
            # 1 = exact pair; 0 = maximum possible total deviation.
            pair_score = 1.0 - (
                pair_distance / 2.0
            )

            exact_match = (
                abs(pair_distance) < 1e-12
            )
        else:
            left_distance = math.nan
            right_distance = math.nan
            pair_distance = math.nan
            pair_score = math.nan
            exact_match = False

        comparative_ba = math.nan

        if comparative_ba_col is not None:
            comparative_ba = as_number(
                site_data.iloc[0][comparative_ba_col]
            )

        rows.append(
            {
                "focal_class": focal_class,
                "source_module": source_module,
                "chrom": chrom,
                "pos": int(float(pos)),
                "ref": ref,
                "alt": alt,
                "comparative_minimum_BA": comparative_ba,
                "selection_used_focal_data": False,
                "left_expected_GT": expectation["left_gt"],
                "left_expected_dosage": expectation[
                    "left_dosage"
                ],
                "left_GT": left_call["GT"],
                "left_biological_GT": left_call[
                    "biological_gt"
                ],
                "left_observed_dosage": left_call[
                    "dosage"
                ],
                "left_pattern_allele_fraction": left_call[
                    "pattern_fraction"
                ],
                "left_AD": left_call["AD"],
                "left_DP": left_call["DP"],
                "left_GQ": left_call["GQ"],
                "left_DP2_GQ20_pass": left_call["QC"],
                "left_zygosity_distance": left_distance,
                "right_expected_GT": expectation["right_gt"],
                "right_expected_dosage": expectation[
                    "right_dosage"
                ],
                "right_GT": right_call["GT"],
                "right_biological_GT": right_call[
                    "biological_gt"
                ],
                "right_observed_dosage": right_call[
                    "dosage"
                ],
                "right_pattern_allele_fraction": right_call[
                    "pattern_fraction"
                ],
                "right_AD": right_call["AD"],
                "right_DP": right_call["DP"],
                "right_GQ": right_call["GQ"],
                "right_DP2_GQ20_pass": right_call["QC"],
                "right_zygosity_distance": right_distance,
                "both_called": both_called,
                "both_DP2_GQ20_pass": (
                    left_call["QC"]
                    and right_call["QC"]
                ),
                "pair_zygosity_distance": pair_distance,
                "pair_zygosity_score": pair_score,
                "exact_expected_zygosity_match": exact_match,
            }
        )

    evaluation = pd.DataFrame(rows)

    if evaluation.empty:
        raise RuntimeError(
            "No B or DR focal sites were found in the existing "
            "genotype table."
        )

    evaluation = evaluation.sort_values(
        [
            "focal_class",
            "source_module",
            "chrom",
            "pos",
        ]
    ).reset_index(drop=True)

    exact_called = evaluation[
        evaluation[
            "exact_expected_zygosity_match"
        ]
    ].copy()

    exact_qc = exact_called[
        exact_called[
            "both_DP2_GQ20_pass"
        ]
    ].copy()

    evaluation.to_csv(
        OUTPUT_ALL_EVALUATED,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    exact_called.to_csv(
        OUTPUT_EXACT_CALLED,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    exact_qc.to_csv(
        OUTPUT_EXACT_QC,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary_rows = []

    for focal_class in ["B", "DR"]:
        class_rows = evaluation[
            evaluation["focal_class"] == focal_class
        ]

        class_exact = exact_called[
            exact_called["focal_class"] == focal_class
        ]

        class_exact_qc = exact_qc[
            exact_qc["focal_class"] == focal_class
        ]

        summary_rows.append(
            {
                "focal_class": focal_class,
                "comparative_selected_SNPs": len(
                    class_rows
                ),
                "both_focal_sides_called": int(
                    class_rows["both_called"].sum()
                ),
                "exact_zygosity_matches_all_called": len(
                    class_exact
                ),
                "exact_zygosity_matches_DP2_GQ20": len(
                    class_exact_qc
                ),
                "expected_left": (
                    "A/A"
                    if focal_class == "B"
                    else "P/P"
                ),
                "expected_right": "P/A",
            }
        )

    summary = pd.DataFrame(summary_rows)

    summary.to_csv(
        OUTPUT_SUMMARY,
        sep="\t",
        index=False,
    )

    print("EXACT B VERSUS DR FOCAL ZYGOSITY MATCHES")
    print()
    print(
        "The starting SNP set was selected using comparative "
        "butterflies only."
    )
    print(
        "Focal genotypes were evaluated only after site "
        "ascertainment."
    )
    print()
    print("EXPECTED PATTERNS")
    print("B:  Left A/A (0); Right P/A (0.5)")
    print("DR: Left P/P (1); Right P/A (0.5)")
    print()
    print("SUMMARY")
    print(
        summary.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )
    print()
    print("DP2/GQ20 EXACT-MATCH SITES")
    print(
        "class\tsource_module\tposition\tcomparative_BA\t"
        "left_GT\tright_GT\tleft_AD\tright_AD"
    )

    for _, row in exact_qc.iterrows():
        comparative_ba = row[
            "comparative_minimum_BA"
        ]

        if pd.isna(comparative_ba):
            ba_text = "NA"
        else:
            ba_text = f"{comparative_ba:.4f}"

        print(
            f'{row["focal_class"]}\t'
            f'{row["source_module"]}\t'
            f'{row["chrom"]}:{row["pos"]}\t'
            f'{ba_text}\t'
            f'{row["left_biological_GT"]}\t'
            f'{row["right_biological_GT"]}\t'
            f'{row["left_AD"]}\t'
            f'{row["right_AD"]}'
        )

    print()
    print(f"All evaluated sites: {OUTPUT_ALL_EVALUATED}")
    print(f"All called exact matches: {OUTPUT_EXACT_CALLED}")
    print(f"QC-passing exact matches: {OUTPUT_EXACT_QC}")
    print(f"Summary: {OUTPUT_SUMMARY}")


if __name__ == "__main__":
    main()
