#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
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
    / "B_DR_three_way_physical_SNP_classification.tsv"
)

OUT_EXACT = (
    CHECKPOINT_DIR
    / "B_DR_three_way_exact_focal_zygosity.tsv"
)

OUT_FINAL = (
    CHECKPOINT_DIR
    / "B_DR_three_way_high_confidence_exact_focal.tsv"
)

OUT_AMBIGUOUS = (
    CHECKPOINT_DIR
    / "B_DR_three_way_ambiguous_exact_focal.tsv"
)

OUT_SUMMARY = (
    CHECKPOINT_DIR
    / "B_DR_three_way_classification.summary.tsv"
)

MINIMUM_GROUP_CALL_FRACTION = 0.50
MINIMUM_MACRO_ACCURACY = 0.80
MINIMUM_SCORE_MARGIN = 0.10

GROUPS = [
    "besckei_B",
    "meriana_D",
    "contigua_R",
]

HYPOTHESES = [
    "B",
    "D",
    "R",
]

EXPECTED_CARRIER = {
    "B": {
        "besckei_B": True,
        "meriana_D": False,
        "contigua_R": False,
    },
    "D": {
        "besckei_B": False,
        "meriana_D": True,
        "contigua_R": False,
    },
    "R": {
        "besckei_B": False,
        "meriana_D": False,
        "contigua_R": True,
    },
}


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


def find_ba_column(dataframe):
    direct = find_column(
        dataframe,
        [
            "minimum_BA",
            "minimum_comparative_BA",
            "minimum_BA_all_cohorts_schemes",
            "minimum_BA_across_cohorts_and_schemes",
            "strict_minimum_BA",
            "minimum_BA_across_schemes_and_cohorts",
        ],
        required=False,
    )

    if direct is not None:
        return direct

    candidates = []

    for column in dataframe.columns:
        cleaned = clean_header(column).lower()

        if (
            "ba" in cleaned
            and (
                "minimum" in cleaned
                or "min_" in cleaned
                or cleaned.startswith("min")
            )
        ):
            candidates.append(column)

    if len(candidates) == 1:
        return candidates[0]

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


def binary(value):
    result = pd.to_numeric(
        value,
        errors="coerce",
    )

    if pd.isna(result):
        return np.nan

    result = float(result)

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
        alleles = sorted(
            int(allele)
            for allele in alleles
        )
    except ValueError:
        return None

    if any(
        allele not in {0, 1}
        for allele in alleles
    ):
        return None

    return "/".join(
        str(allele)
        for allele in alleles
    )


def raw_alt_dosage(gt):
    gt = canonical_gt(gt)

    if gt is None:
        return np.nan

    alleles = [
        int(allele)
        for allele in gt.split("/")
    ]

    return sum(alleles) / 2.0


def biological_dosage(gt, hypothesis):
    alt_dosage = raw_alt_dosage(gt)

    if not np.isfinite(alt_dosage):
        return np.nan

    if hypothesis == "B":
        # HMR reference carries Band.
        return 1.0 - alt_dosage

    # HMR reference lacks Dennis and Ray.
    return alt_dosage


def biological_gt(dosage):
    if not np.isfinite(dosage):
        return "missing"

    if abs(dosage) < 1e-9:
        return "A/A"

    if abs(dosage - 0.5) < 1e-9:
        return "P/A"

    if abs(dosage - 1.0) < 1e-9:
        return "P/P"

    return "other"


def hypothesis_from_module(module):
    text = str(module).strip().lower()

    if text == "band":
        return "B"

    if text == "dennis":
        return "D"

    if text == "ray":
        return "R"

    raise RuntimeError(
        f"Unexpected comparative module: {module}"
    )


def explanatory_group(row):
    b = row["B"]
    d = row["D"]
    r = row["R"]

    if (b, d, r) == (1, 0, 0):
        return "besckei_B"

    if (b, d, r) == (0, 1, 0):
        return "meriana_D"

    if (b, d, r) == (0, 0, 1):
        return "contigua_R"

    return ""


def safe_accuracy(correct, called):
    if called <= 0:
        return np.nan

    return correct / called


def focal_expectation(hypothesis):
    if hypothesis == "B":
        return {
            "Left": 0.0,
            "Right": 0.5,
        }

    return {
        "Left": 1.0,
        "Right": 0.5,
    }


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

    selected_raw = pd.read_csv(
        SELECTED_SITES,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    genotype_raw = pd.read_csv(
        GENOTYPES,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    metadata_raw = pd.read_csv(
        METADATA,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    metadata_raw.columns = [
        clean_header(column)
        for column in metadata_raw.columns
    ]

    # ------------------------------------------------------------
    # Comparative-selected hypothesis table
    # ------------------------------------------------------------

    s_module = find_column(
        selected_raw,
        ["module", "source_module"],
    )

    s_chrom = find_column(
        selected_raw,
        ["chrom", "CHROM", "chromosome"],
    )

    s_pos = find_column(
        selected_raw,
        ["pos", "POS", "position"],
    )

    s_ref = find_column(
        selected_raw,
        ["ref", "REF"],
    )

    s_alt = find_column(
        selected_raw,
        ["alt", "ALT"],
    )

    s_ba = find_ba_column(
        selected_raw
    )

    selected = pd.DataFrame(
        {
            "hypothesis": selected_raw[
                s_module
            ].map(hypothesis_from_module),
            "chrom": selected_raw[
                s_chrom
            ].str.strip(),
            "pos": pd.to_numeric(
                selected_raw[s_pos],
                errors="coerce",
            ),
            "ref": selected_raw[
                s_ref
            ].str.strip(),
            "alt": selected_raw[
                s_alt
            ].str.strip(),
            "comparative_BA": (
                pd.to_numeric(
                    selected_raw[s_ba],
                    errors="coerce",
                )
                if s_ba is not None
                else np.nan
            ),
        }
    )

    selected = selected.dropna(
        subset=["pos"]
    )

    selected["pos"] = selected[
        "pos"
    ].astype(int)

    selected = selected.drop_duplicates(
        [
            "hypothesis",
            "chrom",
            "pos",
            "ref",
            "alt",
        ]
    )

    physical_keys = [
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    supported_hypotheses = defaultdict(set)
    comparative_ba = {}

    for _, row in selected.iterrows():
        key = tuple(
            row[column]
            for column in physical_keys
        )

        hypothesis = row["hypothesis"]

        supported_hypotheses[key].add(
            hypothesis
        )

        comparative_ba[
            (key, hypothesis)
        ] = row["comparative_BA"]

    physical_sites = (
        selected[physical_keys]
        .drop_duplicates()
        .copy()
    )

    # ------------------------------------------------------------
    # Raw genotype calls
    # ------------------------------------------------------------

    g_module = find_column(
        genotype_raw,
        ["module", "source_module"],
    )

    g_chrom = find_column(
        genotype_raw,
        ["chrom", "CHROM", "chromosome"],
    )

    g_pos = find_column(
        genotype_raw,
        ["pos", "POS", "position"],
    )

    g_ref = find_column(
        genotype_raw,
        ["ref", "REF"],
    )

    g_alt = find_column(
        genotype_raw,
        ["alt", "ALT"],
    )

    g_sample = find_column(
        genotype_raw,
        [
            "sample",
            "Sample",
            "accession",
            "Accession",
        ],
    )

    g_gt = find_column(
        genotype_raw,
        [
            "GT",
            "GT_raw",
            "canonical_GT",
            "genotype",
        ],
    )

    g_ad = find_column(
        genotype_raw,
        ["AD", "AD_raw"],
        required=False,
    )

    g_dp = find_column(
        genotype_raw,
        ["DP", "depth"],
        required=False,
    )

    g_gq = find_column(
        genotype_raw,
        ["GQ", "genotype_quality"],
        required=False,
    )

    g_qc = find_column(
        genotype_raw,
        [
            "DP2_GQ20_pass",
            "QC_pass_DP2_GQ20",
            "qc_pass",
        ],
        required=False,
    )

    calls = pd.DataFrame(
        {
            "source_module": genotype_raw[
                g_module
            ].str.strip(),
            "chrom": genotype_raw[
                g_chrom
            ].str.strip(),
            "pos": pd.to_numeric(
                genotype_raw[g_pos],
                errors="coerce",
            ),
            "ref": genotype_raw[
                g_ref
            ].str.strip(),
            "alt": genotype_raw[
                g_alt
            ].str.strip(),
            "sample": genotype_raw[
                g_sample
            ].str.strip(),
            "GT": genotype_raw[g_gt].map(
                canonical_gt
            ),
            "AD": (
                genotype_raw[g_ad]
                if g_ad is not None
                else ""
            ),
            "DP": (
                pd.to_numeric(
                    genotype_raw[g_dp],
                    errors="coerce",
                )
                if g_dp is not None
                else np.nan
            ),
            "GQ": (
                pd.to_numeric(
                    genotype_raw[g_gq],
                    errors="coerce",
                )
                if g_gq is not None
                else np.nan
            ),
            "QC_input": (
                genotype_raw[g_qc]
                if g_qc is not None
                else ""
            ),
        }
    )

    calls = calls.dropna(
        subset=["pos"]
    )

    calls["pos"] = calls[
        "pos"
    ].astype(int)

    calls = calls.merge(
        physical_sites,
        on=physical_keys,
        how="inner",
    )

    if g_qc is not None:
        calls["DP2_GQ20_pass"] = calls[
            "QC_input"
        ].map(boolean)
    else:
        calls["DP2_GQ20_pass"] = (
            calls["GT"].notna()
            & calls["DP"].ge(2)
            & calls["GQ"].ge(20)
        )

    # The same physical SNP may have entered through multiple
    # comparative module hypotheses. Raw calls must agree.
    consistency = (
        calls
        .dropna(subset=["GT"])
        .groupby(
            physical_keys + ["sample"]
        )["GT"]
        .nunique()
    )

    inconsistent = int(
        (consistency > 1).sum()
    )

    if inconsistent > 0:
        raise RuntimeError(
            f"{inconsistent} physical-site/sample combinations "
            "have inconsistent raw GTs across module copies."
        )

    calls = (
        calls
        .sort_values(
            physical_keys
            + ["sample", "source_module"]
        )
        .drop_duplicates(
            physical_keys + ["sample"],
            keep="first",
        )
    )

    # ------------------------------------------------------------
    # Canonical metadata and explanatory groups
    # ------------------------------------------------------------

    m_sample = find_column(
        metadata_raw,
        ["Accession"],
    )

    m_taxon = find_column(
        metadata_raw,
        ["Taxon Name"],
    )

    m_b = find_column(
        metadata_raw,
        ["B"],
    )

    m_d = find_column(
        metadata_raw,
        ["D"],
    )

    m_r = find_column(
        metadata_raw,
        ["R"],
    )

    metadata = pd.DataFrame(
        {
            "sample": metadata_raw[
                m_sample
            ].str.strip(),
            "taxon": metadata_raw[
                m_taxon
            ].map(clean_header),
            "B": metadata_raw[m_b].map(
                binary
            ),
            "D": metadata_raw[m_d].map(
                binary
            ),
            "R": metadata_raw[m_r].map(
                binary
            ),
        }
    )

    metadata = metadata.drop_duplicates(
        "sample"
    )

    calls = calls.merge(
        metadata,
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

    calls["explanatory_group"] = calls.apply(
        explanatory_group,
        axis=1,
    )

    group_totals = (
        calls[
            calls["explanatory_group"].isin(
                GROUPS
            )
            & ~calls["piano_key"]
        ][
            [
                "sample",
                "explanatory_group",
            ]
        ]
        .drop_duplicates()
        .groupby(
            "explanatory_group"
        )
        .size()
        .to_dict()
    )

    for group in GROUPS:
        if group_totals.get(group, 0) == 0:
            raise RuntimeError(
                f"No samples found for explanatory group {group}"
            )

    minimum_group_calls = {
        group: max(
            1,
            math.ceil(
                group_totals[group]
                * MINIMUM_GROUP_CALL_FRACTION
            ),
        )
        for group in GROUPS
    }

    # ------------------------------------------------------------
    # Three-way classification
    # ------------------------------------------------------------

    result_rows = []

    for key, site_calls in calls.groupby(
        physical_keys,
        sort=False,
        dropna=False,
    ):
        key = tuple(key)
        chrom, pos, ref, alt = key

        eligible_hypotheses = (
            supported_hypotheses[key]
        )

        hypothesis_results = {}

        for hypothesis in HYPOTHESES:
            test_calls = site_calls[
                site_calls[
                    "explanatory_group"
                ].isin(GROUPS)
                & ~site_calls["piano_key"]
                & site_calls["DP2_GQ20_pass"]
                & site_calls["GT"].notna()
            ].copy()

            test_calls["pattern_dosage"] = (
                test_calls["GT"].map(
                    lambda gt: biological_dosage(
                        gt,
                        hypothesis,
                    )
                )
            )

            test_calls["predicted_carrier"] = (
                test_calls[
                    "pattern_dosage"
                ] > 0
            )

            group_metrics = {}

            for group in GROUPS:
                group_calls = test_calls[
                    test_calls[
                        "explanatory_group"
                    ] == group
                ]

                expected = EXPECTED_CARRIER[
                    hypothesis
                ][group]

                called = len(group_calls)

                correct = int(
                    (
                        group_calls[
                            "predicted_carrier"
                        ]
                        == expected
                    ).sum()
                )

                accuracy = safe_accuracy(
                    correct,
                    called,
                )

                group_metrics[group] = {
                    "called": called,
                    "correct": correct,
                    "errors": called - correct,
                    "accuracy": accuracy,
                }

            callability_pass = all(
                group_metrics[group]["called"]
                >= minimum_group_calls[group]
                for group in GROUPS
            )

            macro_accuracy = (
                np.mean(
                    [
                        group_metrics[group][
                            "accuracy"
                        ]
                        for group in GROUPS
                    ]
                )
                if all(
                    np.isfinite(
                        group_metrics[group][
                            "accuracy"
                        ]
                    )
                    for group in GROUPS
                )
                else np.nan
            )

            hypothesis_results[hypothesis] = {
                "groups": group_metrics,
                "callability_pass": callability_pass,
                "macro_accuracy": macro_accuracy,
                "comparative_supported": (
                    hypothesis
                    in eligible_hypotheses
                ),
                "comparative_BA": comparative_ba.get(
                    (key, hypothesis),
                    np.nan,
                ),
            }

        ranked = sorted(
            HYPOTHESES,
            key=lambda hypothesis: (
                hypothesis_results[hypothesis][
                    "macro_accuracy"
                ]
                if np.isfinite(
                    hypothesis_results[
                        hypothesis
                    ]["macro_accuracy"]
                )
                else -1.0
            ),
            reverse=True,
        )

        best = ranked[0]
        second = ranked[1]

        best_score = hypothesis_results[
            best
        ]["macro_accuracy"]

        second_score = hypothesis_results[
            second
        ]["macro_accuracy"]

        score_margin = (
            best_score - second_score
            if np.isfinite(best_score)
            and np.isfinite(second_score)
            else np.nan
        )

        tied = (
            np.isfinite(score_margin)
            and abs(score_margin) < 1e-12
        )

        best_supported = hypothesis_results[
            best
        ]["comparative_supported"]

        best_callability = hypothesis_results[
            best
        ]["callability_pass"]

        classification_pass = (
            best_callability
            and best_supported
            and np.isfinite(best_score)
            and best_score
            >= MINIMUM_MACRO_ACCURACY
            and np.isfinite(score_margin)
            and score_margin
            >= MINIMUM_SCORE_MARGIN
            and not tied
        )

        focal_rows = {
            side: site_calls[
                site_calls["sample"] == sample
            ]
            for side, sample in {
                "Left": "TN_Left",
                "Right": "TN_Right",
            }.items()
        }

        focal_values = {}

        for side in ["Left", "Right"]:
            subset = focal_rows[side]

            if subset.empty:
                focal_values[side] = {
                    "GT": "missing",
                    "AD": "",
                    "DP": np.nan,
                    "GQ": np.nan,
                    "QC": False,
                    "dosage": np.nan,
                    "biological_GT": "missing",
                }
            else:
                row = subset.iloc[0]

                focal_dosage = biological_dosage(
                    row["GT"],
                    best,
                )

                focal_values[side] = {
                    "GT": (
                        row["GT"]
                        if row["GT"] is not None
                        else "missing"
                    ),
                    "AD": row["AD"],
                    "DP": row["DP"],
                    "GQ": row["GQ"],
                    "QC": bool(
                        row["DP2_GQ20_pass"]
                    ),
                    "dosage": focal_dosage,
                    "biological_GT": biological_gt(
                        focal_dosage
                    ),
                }

        expectation = focal_expectation(
            best
        )

        focal_exact = (
            np.isfinite(
                focal_values["Left"]["dosage"]
            )
            and np.isfinite(
                focal_values["Right"]["dosage"]
            )
            and focal_values["Left"]["dosage"]
            == expectation["Left"]
            and focal_values["Right"]["dosage"]
            == expectation["Right"]
        )

        focal_exact_qc = (
            focal_exact
            and focal_values["Left"]["QC"]
            and focal_values["Right"]["QC"]
        )

        output = {
            "chrom": chrom,
            "pos": int(pos),
            "ref": ref,
            "alt": alt,
            "comparative_supported_hypotheses": ",".join(
                sorted(eligible_hypotheses)
            ),
            "best_hypothesis": best,
            "final_class": (
                "B"
                if best == "B"
                else "DR"
            ),
            "second_hypothesis": second,
            "best_macro_accuracy": best_score,
            "second_macro_accuracy": second_score,
            "score_margin": score_margin,
            "best_hypothesis_comparative_supported": (
                best_supported
            ),
            "best_hypothesis_comparative_BA": (
                hypothesis_results[best][
                    "comparative_BA"
                ]
            ),
            "core_callability_pass": (
                best_callability
            ),
            "three_way_classification_pass": (
                classification_pass
            ),
            "focal_expected_left": (
                "A/A"
                if best == "B"
                else "P/P"
            ),
            "focal_left_GT": focal_values[
                "Left"
            ]["GT"],
            "focal_left_biological_GT": (
                focal_values["Left"][
                    "biological_GT"
                ]
            ),
            "focal_left_AD": focal_values[
                "Left"
            ]["AD"],
            "focal_left_DP": focal_values[
                "Left"
            ]["DP"],
            "focal_left_GQ": focal_values[
                "Left"
            ]["GQ"],
            "focal_left_QC": focal_values[
                "Left"
            ]["QC"],
            "focal_expected_right": "P/A",
            "focal_right_GT": focal_values[
                "Right"
            ]["GT"],
            "focal_right_biological_GT": (
                focal_values["Right"][
                    "biological_GT"
                ]
            ),
            "focal_right_AD": focal_values[
                "Right"
            ]["AD"],
            "focal_right_DP": focal_values[
                "Right"
            ]["DP"],
            "focal_right_GQ": focal_values[
                "Right"
            ]["GQ"],
            "focal_right_QC": focal_values[
                "Right"
            ]["QC"],
            "exact_focal_zygosity": focal_exact,
            "exact_focal_zygosity_DP2_GQ20": (
                focal_exact_qc
            ),
            "high_confidence_exact_focal": (
                classification_pass
                and focal_exact_qc
            ),
        }

        for hypothesis in HYPOTHESES:
            hypothesis_result = (
                hypothesis_results[hypothesis]
            )

            output[
                f"{hypothesis}_macro_accuracy"
            ] = hypothesis_result[
                "macro_accuracy"
            ]

            output[
                f"{hypothesis}_comparative_supported"
            ] = hypothesis_result[
                "comparative_supported"
            ]

            output[
                f"{hypothesis}_comparative_BA"
            ] = hypothesis_result[
                "comparative_BA"
            ]

            for group in GROUPS:
                metrics = hypothesis_result[
                    "groups"
                ][group]

                output[
                    f"{hypothesis}_{group}_called"
                ] = metrics["called"]

                output[
                    f"{hypothesis}_{group}_correct"
                ] = metrics["correct"]

                output[
                    f"{hypothesis}_{group}_accuracy"
                ] = metrics["accuracy"]

        result_rows.append(output)

    result = pd.DataFrame(
        result_rows
    )

    result = result.sort_values(
        [
            "final_class",
            "best_hypothesis",
            "pos",
        ]
    ).reset_index(drop=True)

    exact = result[
        result[
            "exact_focal_zygosity_DP2_GQ20"
        ]
    ].copy()

    final = result[
        result[
            "high_confidence_exact_focal"
        ]
    ].copy()

    ambiguous = exact[
        ~exact[
            "three_way_classification_pass"
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

    final.to_csv(
        OUT_FINAL,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    ambiguous.to_csv(
        OUT_AMBIGUOUS,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary_rows = []

    for classification in [
        "B",
        "D",
        "R",
    ]:
        subset = result[
            result[
                "best_hypothesis"
            ] == classification
        ]

        summary_rows.append(
            {
                "best_hypothesis": classification,
                "final_class": (
                    "B"
                    if classification == "B"
                    else "DR"
                ),
                "physical_SNPs_best_classified": len(
                    subset
                ),
                "three_way_classification_pass": int(
                    subset[
                        "three_way_classification_pass"
                    ].sum()
                ),
                "exact_focal_zygosity_QC": int(
                    subset[
                        "exact_focal_zygosity_DP2_GQ20"
                    ].sum()
                ),
                "high_confidence_exact_focal": int(
                    subset[
                        "high_confidence_exact_focal"
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
        "THREE-WAY B / DENNIS / RAY CLASSIFICATION"
    )
    print()
    print("EXPLANATORY GROUPS")

    for group in GROUPS:
        print(
            f"{group}: "
            f"{group_totals[group]} samples; "
            f"minimum calls={minimum_group_calls[group]}"
        )

    print()
    print("CLASSIFICATION RULE")
    print(
        f"Macro accuracy >= {MINIMUM_MACRO_ACCURACY:.2f}"
    )
    print(
        f"Margin over second hypothesis >= "
        f"{MINIMUM_SCORE_MARGIN:.2f}"
    )
    print(
        "Winning hypothesis must have independent broad "
        "comparative support."
    )

    print()
    print(
        summary.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print()
    print("HIGH-CONFIDENCE EXACT FOCAL SITES")
    print(
        "class\thypothesis\tposition\tbest_score\t"
        "second\tmargin\tleft\tright"
    )

    for _, row in final.iterrows():
        print(
            f'{row["final_class"]}\t'
            f'{row["best_hypothesis"]}\t'
            f'{row["chrom"]}:{row["pos"]}\t'
            f'{row["best_macro_accuracy"]:.4f}\t'
            f'{row["second_hypothesis"]}:'
            f'{row["second_macro_accuracy"]:.4f}\t'
            f'{row["score_margin"]:.4f}\t'
            f'{row["focal_left_biological_GT"]}\t'
            f'{row["focal_right_biological_GT"]}'
        )

    print()
    print(f"All physical SNPs: {OUT_ALL}")
    print(f"Every exact focal SNP: {OUT_EXACT}")
    print(f"High-confidence final list: {OUT_FINAL}")
    print(f"Ambiguous exact focal sites: {OUT_AMBIGUOUS}")
    print(f"Summary: {OUT_SUMMARY}")

    if s_ba is None:
        print()
        print(
            "NOTE: no comparative BA column was recognized in "
            "the strict selected-site table."
        )
        print(
            "Independent comparative support was retained, but "
            "its numerical BA was not carried into this output."
        )
        print(
            "Selected-site columns were: "
            + ", ".join(
                map(str, selected_raw.columns)
            )
        )
    else:
        print()
        print(
            f"Comparative BA source column: {s_ba}"
        )


if __name__ == "__main__":
    main()
