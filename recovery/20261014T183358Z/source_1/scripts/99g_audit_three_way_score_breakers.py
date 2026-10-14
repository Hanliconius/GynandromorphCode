#!/usr/bin/env python3

from pathlib import Path
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

CLASSIFICATION = (
    CHECKPOINT_DIR
    / "B_DR_three_way_physical_SNP_classification.tsv"
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

OUT_ALL_CALLS = (
    CHECKPOINT_DIR
    / "B_DR_exact_focal_individual_call_audit.tsv"
)

OUT_ALL_BREAKERS = (
    CHECKPOINT_DIR
    / "B_DR_exact_focal_breakers_and_N.tsv"
)

OUT_FINAL_BREAKERS = (
    CHECKPOINT_DIR
    / "B_DR_high_confidence_breakers_and_N.tsv"
)

OUT_SITE_SUMMARY = (
    CHECKPOINT_DIR
    / "B_DR_exact_focal_site_breaker_summary.tsv"
)

OUT_INDIVIDUAL_SUMMARY = (
    CHECKPOINT_DIR
    / "B_DR_individual_breaker_summary.tsv"
)


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

GROUPS = [
    "besckei_B",
    "meriana_D",
    "contigua_R",
]


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
            "Required column was not found.\n"
            f"Tried: {', '.join(candidates)}\n"
            "Available:\n"
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


def biological_dosage(gt, hypothesis):
    gt = canonical_gt(gt)

    if gt is None:
        return np.nan

    alleles = [
        int(allele)
        for allele in gt.split("/")
    ]

    alt_dosage = sum(alleles) / 2.0

    if hypothesis == "B":
        return 1.0 - alt_dosage

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


def explanatory_group(row):
    combination = (
        row["B"],
        row["D"],
        row["R"],
    )

    if combination == (1, 0, 0):
        return "besckei_B"

    if combination == (0, 1, 0):
        return "meriana_D"

    if combination == (0, 0, 1):
        return "contigua_R"

    return ""


def call_status(
    gt,
    dp,
    gq,
    expected_carrier,
    hypothesis,
):
    gt = canonical_gt(gt)

    if gt is None:
        return {
            "status": "N_missing_GT",
            "dosage": np.nan,
            "biological_GT": "missing",
            "observed_carrier": np.nan,
            "callable": False,
            "correct": False,
            "discordant": False,
            "N_or_failed_QC": True,
        }

    failure_reasons = []

    if not np.isfinite(dp):
        failure_reasons.append("missing_DP")
    elif dp < 2:
        failure_reasons.append("low_DP")

    if not np.isfinite(gq):
        failure_reasons.append("missing_GQ")
    elif gq < 20:
        failure_reasons.append("low_GQ")

    dosage = biological_dosage(
        gt,
        hypothesis,
    )

    observed_carrier = dosage > 0

    if failure_reasons:
        return {
            "status": (
                "N_"
                + "_and_".join(failure_reasons)
            ),
            "dosage": dosage,
            "biological_GT": biological_gt(
                dosage
            ),
            "observed_carrier": observed_carrier,
            "callable": False,
            "correct": False,
            "discordant": False,
            "N_or_failed_QC": True,
        }

    correct = (
        observed_carrier
        == expected_carrier
    )

    return {
        "status": (
            "correct"
            if correct
            else "callable_genotype_discordance"
        ),
        "dosage": dosage,
        "biological_GT": biological_gt(
            dosage
        ),
        "observed_carrier": observed_carrier,
        "callable": True,
        "correct": correct,
        "discordant": not correct,
        "N_or_failed_QC": False,
    }


def main():
    for path in [
        CLASSIFICATION,
        GENOTYPES,
        METADATA,
    ]:
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing input: {path}"
            )

    classifications = pd.read_csv(
        CLASSIFICATION,
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

    classifications[
        "exact_focal_zygosity_DP2_GQ20"
    ] = classifications[
        "exact_focal_zygosity_DP2_GQ20"
    ].map(boolean)

    classifications[
        "high_confidence_exact_focal"
    ] = classifications[
        "high_confidence_exact_focal"
    ].map(boolean)

    exact_sites = classifications[
        classifications[
            "exact_focal_zygosity_DP2_GQ20"
        ]
    ].copy()

    site_keys = [
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    exact_sites["pos"] = pd.to_numeric(
        exact_sites["pos"],
        errors="raise",
    ).astype(int)

    # Genotype columns
    g_chrom = find_column(
        genotypes,
        ["chrom", "CHROM", "chromosome"],
    )

    g_pos = find_column(
        genotypes,
        ["pos", "POS", "position"],
    )

    g_ref = find_column(
        genotypes,
        ["ref", "REF"],
    )

    g_alt = find_column(
        genotypes,
        ["alt", "ALT"],
    )

    g_sample = find_column(
        genotypes,
        [
            "sample",
            "Sample",
            "accession",
            "Accession",
        ],
    )

    g_gt = find_column(
        genotypes,
        [
            "GT",
            "GT_raw",
            "canonical_GT",
            "genotype",
        ],
    )

    g_ad = find_column(
        genotypes,
        ["AD", "AD_raw"],
        required=False,
    )

    g_dp = find_column(
        genotypes,
        ["DP", "depth"],
        required=False,
    )

    g_gq = find_column(
        genotypes,
        ["GQ", "genotype_quality"],
        required=False,
    )

    calls = pd.DataFrame(
        {
            "chrom": genotypes[
                g_chrom
            ].str.strip(),
            "pos": pd.to_numeric(
                genotypes[g_pos],
                errors="coerce",
            ),
            "ref": genotypes[
                g_ref
            ].str.strip(),
            "alt": genotypes[
                g_alt
            ].str.strip(),
            "sample": genotypes[
                g_sample
            ].str.strip(),
            "GT": genotypes[g_gt].map(
                canonical_gt
            ),
            "AD": (
                genotypes[g_ad]
                if g_ad is not None
                else ""
            ),
            "DP": (
                pd.to_numeric(
                    genotypes[g_dp],
                    errors="coerce",
                )
                if g_dp is not None
                else np.nan
            ),
            "GQ": (
                pd.to_numeric(
                    genotypes[g_gq],
                    errors="coerce",
                )
                if g_gq is not None
                else np.nan
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
        exact_sites[
            site_keys
            + [
                "best_hypothesis",
                "final_class",
                "best_macro_accuracy",
                "second_hypothesis",
                "second_macro_accuracy",
                "score_margin",
                "best_hypothesis_comparative_BA",
                "three_way_classification_pass",
                "high_confidence_exact_focal",
            ]
        ],
        on=site_keys,
        how="inner",
        validate="many_to_one",
    )

    # Remove module-duplicate copies of identical raw calls.
    calls = (
        calls
        .sort_values(
            site_keys + ["sample"]
        )
        .drop_duplicates(
            site_keys + ["sample"],
            keep="first",
        )
    )

    # Metadata
    m_sample = find_column(
        metadata,
        ["Accession"],
    )

    m_taxon = find_column(
        metadata,
        ["Taxon Name"],
    )

    m_b = find_column(
        metadata,
        ["B"],
    )

    m_d = find_column(
        metadata,
        ["D"],
    )

    m_r = find_column(
        metadata,
        ["R"],
    )

    meta = pd.DataFrame(
        {
            "sample": metadata[
                m_sample
            ].str.strip(),
            "taxon": metadata[
                m_taxon
            ].map(clean_header),
            "B": metadata[m_b].map(
                binary
            ),
            "D": metadata[m_d].map(
                binary
            ),
            "R": metadata[m_r].map(
                binary
            ),
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

    calls["explanatory_group"] = calls.apply(
        explanatory_group,
        axis=1,
    )

    calls = calls[
        calls["explanatory_group"].isin(
            GROUPS
        )
    ].copy()

    diagnostic_rows = []

    for _, row in calls.iterrows():
        hypothesis = row[
            "best_hypothesis"
        ]

        group = row[
            "explanatory_group"
        ]

        expected_carrier = (
            EXPECTED_CARRIER[
                hypothesis
            ][group]
        )

        status = call_status(
            row["GT"],
            row["DP"],
            row["GQ"],
            expected_carrier,
            hypothesis,
        )

        diagnostic_rows.append(
            {
                "chrom": row["chrom"],
                "pos": row["pos"],
                "ref": row["ref"],
                "alt": row["alt"],
                "final_class": row[
                    "final_class"
                ],
                "best_hypothesis": hypothesis,
                "best_macro_accuracy": row[
                    "best_macro_accuracy"
                ],
                "second_hypothesis": row[
                    "second_hypothesis"
                ],
                "second_macro_accuracy": row[
                    "second_macro_accuracy"
                ],
                "score_margin": row[
                    "score_margin"
                ],
                "comparative_BA": row[
                    "best_hypothesis_comparative_BA"
                ],
                "high_confidence_exact_focal": boolean(
                    row[
                        "high_confidence_exact_focal"
                    ]
                ),
                "sample": row["sample"],
                "taxon": row["taxon"],
                "explanatory_group": group,
                "expected_state": (
                    "P-carrier"
                    if expected_carrier
                    else "A/A"
                ),
                "GT": row["GT"] or "missing",
                "biological_GT": status[
                    "biological_GT"
                ],
                "pattern_dosage": status[
                    "dosage"
                ],
                "AD": row["AD"],
                "DP": row["DP"],
                "GQ": row["GQ"],
                "observed_carrier": status[
                    "observed_carrier"
                ],
                "call_status": status[
                    "status"
                ],
                "callable": status[
                    "callable"
                ],
                "correct": status[
                    "correct"
                ],
                "callable_genotype_discordance": status[
                    "discordant"
                ],
                "N_or_failed_QC": status[
                    "N_or_failed_QC"
                ],
            }
        )

    diagnostic = pd.DataFrame(
        diagnostic_rows
    )

    breakers = diagnostic[
        diagnostic["call_status"] != "correct"
    ].copy()

    final_breakers = breakers[
        breakers[
            "high_confidence_exact_focal"
        ]
    ].copy()

    site_summary_rows = []

    for key, subset in diagnostic.groupby(
        site_keys,
        sort=False,
    ):
        first = subset.iloc[0]

        mismatches = subset[
            subset[
                "callable_genotype_discordance"
            ]
        ]

        missing = subset[
            subset[
                "N_or_failed_QC"
            ]
        ]

        output = {
            "chrom": key[0],
            "pos": key[1],
            "ref": key[2],
            "alt": key[3],
            "final_class": first[
                "final_class"
            ],
            "best_hypothesis": first[
                "best_hypothesis"
            ],
            "best_macro_accuracy": first[
                "best_macro_accuracy"
            ],
            "score_margin": first[
                "score_margin"
            ],
            "comparative_BA": first[
                "comparative_BA"
            ],
            "high_confidence_exact_focal": first[
                "high_confidence_exact_focal"
            ],
            "correct_calls": int(
                subset["correct"].sum()
            ),
            "callable_discordances": len(
                mismatches
            ),
            "N_or_failed_QC": len(
                missing
            ),
            "discordant_samples": ",".join(
                mismatches["sample"]
                .astype(str)
                .tolist()
            ),
            "N_or_QC_samples": ",".join(
                missing["sample"]
                .astype(str)
                .tolist()
            ),
        }

        for group in GROUPS:
            group_rows = subset[
                subset[
                    "explanatory_group"
                ] == group
            ]

            output[
                f"{group}_correct"
            ] = int(
                group_rows["correct"].sum()
            )

            output[
                f"{group}_discordant"
            ] = int(
                group_rows[
                    "callable_genotype_discordance"
                ].sum()
            )

            output[
                f"{group}_N_or_QC"
            ] = int(
                group_rows[
                    "N_or_failed_QC"
                ].sum()
            )

        site_summary_rows.append(
            output
        )

    site_summary = pd.DataFrame(
        site_summary_rows
    )

    individual_rows = []

    for sample, subset in diagnostic.groupby(
        "sample",
        sort=False,
    ):
        callable_rows = subset[
            subset["callable"]
        ]

        first = subset.iloc[0]

        discordances = int(
            subset[
                "callable_genotype_discordance"
            ].sum()
        )

        individual_rows.append(
            {
                "sample": sample,
                "taxon": first["taxon"],
                "explanatory_group": first[
                    "explanatory_group"
                ],
                "sites_evaluated": len(
                    subset
                ),
                "correct_calls": int(
                    subset["correct"].sum()
                ),
                "callable_discordances": discordances,
                "N_or_failed_QC": int(
                    subset[
                        "N_or_failed_QC"
                    ].sum()
                ),
                "callable_discordance_fraction": (
                    discordances
                    / len(callable_rows)
                    if len(callable_rows) > 0
                    else np.nan
                ),
            }
        )

    individual_summary = pd.DataFrame(
        individual_rows
    ).sort_values(
        [
            "callable_discordances",
            "N_or_failed_QC",
            "sample",
        ],
        ascending=[
            False,
            False,
            True,
        ],
    )

    diagnostic.to_csv(
        OUT_ALL_CALLS,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    breakers.to_csv(
        OUT_ALL_BREAKERS,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    final_breakers.to_csv(
        OUT_FINAL_BREAKERS,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    site_summary.to_csv(
        OUT_SITE_SUMMARY,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    individual_summary.to_csv(
        OUT_INDIVIDUAL_SUMMARY,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print(
        "THREE-WAY SCORE-BREAKER AUDIT"
    )
    print()
    print(
        "A callable genotype discordance reduces the score."
    )
    print(
        "N/low-DP/low-GQ calls do not reduce the score; "
        "they reduce callability."
    )

    print()
    print(
        "HIGH-CONFIDENCE SITE BREAKERS"
    )
    print(
        "hypothesis\tposition\tscore\t"
        "discordances\tN_or_QC\t"
        "discordant_samples\tN_or_QC_samples"
    )

    final_site_summary = site_summary[
        site_summary[
            "high_confidence_exact_focal"
        ].map(boolean)
    ]

    for _, row in final_site_summary.iterrows():
        print(
            f'{row["best_hypothesis"]}\t'
            f'{row["chrom"]}:{row["pos"]}\t'
            f'{float(row["best_macro_accuracy"]):.4f}\t'
            f'{row["callable_discordances"]}\t'
            f'{row["N_or_failed_QC"]}\t'
            f'{row["discordant_samples"] or "-"}\t'
            f'{row["N_or_QC_samples"] or "-"}'
        )

    print()
    print("INDIVIDUALS WITH CALLABLE DISCORDANCES")
    print(
        "sample\tgroup\tdiscordances\t"
        "callable_discordance_fraction\tN_or_QC\ttaxon"
    )

    recurrent = individual_summary[
        individual_summary[
            "callable_discordances"
        ] > 0
    ]

    for _, row in recurrent.iterrows():
        print(
            f'{row["sample"]}\t'
            f'{row["explanatory_group"]}\t'
            f'{row["callable_discordances"]}\t'
            f'{row["callable_discordance_fraction"]:.4f}\t'
            f'{row["N_or_failed_QC"]}\t'
            f'{row["taxon"]}'
        )

    print()
    print(f"Every individual call: {OUT_ALL_CALLS}")
    print(f"All exact-site breakers/N: {OUT_ALL_BREAKERS}")
    print(f"Final-site breakers/N: {OUT_FINAL_BREAKERS}")
    print(f"Per-site summary: {OUT_SITE_SUMMARY}")
    print(f"Per-individual summary: {OUT_INDIVIDUAL_SUMMARY}")


if __name__ == "__main__":
    main()
