#!/usr/bin/env python3

from pathlib import Path
import importlib.util

import numpy as np
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

HELPER_SCRIPT = (
    PROJECT
    / "scripts"
    / "110b_base_scorer_B_no_DR_positive.py"
)

LONG_VCF_AUDIT = (
    PROJECT
    / "analysis"
    / "checkpoints"
    / "longest_optix_ivory_interval_coverage.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis"
    / "B_heurippa_besckei_no_DR_positive_reanalysis"
    / "checkpoints"
)

OUT_RAW = (
    OUTDIR
    / "longBand_original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.tsv"
)

OUT_PRIMARY = (
    OUTDIR
    / "longBand_original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.DP2_GQ20.tsv"
)

OUT_ACROSS_SCHEMES = (
    OUTDIR
    / "longBand_original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.minimum_across_schemes.tsv"
)

OUT_MATCHING = (
    OUTDIR
    / "longBand_original502_expected_pattern_sites.WGS_vs_WGSplusCS.tsv"
)

OUT_BA80 = (
    OUTDIR
    / "longBand_original502_expected_pattern_BA80_sites.WGS_vs_WGSplusCS.tsv"
)

OUT_SUMMARY = (
    OUTDIR
    / "longBand_original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.summary.tsv"
)


COHORTS = [
    "WGS_only",
    "WGS_plus_CS",
]

PRIMARY_SCHEME = "DP2_GQ20"

# Previously specified a priori Band window.
OLD_BAND_START = 1112651
OLD_BAND_END = 1114314


def fail(message):
    raise RuntimeError(message)


def load_helper():

    if not HELPER_SCRIPT.exists():
        fail(
            f"Missing helper script: {HELPER_SCRIPT}"
        )

    spec = importlib.util.spec_from_file_location(
        "all_snp_helper",
        HELPER_SCRIPT,
    )

    if spec is None or spec.loader is None:
        fail(
            f"Could not import {HELPER_SCRIPT}"
        )

    module = importlib.util.module_from_spec(
        spec
    )

    spec.loader.exec_module(
        module
    )

    return module


def resolve_path(value):

    value = str(
        value
    ).strip()

    if not value:
        fail(
            "Empty source_vcf value in long-region audit."
        )

    p = Path(
        value
    )

    candidates = [
        p,
        PROJECT / p,
        LONG_VCF_AUDIT.parent / p,
    ]

    for candidate in candidates:

        if candidate.exists():
            return candidate.resolve()

    fail(
        "Could not resolve long-region VCF path from "
        f"source_vcf={value!r}"
    )


def choose_long_vcf():

    if not LONG_VCF_AUDIT.exists():
        fail(
            f"Missing long-region audit: {LONG_VCF_AUDIT}"
        )

    d = pd.read_csv(
        LONG_VCF_AUDIT,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "source_vcf",
        "interval_length_bp",
    }

    missing = required - set(
        d.columns
    )

    if missing:
        fail(
            "Long-region audit is missing columns: "
            + ", ".join(
                sorted(
                    missing
                )
            )
        )

    d[
        "interval_length_bp_numeric"
    ] = pd.to_numeric(
        d[
            "interval_length_bp"
        ],
        errors="coerce",
    )

    usable = d.loc[
        d[
            "source_vcf"
        ].astype(str).str.len()
        > 0
    ].copy()

    if usable.empty:
        fail(
            "No source_vcf entries in long-region audit."
        )

    if usable[
        "interval_length_bp_numeric"
    ].notna().any():

        usable = usable.sort_values(
            "interval_length_bp_numeric",
            ascending=False,
        )

    chosen = usable.iloc[
        0
    ]

    path = resolve_path(
        chosen[
            "source_vcf"
        ]
    )

    print(
        "\n=============================================="
    )
    print(
        "LONG BAND SOURCE"
    )
    print(
        "=============================================="
    )

    print(
        "VCF:",
        path,
    )

    for field in [
        "interval",
        "region",
        "interval_length_bp",
        "source_samples",
        "records",
        "biallelic_records",
        "first_record",
        "last_record",
    ]:

        if field in chosen.index:

            print(
                f"{field}:",
                chosen[
                    field
                ],
            )

    if str(
        path
    ).endswith(
        ".gz"
    ):

        fail(
            "The selected long VCF is gzip-compressed, while "
            "the original scorer expects a plain-text VCF. "
            "No output has been modified."
        )

    return path


def location_class(pos):

    pos = int(
        pos
    )

    if pos < OLD_BAND_START:
        return "upstream_of_prior_window"

    if pos > OLD_BAND_END:
        return "downstream_of_prior_window"

    return "within_prior_window"


def make_summary(
    primary,
    across,
):

    rows = []

    for cohort in COHORTS:

        p = primary.loc[
            primary[
                "analysis_cohort"
            ]
            == cohort
        ].copy()

        q = across.loc[
            across[
                "analysis_cohort"
            ]
            == cohort
        ].copy()

        p_pattern = (
            p[
                "score_eligible_bool"
            ]
            &
            p[
                "matches_reference_expectation_bool"
            ]
        )

        q_pattern = (
            q[
                "all_schemes_expected_pattern_eligible"
            ]
            .astype(bool)
        )

        row = {
            "cohort":
                cohort,

            "candidate_SNPs":
                len(
                    p
                ),

            "score_eligible":
                int(
                    p[
                        "score_eligible_bool"
                    ].sum()
                ),

            "expected_pattern_eligible":
                int(
                    p_pattern.sum()
                ),

            "expected_pattern_outside_prior_window":
                int(
                    (
                        p_pattern
                        &
                        (
                            p[
                                "location_class"
                            ]
                            != "within_prior_window"
                        )
                    ).sum()
                ),

            "BA_max_expected_pattern":
                p.loc[
                    p_pattern,
                    "balanced_accuracy",
                ].max(),

            "all_schemes_expected_pattern":
                int(
                    q_pattern.sum()
                ),
        }

        for threshold in [
            0.80,
            0.85,
            0.90,
            0.95,
            0.99,
        ]:

            label = (
                f"{threshold:.2f}"
            )

            passing = (
                p_pattern
                &
                (
                    p[
                        "balanced_accuracy"
                    ]
                    >= threshold
                )
            )

            outside = (
                passing
                &
                (
                    p[
                        "location_class"
                    ]
                    != "within_prior_window"
                )
            )

            row[
                "BA_ge_"
                + label
            ] = int(
                passing.sum()
            )

            row[
                "BA_ge_"
                + label
                + "_outside_prior_window"
            ] = int(
                outside.sum()
            )

            conservative = (
                q_pattern
                &
                (
                    q[
                        "minimum_BA_across_schemes"
                    ]
                    >= threshold
                )
            )

            row[
                "minBA_all_schemes_ge_"
                + label
            ] = int(
                conservative.sum()
            )

        rows.append(
            row
        )

    return pd.DataFrame(
        rows
    )


def main():

    helper = load_helper()

    original = (
        helper.load_original_module()
    )

    metadata = helper.load_metadata(
        original
    )

    manifest = (
        helper.load_final_filtered_manifest()
    )

    retained = manifest.loc[
        manifest[
            "keep_bool"
        ]
        &
        (
            manifest[
                "module"
            ]
            == "Band"
        )
    ].copy()

    print(
        "\nFinal Band cohort sizes from 94b:"
    )

    print(
        retained.groupby(
            "cohort"
        )
        .size()
        .to_string()
    )

    long_vcf = (
        PROJECT
        / "external"
        / "krista_genotype_plots_original"
        / "Genotype_Plots"
        / "Optix_Hmr18-938195-1239421_biallelic.vcf"
    )

    if not long_vcf.exists():
        fail(
            f"Missing long Band VCF: {long_vcf}"
        )

    print(
        "\n=============================================="
    )
    print(
        "LONG BAND SOURCE"
    )
    print(
        "=============================================="
    )
    print(
        "VCF:",
        long_vcf,
    )

    # Start with the exact original Band scoring specification,
    # but substitute the longer VCF.
    specification = dict(
        original.PANELS[
            "band"
        ]
    )

    specification[
        "vcf"
    ] = long_vcf

    print(
        "\nBand scoring specification:"
    )

    print(
        "  trait:",
        specification[
            "trait"
        ],
    )

    print(
        "  expected orientation:",
        specification[
            "expected_orientation"
        ],
    )

    raw_rows = []

    original_make_cohort = (
        original.make_cohort
    )

    try:

        for cohort_name in COHORTS:

            allowed = set(
                retained.loc[
                    retained[
                        "cohort"
                    ]
                    == cohort_name,
                    "sample",
                ]
            )

            print(
                "\n=============================================="
            )

            print(
                cohort_name,
                "/ LONG BAND"
            )

            print(
                "Allowed comparative samples:",
                len(
                    allowed
                ),
            )

            original.make_cohort = (
                helper.make_filtered_cohort_function(
                    allowed
                )
            )

            result = original.score_panel(
                "band_long",
                specification,
                metadata,
            )

            (
                rows,
                score_manifest,
                cohort_summary,
            ) = helper.unpack_score_panel_result(
                result
            )

            print(
                "Scored rows across QC schemes:",
                len(
                    rows
                ),
            )

            for row in rows:

                raw_rows.append(
                    {
                        "analysis_cohort":
                            cohort_name,

                        "module":
                            "Band",

                        "source_vcf":
                            str(
                                long_vcf
                            ),

                        **row,
                    }
                )

    finally:

        original.make_cohort = (
            original_make_cohort
        )

    raw = pd.DataFrame(
        raw_rows
    )

    if raw.empty:
        fail(
            "No SNP scores produced."
        )

    raw = helper.prepare_raw_scores(
        raw
    )

    raw[
        "location_class"
    ] = raw[
        "pos"
    ].map(
        location_class
    )

    raw[
        "distance_from_prior_window_bp"
    ] = np.where(
        raw[
            "pos"
        ]
        < OLD_BAND_START,
        raw[
            "pos"
        ]
        - OLD_BAND_START,
        np.where(
            raw[
                "pos"
            ]
            > OLD_BAND_END,
            raw[
                "pos"
            ]
            - OLD_BAND_END,
            0,
        ),
    )

    raw.to_csv(
        OUT_RAW,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    primary = raw.loc[
        raw[
            "scheme"
        ]
        == PRIMARY_SCHEME
    ].copy()

    if primary.empty:
        fail(
            f"No {PRIMARY_SCHEME} rows produced."
        )

    primary.to_csv(
        OUT_PRIMARY,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    across = (
        helper.make_across_scheme_table(
            raw
        )
    )

    across[
        "location_class"
    ] = across[
        "pos"
    ].map(
        location_class
    )

    across.to_csv(
        OUT_ACROSS_SCHEMES,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    matching = primary.loc[
        primary[
            "score_eligible_bool"
        ]
        &
        primary[
            "matches_reference_expectation_bool"
        ]
    ].copy()

    matching = matching.sort_values(
        [
            "analysis_cohort",
            "balanced_accuracy",
            "pos",
        ],
        ascending=[
            True,
            False,
            True,
        ],
    )

    matching.to_csv(
        OUT_MATCHING,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    ba80 = matching.loc[
        matching[
            "balanced_accuracy"
        ]
        >= 0.80
    ].copy()

    ba80.to_csv(
        OUT_BA80,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary = make_summary(
        primary,
        across,
    )

    summary.to_csv(
        OUT_SUMMARY,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print(
        "\n=============================================="
    )

    print(
        "LONG-BAND SNP SCAN SUMMARY"
    )

    print(
        "=============================================="
    )

    print(
        summary.to_string(
            index=False
        )
    )

    print(
        "\n=============================================="
    )

    print(
        "TOP WGS+CS MATCHING SNPs OUTSIDE PRIOR BAND WINDOW"
    )

    print(
        "=============================================="
    )

    outside = matching.loc[
        (
            matching[
                "analysis_cohort"
            ]
            == "WGS_plus_CS"
        )
        &
        (
            matching[
                "location_class"
            ]
            != "within_prior_window"
        )
    ].copy()

    columns = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "location_class",
        "distance_from_prior_window_bp",
        "comparative_cohort_n",
        "phenotype_0_called",
        "phenotype_1_called",
        "call_rate",
        "balanced_accuracy",
        "sensitivity",
        "specificity",
        "misclassified_n",
        "MCC",
    ]

    if outside.empty:

        print(
            "None."
        )

    else:

        print(
            outside[
                columns
            ]
            .head(
                40
            )
            .to_string(
                index=False
            )
        )

    print(
        "\nOutputs:"
    )

    for path in [
        OUT_RAW,
        OUT_PRIMARY,
        OUT_ACROSS_SCHEMES,
        OUT_MATCHING,
        OUT_BA80,
        OUT_SUMMARY,
    ]:

        print(
            f"  {path}"
        )


if __name__ == "__main__":
    main()
