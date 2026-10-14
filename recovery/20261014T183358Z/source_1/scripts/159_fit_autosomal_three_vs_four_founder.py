#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import math
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import minimize
from scipy.special import betaln, expit, gammaln, logit, logsumexp


PROJECT = Path(__file__).resolve().parents[1]

MIXTURE_DIR = (
    PROJECT
    / "analysis/original_502/reciprocal_L_R_mixture"
)

COUNT_INPUT = (
    MIXTURE_DIR
    / "optix_reciprocal_L_R_genotype_pairs.tsv.gz"
)

SUMMARY_INPUT = (
    MIXTURE_DIR
    / "optix_reciprocal_L_R_mixture.summary.tsv"
)

CONTROL_INPUT = (
    PROJECT
    / "data/derived/region_genotypes"
    / "optix_51sample.genotypes.tsv.gz"
)

EXPECTED_COUNT_HASH = (
    "8dd2ebecdcd039ca643e9b68c122407b"
    "fa8f98745df2aac1e00a61901e0aa9a8"
)

EXPECTED_SUMMARY_HASH = (
    "8bfef092d619898612c87fa7ad6702c5"
    "44b79a03ce5f8c558a9bc0acbca158ea"
)

OUTDIR = (
    PROJECT
    / "analysis/autosomal_three_vs_four_founder"
)

MATCHED_OUT = (
    OUTDIR
    / "matched_autosomal_focal_counts.tsv.gz"
)

CONTROL_SUMMARY_OUT = (
    OUTDIR
    / "eligible_control_callability_and_heterozygosity.tsv"
)

SITE_OUT = (
    OUTDIR
    / "observed_site_model_evidence.tsv.gz"
)

BLOCK_OUT = (
    OUTDIR
    / "observed_10kb_block_model_evidence.tsv"
)

FIT_OUT = (
    OUTDIR
    / "observed_model_fits.tsv"
)

OPTIMIZER_OUT = (
    OUTDIR
    / "observed_optimizer_multistart.tsv"
)

AUDIT_OUT = (
    OUTDIR
    / "observed_input_and_model_audit.tsv"
)

CHECKSUM_OUT = (
    OUTDIR
    / "observed_output_checksums.sha256"
)

KEYS = [
    "chrom",
    "pos",
    "ref",
    "alt",
]

FORBIDDEN = {
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
}

SEX_CHROMS = {
    "Z",
    "chrZ",
    "Hmr21",
    "W",
    "chrW",
}

MIN_FIT_DP = 10
BLOCK_SIZE = 10_000
FLOAT_FORMAT = "%.10g"

DOSAGE = np.array(
    [0.0, 0.5, 1.0],
    dtype=float,
)

GT = np.array(
    ["0/0", "0/1", "1/1"],
    dtype=object,
)

GL = np.repeat(DOSAGE, 3)
GR = np.tile(DOSAGE, 3)

GTL = np.repeat(GT, 3)
GTR = np.tile(GT, 3)

STATE_LABELS = np.array(
    [
        f"{left}|{right}"
        for left, right
        in zip(GTL, GTR)
    ],
    dtype=object,
)

M4_STATES = np.arange(
    9,
    dtype=int,
)

M3_STATES = np.array(
    [
        index
        for index, label
        in enumerate(STATE_LABELS)
        if label
        not in {
            "0/0|1/1",
            "1/1|0/0",
        }
    ],
    dtype=int,
)

STARTS = [
    (0.05, 0.05),
    (0.15, 0.15),
    (0.25, 0.25),
    (0.35, 0.35),
    (0.15, 0.30),
    (0.30, 0.15),
]


def arguments():

    parser = argparse.ArgumentParser(
        description=(
            "Fit reciprocal mixing under autosomal "
            "three-founder and four-founder models."
        )
    )

    parser.add_argument(
        "--prepare-only",
        action="store_true",
    )

    return parser.parse_args()


def sha256(path):

    digest = hashlib.sha256()

    with open(path, "rb") as stream:

        for block in iter(
            lambda: stream.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def beta_binomial_logpmf(
    alt,
    total,
    mean,
    concentration,
):

    alt = np.asarray(
        alt,
        dtype=float,
    )

    total = np.asarray(
        total,
        dtype=float,
    )

    mean = np.clip(
        np.asarray(
            mean,
            dtype=float,
        ),
        1e-6,
        1.0 - 1e-6,
    )

    alpha = mean * concentration
    beta = (1.0 - mean) * concentration

    return (
        gammaln(total + 1.0)
        - gammaln(alt + 1.0)
        - gammaln(total - alt + 1.0)
        + betaln(
            alt + alpha,
            total - alt + beta,
        )
        - betaln(alpha, beta)
    )


def map_fraction(
    fraction,
    error_rate,
    bias_log_odds,
):

    fraction = np.clip(
        np.asarray(
            fraction,
            dtype=float,
        ),
        1e-6,
        1.0 - 1e-6,
    )

    biased = expit(
        logit(fraction)
        + bias_log_odds
    )

    return (
        error_rate
        + (
            1.0
            - 2.0 * error_rate
        )
        * biased
    )


def state_log_likelihoods(
    arrays,
    c_left,
    c_right,
    emission,
):

    biological_left = (
        (1.0 - c_left) * GL
        + c_left * GR
    )

    biological_right = (
        c_right * GL
        + (1.0 - c_right) * GR
    )

    expected_left = map_fraction(
        biological_left,
        emission["error_rate"],
        emission["bias_log_odds"],
    )

    expected_right = map_fraction(
        biological_right,
        emission["error_rate"],
        emission["bias_log_odds"],
    )

    return (
        beta_binomial_logpmf(
            arrays["left_alt"][:, None],
            arrays["left_total"][:, None],
            expected_left[None, :],
            emission["concentration"],
        )
        + beta_binomial_logpmf(
            arrays["right_alt"][:, None],
            arrays["right_total"][:, None],
            expected_right[None, :],
            emission["concentration"],
        )
    )


def model_evidence(
    likelihood,
    states,
):

    return logsumexp(
        likelihood[:, states]
        - math.log(len(states)),
        axis=1,
    )


def load_emission():

    summary = pd.read_csv(
        SUMMARY_INPUT,
        sep="\t",
    )

    if set(summary.columns) != {
        "statistic",
        "value",
    }:
        raise RuntimeError(
            "Unexpected mixture-summary columns: "
            f"{list(summary.columns)}"
        )

    values = dict(
        zip(
            summary["statistic"],
            summary["value"],
        )
    )

    required = {
        "estimated_error_rate",
        "control_heterozygote_mean_AF",
        "ALT_mapping_bias_log_odds",
        "beta_binomial_concentration",
    }

    missing = required - set(values)

    if missing:
        raise RuntimeError(
            "Mixture summary missing statistics: "
            f"{sorted(missing)}"
        )

    return {
        "error_rate":
            float(
                values[
                    "estimated_error_rate"
                ]
            ),
        "heterozygote_mean_AF":
            float(
                values[
                    "control_heterozygote_mean_AF"
                ]
            ),
        "bias_log_odds":
            float(
                values[
                    "ALT_mapping_bias_log_odds"
                ]
            ),
        "concentration":
            float(
                values[
                    "beta_binomial_concentration"
                ]
            ),
    }


def prepare_matched():

    counts = pd.read_csv(
        COUNT_INPUT,
        sep="\t",
        low_memory=False,
    )

    required_counts = set(
        KEYS
        + [
            "left_ref_depth",
            "left_alt_depth",
            "left_DP",
            "right_ref_depth",
            "right_alt_depth",
            "right_DP",
        ]
    )

    missing = (
        required_counts
        - set(counts.columns)
    )

    if missing:
        raise RuntimeError(
            "Count table missing columns: "
            f"{sorted(missing)}"
        )

    wanted_control_columns = set(
        KEYS
        + [
            "sample",
            "GT",
            "DP",
            "GQ",
        ]
    )

    controls = pd.read_csv(
        CONTROL_INPUT,
        sep="\t",
        usecols=lambda column:
            column in wanted_control_columns,
        low_memory=False,
    )

    missing_control_columns = (
        wanted_control_columns
        - set(controls.columns)
    )

    if missing_control_columns:
        raise RuntimeError(
            "Control table missing columns: "
            f"{sorted(missing_control_columns)}"
        )

    all_samples = set(
        controls["sample"]
        .dropna()
        .astype(str)
    )

    forbidden_present = sorted(
        all_samples
        & FORBIDDEN
    )

    controls = controls.loc[
        (
            ~controls["sample"]
            .isin(FORBIDDEN)
        )
        & (
            ~controls["chrom"]
            .isin(SEX_CHROMS)
        )
    ].copy()

    eligible_samples = sorted(
        controls["sample"]
        .dropna()
        .astype(str)
        .unique()
    )

    if set(eligible_samples) & FORBIDDEN:
        raise RuntimeError(
            "A forbidden T2/T3 sample survived exclusion."
        )

    if len(eligible_samples) != 47:
        raise RuntimeError(
            "Expected 47 eligible controls after "
            "T2/T3 exclusion; found "
            f"{len(eligible_samples)}."
        )

    site_keys = (
        controls[KEYS]
        .dropna()
        .drop_duplicates()
    )

    sex_rows_removed = int(
        counts["chrom"]
        .isin(SEX_CHROMS)
        .sum()
    )

    counts = counts.loc[
        ~counts["chrom"]
        .isin(SEX_CHROMS)
    ].copy()

    matched = counts.merge(
        site_keys,
        on=KEYS,
        how="inner",
        validate="one_to_one",
    )

    numeric = [
        "pos",
        "left_ref_depth",
        "left_alt_depth",
        "left_DP",
        "right_ref_depth",
        "right_alt_depth",
        "right_DP",
    ]

    for column in numeric:

        matched[column] = pd.to_numeric(
            matched[column],
            errors="coerce",
        )

    matched["left_AD_total"] = (
        matched["left_ref_depth"]
        + matched["left_alt_depth"]
    )

    matched["right_AD_total"] = (
        matched["right_ref_depth"]
        + matched["right_alt_depth"]
    )

    matched = matched.loc[
        (
            matched["left_AD_total"]
            >= MIN_FIT_DP
        )
        & (
            matched["right_AD_total"]
            >= MIN_FIT_DP
        )
        & (
            matched["left_DP"]
            >= MIN_FIT_DP
        )
        & (
            matched["right_DP"]
            >= MIN_FIT_DP
        )
        & (
            matched["left_alt_depth"]
            >= 0
        )
        & (
            matched["left_alt_depth"]
            <= matched["left_AD_total"]
        )
        & (
            matched["right_alt_depth"]
            >= 0
        )
        & (
            matched["right_alt_depth"]
            <= matched["right_AD_total"]
        )
    ].copy()

    matched["pos"] = (
        matched["pos"]
        .astype(int)
    )

    matched = (
        matched
        .sort_values(KEYS)
        .reset_index(drop=True)
    )

    matched["block_start"] = (
        matched["pos"]
        // BLOCK_SIZE
    ) * BLOCK_SIZE

    matched["block_id"] = (
        matched["chrom"].astype(str)
        + ":"
        + matched["block_start"].astype(str)
    )

    if len(matched) < 500:
        raise RuntimeError(
            "Only "
            f"{len(matched)} "
            "matched fitting sites remain."
        )

    matched_keys = (
        matched[KEYS]
        .drop_duplicates()
    )

    control_calls = controls.merge(
        matched_keys,
        on=KEYS,
        how="inner",
        validate="many_to_one",
    )

    control_calls["DP"] = pd.to_numeric(
        control_calls["DP"],
        errors="coerce",
    )

    control_calls["GQ"] = pd.to_numeric(
        control_calls["GQ"],
        errors="coerce",
    )

    gt = (
        control_calls["GT"]
        .astype(str)
        .str.replace(
            "|",
            "/",
            regex=False,
        )
        .str.replace(
            "1/0",
            "0/1",
            regex=False,
        )
    )

    control_calls["canonical_GT"] = gt.where(
        gt.isin(
            {
                "0/0",
                "0/1",
                "1/1",
            }
        )
    )

    control_calls["callable"] = (
        control_calls[
            "canonical_GT"
        ].notna()
        & (
            control_calls["DP"]
            >= 2
        )
        & (
            control_calls["GQ"]
            >= 20
        )
    )

    control_calls["heterozygous"] = (
        control_calls["callable"]
        & control_calls[
            "canonical_GT"
        ].eq("0/1")
    )

    control_summary = (
        control_calls
        .groupby(
            "sample",
            as_index=False,
        )
        .agg(
            matched_sites=(
                "pos",
                "size",
            ),
            callable_sites=(
                "callable",
                "sum",
            ),
            heterozygous_sites=(
                "heterozygous",
                "sum",
            ),
        )
    )

    control_summary["call_rate"] = (
        control_summary["callable_sites"]
        / control_summary["matched_sites"]
    )

    control_summary[
        "heterozygosity_among_callable"
    ] = (
        control_summary[
            "heterozygous_sites"
        ]
        / control_summary[
            "callable_sites"
        ]
    )

    control_summary = (
        control_summary
        .sort_values(
            [
                "heterozygosity_among_callable",
                "sample",
            ]
        )
        .reset_index(drop=True)
    )

    details = {
        "eligible_control_samples":
            len(eligible_samples),
        "forbidden_present_before_exclusion":
            ",".join(forbidden_present),
        "forbidden_after_exclusion":
            0,
        "control_sites":
            len(site_keys),
        "sex_chromosome_rows_removed":
            sex_rows_removed,
        "matched_fitting_sites":
            len(matched),
        "matched_chromosomes":
            ",".join(
                sorted(
                    matched["chrom"]
                    .astype(str)
                    .unique()
                )
            ),
        "ten_kb_blocks":
            matched["block_id"].nunique(),
        "minimum_control_call_rate":
            control_summary[
                "call_rate"
            ].min(),
        "maximum_control_call_rate":
            control_summary[
                "call_rate"
            ].max(),
        "minimum_control_heterozygosity":
            control_summary[
                "heterozygosity_among_callable"
            ].min(),
        "maximum_control_heterozygosity":
            control_summary[
                "heterozygosity_among_callable"
            ].max(),
    }

    return (
        matched,
        control_summary,
        details,
    )


def fit_model(
    name,
    states,
    arrays,
    emission,
):

    def objective(parameters):

        likelihood = state_log_likelihoods(
            arrays,
            float(parameters[0]),
            float(parameters[1]),
            emission,
        )

        return -float(
            model_evidence(
                likelihood,
                states,
            ).sum()
        )

    results = []
    optimizer_rows = []

    for start_left, start_right in STARTS:

        fitted = minimize(
            objective,
            np.array(
                [
                    start_left,
                    start_right,
                ],
                dtype=float,
            ),
            method="L-BFGS-B",
            bounds=(
                (0.0, 0.45),
                (0.0, 0.45),
            ),
            options={
                "maxiter": 250,
                "ftol": 1e-10,
            },
        )

        results.append(fitted)

        optimizer_rows.append(
            {
                "model":
                    name,
                "start_c_left":
                    start_left,
                "start_c_right":
                    start_right,
                "success":
                    bool(fitted.success),
                "fitted_c_left":
                    float(fitted.x[0]),
                "fitted_c_right":
                    float(fitted.x[1]),
                "negative_log_evidence":
                    float(fitted.fun),
                "iterations":
                    int(
                        getattr(
                            fitted,
                            "nit",
                            -1,
                        )
                    ),
                "function_evaluations":
                    int(
                        getattr(
                            fitted,
                            "nfev",
                            -1,
                        )
                    ),
                "message":
                    str(fitted.message),
            }
        )

    successful = [
        fitted
        for fitted in results
        if (
            fitted.success
            and np.isfinite(
                fitted.fun
            )
        )
    ]

    if not successful:
        raise RuntimeError(
            "Every optimizer start failed "
            f"for {name}."
        )

    best = min(
        successful,
        key=lambda fitted:
            fitted.fun,
    )

    c_left = float(
        best.x[0]
    )

    c_right = float(
        best.x[1]
    )

    likelihood = state_log_likelihoods(
        arrays,
        c_left,
        c_right,
        emission,
    )

    evidence = model_evidence(
        likelihood,
        states,
    )

    restricted = (
        likelihood[:, states]
        - math.log(len(states))
    )

    posterior = np.exp(
        restricted
        - logsumexp(
            restricted,
            axis=1,
        )[:, None]
    )

    local_best = np.argmax(
        posterior,
        axis=1,
    )

    full_posterior = np.zeros(
        (
            len(likelihood),
            9,
        ),
        dtype=float,
    )

    full_posterior[
        :,
        states,
    ] = posterior

    fit = {
        "model":
            name,
        "states":
            len(states),
        "state_prior":
            f"equal_1_over_{len(states)}",
        "c_left_from_right":
            c_left,
        "c_right_from_left":
            c_right,
        "log_marginal_evidence":
            float(evidence.sum()),
        "mean_log_marginal_evidence_per_site":
            float(evidence.mean()),
    }

    return {
        "fit":
            fit,
        "optimizer":
            optimizer_rows,
        "evidence":
            evidence,
        "best_state":
            STATE_LABELS[
                states[local_best]
            ],
        "best_posterior":
            posterior[
                np.arange(
                    len(posterior)
                ),
                local_best,
            ],
        "posterior":
            full_posterior,
    }


def write_checksums(paths):

    lines = [
        (
            f"{sha256(path)}  "
            f"{path.relative_to(PROJECT)}"
        )
        for path in paths
    ]

    CHECKSUM_OUT.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def main():

    args = arguments()

    for path in (
        COUNT_INPUT,
        SUMMARY_INPUT,
        CONTROL_INPUT,
    ):
        if not path.is_file():
            raise RuntimeError(
                f"Missing input: {path}"
            )

    count_hash = sha256(
        COUNT_INPUT
    )

    summary_hash = sha256(
        SUMMARY_INPUT
    )

    if count_hash != EXPECTED_COUNT_HASH:
        raise RuntimeError(
            "Locked raw-count table hash mismatch: "
            f"{count_hash}"
        )

    if summary_hash != EXPECTED_SUMMARY_HASH:
        raise RuntimeError(
            "Locked emission-summary hash mismatch: "
            f"{summary_hash}"
        )

    OUTDIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    emission = load_emission()

    (
        matched,
        control_summary,
        details,
    ) = prepare_matched()

    matched.to_csv(
        MATCHED_OUT,
        sep="\t",
        index=False,
        compression="gzip",
        float_format=FLOAT_FORMAT,
    )

    control_summary.to_csv(
        CONTROL_SUMMARY_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
    )

    audit = [
        (
            "count_input",
            COUNT_INPUT,
        ),
        (
            "count_input_sha256",
            count_hash,
        ),
        (
            "emission_summary_input",
            SUMMARY_INPUT,
        ),
        (
            "emission_summary_sha256",
            summary_hash,
        ),
        (
            "control_input",
            CONTROL_INPUT,
        ),
        (
            "control_input_sha256",
            sha256(CONTROL_INPUT),
        ),
        (
            "focal_likelihood_input",
            "raw merged Left/Right ref and alt counts",
        ),
        (
            "original_focal_GT_used",
            False,
        ),
        (
            "original_focal_GQ_used",
            False,
        ),
        (
            "sex_chromosomes_included",
            False,
        ),
        (
            "two_founder_model_included",
            False,
        ),
        (
            "three_founder_disallowed_states",
            "0/0|1/1,1/1|0/0",
        ),
        (
            "four_founder_states",
            "all_nine_diploid_pairs",
        ),
        (
            "block_size_bp",
            BLOCK_SIZE,
        ),
    ]

    audit += list(
        details.items()
    )

    audit += list(
        emission.items()
    )

    audit_table = pd.DataFrame(
        audit,
        columns=[
            "item",
            "value",
        ],
    )

    audit_table.to_csv(
        AUDIT_OUT,
        sep="\t",
        index=False,
    )

    print(
        "AUTOSOMAL THREE- VERSUS "
        "FOUR-FOUNDER INPUT PREPARATION"
    )

    print(
        audit_table.to_string(
            index=False
        )
    )

    print(
        f"\nMatched table: {MATCHED_OUT}"
    )

    print(
        "Control summary: "
        f"{CONTROL_SUMMARY_OUT}"
    )

    print(
        f"Audit: {AUDIT_OUT}"
    )

    if args.prepare_only:

        write_checksums(
            [
                MATCHED_OUT,
                CONTROL_SUMMARY_OUT,
                AUDIT_OUT,
            ]
        )

        print(
            "\nPrepare-only validation completed; "
            "model fitting was not run."
        )

        return

    arrays = {
        "left_alt":
            matched[
                "left_alt_depth"
            ].to_numpy(
                dtype=float
            ),
        "left_total":
            matched[
                "left_AD_total"
            ].to_numpy(
                dtype=float
            ),
        "right_alt":
            matched[
                "right_alt_depth"
            ].to_numpy(
                dtype=float
            ),
        "right_total":
            matched[
                "right_AD_total"
            ].to_numpy(
                dtype=float
            ),
    }

    m3 = fit_model(
        "three_founder",
        M3_STATES,
        arrays,
        emission,
    )

    m4 = fit_model(
        "four_founder",
        M4_STATES,
        arrays,
        emission,
    )

    fits = pd.DataFrame(
        [
            m3["fit"],
            m4["fit"],
        ]
    )

    fits.to_csv(
        FIT_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
    )

    optimizer = pd.DataFrame(
        m3["optimizer"]
        + m4["optimizer"]
    )

    optimizer.to_csv(
        OPTIMIZER_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
    )

    site = matched[
        KEYS
        + [
            "block_start",
            "block_id",
        ]
    ].copy()

    site[
        "three_founder_log_evidence"
    ] = m3["evidence"]

    site[
        "four_founder_log_evidence"
    ] = m4["evidence"]

    site[
        "delta_log_evidence_four_minus_three"
    ] = (
        m4["evidence"]
        - m3["evidence"]
    )

    site[
        "three_founder_best_pair"
    ] = m3["best_state"]

    site[
        "three_founder_best_pair_posterior"
    ] = m3["best_posterior"]

    site[
        "four_founder_best_pair"
    ] = m4["best_state"]

    site[
        "four_founder_best_pair_posterior"
    ] = m4["best_posterior"]

    site[
        "four_founder_opposing_homozygote_posterior"
    ] = (
        m4["posterior"][:, 2]
        + m4["posterior"][:, 6]
    )

    site.to_csv(
        SITE_OUT,
        sep="\t",
        index=False,
        compression="gzip",
        float_format=FLOAT_FORMAT,
    )

    block = (
        site
        .groupby(
            [
                "chrom",
                "block_start",
                "block_id",
            ],
            as_index=False,
        )
        .agg(
            sites=(
                "pos",
                "size",
            ),
            delta_log_evidence_four_minus_three=(
                "delta_log_evidence_four_minus_three",
                "sum",
            ),
            mean_delta_per_site=(
                "delta_log_evidence_four_minus_three",
                "mean",
            ),
            opposing_homozygote_posterior_sum=(
                "four_founder_opposing_homozygote_posterior",
                "sum",
            ),
        )
    )

    block.to_csv(
        BLOCK_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
    )

    delta = float(
        site[
            "delta_log_evidence_four_minus_three"
        ].sum()
    )

    if delta > 0:
        favored = "four_founder"
    else:
        favored = "three_founder"

    print(
        "\nOBSERVED MODEL FITS"
    )

    print(
        fits.to_string(
            index=False
        )
    )

    print(
        "\nOBSERVED RAW MODEL COMPARISON"
    )

    print(
        "delta_log_evidence_four_minus_three="
        f"{delta:.10g}"
    )

    print(
        "mean_delta_per_site="
        f"{delta / len(site):.10g}"
    )

    print(
        f"raw_favored_model={favored}"
    )

    print(
        "This is not yet a calibrated founder-model "
        "probability; simulation comes next."
    )

    print(
        f"\nFits: {FIT_OUT}"
    )

    print(
        f"Site evidence: {SITE_OUT}"
    )

    print(
        f"Block evidence: {BLOCK_OUT}"
    )

    write_checksums(
        [
            MATCHED_OUT,
            CONTROL_SUMMARY_OUT,
            AUDIT_OUT,
            FIT_OUT,
            OPTIMIZER_OUT,
            SITE_OUT,
            BLOCK_OUT,
        ]
    )

    print(
        f"Checksums: {CHECKSUM_OUT}"
    )

    print(
        f"Runtime: Python {platform.python_version()}, "
        f"NumPy {np.__version__}, "
        f"pandas {pd.__version__}, "
        f"SciPy {scipy.__version__}"
    )


if __name__ == "__main__":
    main()
