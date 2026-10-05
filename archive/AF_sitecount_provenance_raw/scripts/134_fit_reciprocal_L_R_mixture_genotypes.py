#!/usr/bin/env python3
from pathlib import Path
import hashlib
import math
import platform
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import minimize, minimize_scalar
from scipy.special import betaln, expit, gammaln, logit, logsumexp


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

INPUT = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "optix_fullVCF_focal_and_WGS_control_AF.long.tsv.gz"
)

EXPECTED_SHA256 = (
    "2129355063e80b19588d605a7402a316"
    "782b38da847f0736a880fa8997a98db9"
)

OUTDIR = (
    PROJECT
    / "analysis/original_502/reciprocal_L_R_mixture"
)

ALL_OUT = (
    OUTDIR
    / "optix_reciprocal_L_R_genotype_pairs.tsv.gz"
)

P95_OUT = (
    OUTDIR
    / "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

CENTRES_OUT = (
    OUTDIR
    / "optix_reciprocal_L_R_expected_AF_centres.tsv"
)

SUMMARY_OUT = (
    OUTDIR
    / "optix_reciprocal_L_R_mixture.summary.tsv"
)

MANIFEST_OUT = (
    OUTDIR
    / "optix_reciprocal_L_R_reproducibility_manifest.tsv"
)

OPTIMIZER_OUT = (
    OUTDIR
    / "optix_reciprocal_L_R_optimizer_multistart.tsv"
)

CHECKSUM_OUT = (
    OUTDIR
    / "optix_reciprocal_L_R_output_checksums.sha256"
)

VERSION = "1.0.0"

LEFT = "TN_Left"
RIGHT = "TN_Right"

MIN_CALL_DP = 2
MIN_FIT_DP = 10
CONTROL_MIN_GQ = 20
POSTERIOR_CUTOFF = 0.95
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

REQUIRED = {
    "chrom",
    "pos",
    "ref",
    "alt",
    "sample",
    "raw_GT",
    "ref_depth",
    "alt_depth",
    "AD_total",
    "DP",
    "GQ",
    "QC_pass",
}


def sha256(path):

    digest = hashlib.sha256()

    with open(path, "rb") as stream:

        for block in iter(
            lambda: stream.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def as_bool(values):

    return (
        values
        .astype(str)
        .str.strip()
        .str.lower()
        .isin(
            {
                "true",
                "t",
                "1",
                "yes",
                "y",
            }
        )
    )


def beta_binomial_logpmf(
    alt_depth,
    total_depth,
    mean_fraction,
    concentration,
):

    alt_depth = np.asarray(
        alt_depth,
        dtype=float,
    )

    total_depth = np.asarray(
        total_depth,
        dtype=float,
    )

    mean_fraction = np.clip(
        np.asarray(
            mean_fraction,
            dtype=float,
        ),
        1e-6,
        1.0 - 1e-6,
    )

    concentration = max(
        float(concentration),
        2.0,
    )

    alpha = (
        mean_fraction
        * concentration
    )

    beta = (
        (1.0 - mean_fraction)
        * concentration
    )

    return (
        gammaln(total_depth + 1.0)
        - gammaln(alt_depth + 1.0)
        - gammaln(
            total_depth
            - alt_depth
            + 1.0
        )
        + betaln(
            alt_depth + alpha,
            total_depth - alt_depth + beta,
        )
        - betaln(
            alpha,
            beta,
        )
    )


def map_fraction(
    biological_fraction,
    error_rate,
    bias_log_odds,
):

    biological_fraction = np.clip(
        np.asarray(
            biological_fraction,
            dtype=float,
        ),
        1e-6,
        1.0 - 1e-6,
    )

    biased_fraction = expit(
        logit(
            biological_fraction
        )
        + bias_log_odds
    )

    return (
        error_rate
        + (
            1.0
            - 2.0 * error_rate
        )
        * biased_fraction
    )


def state_log_likelihoods(
    left_alt,
    left_total,
    right_alt,
    right_total,
    contamination_left,
    contamination_right,
    error_rate,
    bias_log_odds,
    concentration,
):

    biological_left = (
        (
            1.0
            - contamination_left
        )
        * GL
        + contamination_left
        * GR
    )

    biological_right = (
        contamination_right
        * GL
        + (
            1.0
            - contamination_right
        )
        * GR
    )

    expected_left = map_fraction(
        biological_left,
        error_rate,
        bias_log_odds,
    )

    expected_right = map_fraction(
        biological_right,
        error_rate,
        bias_log_odds,
    )

    left_likelihood = (
        beta_binomial_logpmf(
            left_alt[:, None],
            left_total[:, None],
            expected_left[None, :],
            concentration,
        )
    )

    right_likelihood = (
        beta_binomial_logpmf(
            right_alt[:, None],
            right_total[:, None],
            expected_right[None, :],
            concentration,
        )
    )

    return (
        left_likelihood
        + right_likelihood
    )


def main():

    if not INPUT.is_file():

        raise RuntimeError(
            f"Missing input: {INPUT}"
        )

    input_hash = sha256(
        INPUT
    )

    if input_hash != EXPECTED_SHA256:

        raise RuntimeError(
            "Locked input SHA256 mismatch. "
            f"Expected {EXPECTED_SHA256}; "
            f"observed {input_hash}."
        )

    OUTDIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    data = pd.read_csv(
        INPUT,
        sep="\t",
        low_memory=False,
    )

    missing_columns = sorted(
        REQUIRED
        - set(data.columns)
    )

    if missing_columns:

        raise RuntimeError(
            "Missing input columns: "
            + ", ".join(
                missing_columns
            )
        )

    numeric_columns = [
        "pos",
        "ref_depth",
        "alt_depth",
        "AD_total",
        "DP",
        "GQ",
    ]

    for column in numeric_columns:

        data[column] = pd.to_numeric(
            data[column],
            errors="coerce",
        )

    data["ad_n"] = (
        data["ref_depth"]
        + data["alt_depth"]
    )

    # --------------------------------------------------------
    # Calibrate the count model using the three ordinary WGS
    # controls. These samples are not believed to be mixtures.
    # --------------------------------------------------------

    controls = data.loc[
        (
            ~data["sample"].isin(
                {
                    LEFT,
                    RIGHT,
                }
            )
        )
        & as_bool(
            data["QC_pass"]
        )
        & (
            data["ad_n"]
            >= MIN_CALL_DP
        )
        & (
            data["DP"]
            >= MIN_CALL_DP
        )
        & (
            data["GQ"]
            >= CONTROL_MIN_GQ
        )
    ].copy()

    controls["raw_GT"] = (
        controls["raw_GT"]
        .astype(str)
        .str.replace(
            "|",
            "/",
            regex=False,
        )
    )

    homozygous_reference = (
        controls.loc[
            controls["raw_GT"]
            == "0/0"
        ]
    )

    homozygous_alternate = (
        controls.loc[
            controls["raw_GT"]
            == "1/1"
        ]
    )

    heterozygous = (
        controls.loc[
            controls["raw_GT"].isin(
                {
                    "0/1",
                    "1/0",
                }
            )
        ]
    )

    if len(heterozygous) < 50:

        raise RuntimeError(
            "Insufficient control heterozygotes: "
            f"{len(heterozygous)}"
        )

    error_denominator = float(
        homozygous_reference[
            "ad_n"
        ].sum()
        + homozygous_alternate[
            "ad_n"
        ].sum()
    )

    error_rate = float(
        (
            homozygous_reference[
                "alt_depth"
            ].sum()
            + homozygous_alternate[
                "ref_depth"
            ].sum()
        )
        / error_denominator
    )

    error_rate = float(
        np.clip(
            error_rate,
            0.0005,
            0.05,
        )
    )

    heterozygote_mean_AF = float(
        heterozygous[
            "alt_depth"
        ].sum()
        / heterozygous[
            "ad_n"
        ].sum()
    )

    bias_target = (
        (
            heterozygote_mean_AF
            - error_rate
        )
        / (
            1.0
            - 2.0 * error_rate
        )
    )

    bias_target = float(
        np.clip(
            bias_target,
            0.01,
            0.99,
        )
    )

    bias_log_odds = float(
        logit(
            bias_target
        )
    )

    heterozygote_expected_AF = float(
        map_fraction(
            0.5,
            error_rate,
            bias_log_odds,
        )
    )

    control_alt = (
        heterozygous[
            "alt_depth"
        ].to_numpy(
            dtype=float
        )
    )

    control_total = (
        heterozygous[
            "ad_n"
        ].to_numpy(
            dtype=float
        )
    )

    def dispersion_objective(
        log_concentration,
    ):

        concentration = math.exp(
            float(
                log_concentration
            )
        )

        return -float(
            beta_binomial_logpmf(
                control_alt,
                control_total,
                heterozygote_expected_AF,
                concentration,
            ).sum()
        )

    concentration_fit = (
        minimize_scalar(
            dispersion_objective,
            bounds=(
                math.log(2.0),
                math.log(5000.0),
            ),
            method="bounded",
            options={
                "xatol": 1e-5,
            },
        )
    )

    if not concentration_fit.success:

        raise RuntimeError(
            "Control dispersion fit failed: "
            f"{concentration_fit.message}"
        )

    concentration = float(
        math.exp(
            concentration_fit.x
        )
    )

    # --------------------------------------------------------
    # Construct paired focal calls.
    #
    # Focal GQ is intentionally retained for reporting but is
    # not used for filtering. It was calculated under the
    # unmixed diploid model that this analysis is correcting.
    # --------------------------------------------------------

    focal = data.loc[
        data["sample"].isin(
            {
                LEFT,
                RIGHT,
            }
        )
    ].copy()

    keys = [
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    if focal.duplicated(
        keys
        + ["sample"]
    ).any():

        raise RuntimeError(
            "Duplicate focal site/sample records found."
        )

    value_columns = [
        "raw_GT",
        "ref_depth",
        "alt_depth",
        "DP",
        "GQ",
        "ad_n",
    ]

    wide = focal.pivot(
        index=keys,
        columns="sample",
        values=value_columns,
    )

    wide.columns = [
        f"{field}_{sample}"
        for field, sample
        in wide.columns
    ]

    wide = wide.reset_index()

    paired = wide.loc[
        (
            wide[f"ad_n_{LEFT}"]
            >= MIN_CALL_DP
        )
        & (
            wide[f"ad_n_{RIGHT}"]
            >= MIN_CALL_DP
        )
        & (
            wide[f"DP_{LEFT}"]
            >= MIN_CALL_DP
        )
        & (
            wide[f"DP_{RIGHT}"]
            >= MIN_CALL_DP
        )
    ].copy()

    if paired.empty:

        raise RuntimeError(
            "No paired focal calls passed "
            "the depth requirement."
        )

    left_alt = (
        paired[
            f"alt_depth_{LEFT}"
        ].to_numpy(
            dtype=float
        )
    )

    left_total = (
        paired[
            f"ad_n_{LEFT}"
        ].to_numpy(
            dtype=float
        )
    )

    right_alt = (
        paired[
            f"alt_depth_{RIGHT}"
        ].to_numpy(
            dtype=float
        )
    )

    right_total = (
        paired[
            f"ad_n_{RIGHT}"
        ].to_numpy(
            dtype=float
        )
    )

    fitting_sites = (
        (
            left_total
            >= MIN_FIT_DP
        )
        & (
            right_total
            >= MIN_FIT_DP
        )
    )

    if fitting_sites.sum() < 500:

        raise RuntimeError(
            "Insufficient fitting sites: "
            f"{fitting_sites.sum()} passed "
            f"DP >= {MIN_FIT_DP} on both sides."
        )

    # --------------------------------------------------------
    # Estimate reciprocal contamination.
    #
    # At each site the likelihood is marginalized over all nine
    # L/R diploid genotype pairs with equal prior probability.
    # Multiple deterministic starting points test convergence.
    # --------------------------------------------------------

    def mixture_objective(
        parameters,
    ):

        contamination_left = float(
            parameters[0]
        )

        contamination_right = float(
            parameters[1]
        )

        likelihood = (
            state_log_likelihoods(
                left_alt[fitting_sites],
                left_total[fitting_sites],
                right_alt[fitting_sites],
                right_total[fitting_sites],
                contamination_left,
                contamination_right,
                error_rate,
                bias_log_odds,
                concentration,
            )
        )

        return -float(
            logsumexp(
                likelihood
                - math.log(9.0),
                axis=1,
            ).sum()
        )

    starting_points = [
        (0.05, 0.05),
        (0.15, 0.15),
        (0.25, 0.25),
        (0.35, 0.35),
        (0.15, 0.30),
        (0.30, 0.15),
    ]

    fitted_models = []
    optimizer_rows = []

    for (
        starting_left,
        starting_right,
    ) in starting_points:

        fitted = minimize(
            mixture_objective,
            np.array(
                [
                    starting_left,
                    starting_right,
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

        fitted_models.append(
            fitted
        )

        optimizer_rows.append(
            {
                "start_c_left":
                    starting_left,
                "start_c_right":
                    starting_right,
                "success":
                    bool(
                        fitted.success
                    ),
                "fitted_c_left":
                    float(
                        fitted.x[0]
                    ),
                "fitted_c_right":
                    float(
                        fitted.x[1]
                    ),
                "negative_log_likelihood":
                    float(
                        fitted.fun
                    ),
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
                    str(
                        fitted.message
                    ),
            }
        )

    successful_models = [
        fitted
        for fitted
        in fitted_models
        if (
            fitted.success
            and np.isfinite(
                fitted.fun
            )
        )
    ]

    if not successful_models:

        raise RuntimeError(
            "All reciprocal-mixing "
            "optimizer starts failed."
        )

    best_model = min(
        successful_models,
        key=lambda fitted:
            fitted.fun,
    )

    contamination_left = float(
        best_model.x[0]
    )

    contamination_right = float(
        best_model.x[1]
    )

    optimizer_table = pd.DataFrame(
        optimizer_rows
    )

    optimizer_table.to_csv(
        OPTIMIZER_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
    )

    # --------------------------------------------------------
    # Calculate posterior support for all nine genotype pairs.
    # --------------------------------------------------------

    likelihood = state_log_likelihoods(
        left_alt,
        left_total,
        right_alt,
        right_total,
        contamination_left,
        contamination_right,
        error_rate,
        bias_log_odds,
        concentration,
    )

    posterior = np.exp(
        likelihood
        - logsumexp(
            likelihood,
            axis=1,
        )[:, None]
    )

    posterior_order = np.argsort(
        posterior,
        axis=1,
    )

    best_state = (
        posterior_order[:, -1]
    )

    second_state = (
        posterior_order[:, -2]
    )

    row_index = np.arange(
        len(paired)
    )

    biological_left_centres = (
        (
            1.0
            - contamination_left
        )
        * GL
        + contamination_left
        * GR
    )

    biological_right_centres = (
        contamination_right
        * GL
        + (
            1.0
            - contamination_right
        )
        * GR
    )

    expected_left_centres = (
        map_fraction(
            biological_left_centres,
            error_rate,
            bias_log_odds,
        )
    )

    expected_right_centres = (
        map_fraction(
            biological_right_centres,
            error_rate,
            bias_log_odds,
        )
    )

    output = paired[
        keys
    ].copy()

    output["left_original_GT"] = (
        paired[
            f"raw_GT_{LEFT}"
        ].astype(str).to_numpy()
    )

    output["right_original_GT"] = (
        paired[
            f"raw_GT_{RIGHT}"
        ].astype(str).to_numpy()
    )

    output["left_original_GQ"] = (
        paired[
            f"GQ_{LEFT}"
        ].to_numpy()
    )

    output["right_original_GQ"] = (
        paired[
            f"GQ_{RIGHT}"
        ].to_numpy()
    )

    output["left_ref_depth"] = (
        paired[
            f"ref_depth_{LEFT}"
        ].to_numpy()
    )

    output["left_alt_depth"] = (
        left_alt
    )

    output["left_DP"] = (
        left_total
    )

    output["left_observed_AF"] = (
        left_alt
        / left_total
    )

    output["right_ref_depth"] = (
        paired[
            f"ref_depth_{RIGHT}"
        ].to_numpy()
    )

    output["right_alt_depth"] = (
        right_alt
    )

    output["right_DP"] = (
        right_total
    )

    output["right_observed_AF"] = (
        right_alt
        / right_total
    )

    output["mixture_left_GT"] = (
        GTL[
            best_state
        ]
    )

    output["mixture_right_GT"] = (
        GTR[
            best_state
        ]
    )

    output["pair_state"] = (
        output["mixture_left_GT"]
        + "|"
        + output["mixture_right_GT"]
    )

    output["pair_posterior"] = (
        posterior[
            row_index,
            best_state,
        ]
    )

    output["second_pair_state"] = (
        GTL[
            second_state
        ]
        + "|"
        + GTR[
            second_state
        ]
    )

    output["second_pair_posterior"] = (
        posterior[
            row_index,
            second_state,
        ]
    )

    output["posterior_margin"] = (
        output["pair_posterior"]
        - output[
            "second_pair_posterior"
        ]
    )

    output["expected_left_AF"] = (
        expected_left_centres[
            best_state
        ]
    )

    output["expected_right_AF"] = (
        expected_right_centres[
            best_state
        ]
    )

    output[
        "estimated_c_left_from_right"
    ] = contamination_left

    output[
        "estimated_c_right_from_left"
    ] = contamination_right

    output["posterior95"] = (
        output["pair_posterior"]
        >= POSTERIOR_CUTOFF
    )

    output["original_pair_GT"] = (
        output["left_original_GT"]
        + "|"
        + output["right_original_GT"]
    )

    output[
        "changed_from_original_pair"
    ] = (
        output["original_pair_GT"]
        != output["pair_state"]
    )

    output = output.sort_values(
        keys
    )

    output.to_csv(
        ALL_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
        compression={
            "method": "gzip",
            "compresslevel": 6,
            "mtime": 0,
        },
    )

    output.loc[
        output["posterior95"]
    ].to_csv(
        P95_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
    )

    centres = pd.DataFrame(
        {
            "left_true_GT":
                GTL,
            "right_true_GT":
                GTR,
            "left_true_dosage":
                GL,
            "right_true_dosage":
                GR,
            "expected_left_AF":
                expected_left_centres,
            "expected_right_AF":
                expected_right_centres,
            "estimated_c_left_from_right":
                contamination_left,
            "estimated_c_right_from_left":
                contamination_right,
        }
    )

    centres.to_csv(
        CENTRES_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
    )

    summary = pd.DataFrame(
        [
            (
                "input_rows",
                len(data),
            ),
            (
                "input_physical_SNPs",
                data[
                    keys
                ].drop_duplicates().shape[0],
            ),
            (
                "paired_focal_depth_pass",
                len(output),
            ),
            (
                f"fit_sites_DP_ge_{MIN_FIT_DP}",
                int(
                    fitting_sites.sum()
                ),
            ),
            (
                "control_samples",
                controls[
                    "sample"
                ].nunique(),
            ),
            (
                "control_heterozygous_calls",
                len(
                    heterozygous
                ),
            ),
            (
                "estimated_error_rate",
                error_rate,
            ),
            (
                "control_heterozygote_mean_AF",
                heterozygote_mean_AF,
            ),
            (
                "ALT_mapping_bias_log_odds",
                bias_log_odds,
            ),
            (
                "beta_binomial_concentration",
                concentration,
            ),
            (
                "estimated_c_left_from_right",
                contamination_left,
            ),
            (
                "estimated_c_right_from_left",
                contamination_right,
            ),
            (
                "negative_log_likelihood",
                float(
                    best_model.fun
                ),
            ),
            (
                "posterior95_sites",
                int(
                    output[
                        "posterior95"
                    ].sum()
                ),
            ),
            (
                "posterior95_changed_from_original",
                int(
                    (
                        output[
                            "posterior95"
                        ]
                        & output[
                            "changed_from_original_pair"
                        ]
                    ).sum()
                ),
            ),
            (
                "posterior95_left11_right01",
                int(
                    (
                        output[
                            "posterior95"
                        ]
                        & (
                            output[
                                "pair_state"
                            ]
                            == "1/1|0/1"
                        )
                    ).sum()
                ),
            ),
        ],
        columns=[
            "statistic",
            "value",
        ],
    )

    summary.to_csv(
        SUMMARY_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
    )

    manifest = pd.DataFrame(
        [
            (
                "analysis_version",
                VERSION,
            ),
            (
                "script",
                str(
                    Path(
                        __file__
                    ).resolve()
                ),
            ),
            (
                "script_SHA256",
                sha256(
                    Path(
                        __file__
                    ).resolve()
                ),
            ),
            (
                "input",
                str(
                    INPUT
                ),
            ),
            (
                "input_SHA256",
                input_hash,
            ),
            (
                "python_version",
                platform.python_version(),
            ),
            (
                "python_executable",
                sys.executable,
            ),
            (
                "numpy_version",
                np.__version__,
            ),
            (
                "pandas_version",
                pd.__version__,
            ),
            (
                "scipy_version",
                scipy.__version__,
            ),
            (
                "minimum_focal_call_DP",
                MIN_CALL_DP,
            ),
            (
                "minimum_fit_DP_per_side",
                MIN_FIT_DP,
            ),
            (
                "minimum_control_GQ",
                CONTROL_MIN_GQ,
            ),
            (
                "posterior_threshold",
                POSTERIOR_CUTOFF,
            ),
            (
                "genotype_pair_prior",
                "equal_1_over_9",
            ),
            (
                "optimizer",
                "scipy_L-BFGS-B_"
                "deterministic_multistart",
            ),
            (
                "contamination_bounds",
                "0_to_0.45_each_direction",
            ),
            (
                "focal_GQ_used",
                False,
            ),
            (
                "phenotype_expectation_used",
                False,
            ),
            (
                "random_numbers_used",
                False,
            ),
            (
                "float_format",
                FLOAT_FORMAT,
            ),
        ],
        columns=[
            "field",
            "value",
        ],
    )

    manifest.to_csv(
        MANIFEST_OUT,
        sep="\t",
        index=False,
    )

    products = [
        ALL_OUT,
        P95_OUT,
        CENTRES_OUT,
        SUMMARY_OUT,
        MANIFEST_OUT,
        OPTIMIZER_OUT,
    ]

    with open(
        CHECKSUM_OUT,
        "w",
        encoding="utf-8",
        newline="\n",
    ) as stream:

        for path in products:

            stream.write(
                f"{sha256(path)}  "
                f"{path.name}\n"
            )

    state_counts = (
        output.loc[
            output[
                "posterior95"
            ],
            "pair_state",
        ]
        .value_counts()
        .rename_axis(
            "pair_state"
        )
        .reset_index(
            name="sites"
        )
    )

    print(
        "RECIPROCAL L/R "
        "MIXTURE GENOTYPING"
    )

    print()

    print(
        "No phenotype or expected focal "
        "state was used."
    )

    print(
        "Focal GQ was retained for reporting "
        "but was not used for filtering."
    )

    print()

    print(
        summary.to_string(
            index=False
        )
    )

    print(
        "\nPOSTERIOR >= 0.95 "
        "PAIR STATES"
    )

    print(
        state_counts.to_string(
            index=False
        )
    )

    print(
        f"\nAll calls: {ALL_OUT}"
    )

    print(
        f"Posterior-0.95 calls: "
        f"{P95_OUT}"
    )

    print(
        f"State centres: {CENTRES_OUT}"
    )

    print(
        f"Summary: {SUMMARY_OUT}"
    )

    print(
        f"Manifest: {MANIFEST_OUT}"
    )

    print(
        f"Optimizer audit: "
        f"{OPTIMIZER_OUT}"
    )

    print(
        f"Checksums: {CHECKSUM_OUT}"
    )


if __name__ == "__main__":
    main()
