#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import itertools
import math
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import scipy
from scipy.optimize import minimize
from scipy.special import expit, logit, logsumexp


PROJECT = Path(__file__).resolve().parents[1]
MODEL_SCRIPT = PROJECT / "scripts/159_fit_autosomal_three_vs_four_founder.py"
OUTDIR = PROJECT / "analysis/autosomal_three_vs_four_founder"
MATCHED_INPUT = OUTDIR / "matched_autosomal_focal_counts.tsv.gz"
FITS_INPUT = OUTDIR / "observed_model_fits.tsv"
CONTROL_INPUT = PROJECT / "data/derived/region_genotypes/optix_51sample.genotypes.tsv.gz"
AUDIT_INPUT = OUTDIR / "observed_input_and_model_audit.tsv"

OBSERVED_OUT = OUTDIR / "observed_nested_founder_fit.tsv"
SIM_OUT = OUTDIR / "founder_simulation_replicates.tsv.gz"
SIM_SUMMARY_OUT = OUTDIR / "founder_simulation_summary.tsv"
MIXING_OUT = OUTDIR / "mixing_calibration.tsv"
BOOTSTRAP_OUT = OUTDIR / "observed_10kb_block_bootstrap.tsv.gz"
FINAL_OUT = OUTDIR / "calibrated_founder_and_mixing_summary.tsv"
DONOR_OUT = OUTDIR / "simulation_all_47_controls.tsv"
FIGURE_OUT = OUTDIR / "three_vs_four_founder_simulation_calibration.pdf"
FIGURE_PNG = OUTDIR / "three_vs_four_founder_simulation_calibration.summary.png"
CHECKSUM_OUT = OUTDIR / "simulation_output_checksums.sha256"

KEYS = ["chrom", "pos", "ref", "alt"]
FORBIDDEN = {"TN2L_S1", "TN2R_S2", "TN3L_S3", "TN3R_S4"}
SEX_CHROMS = {"Z", "chrZ", "Hmr21", "W", "chrW"}
OPPOSING = np.array([2, 6], dtype=int)
COMPATIBLE = np.array([0, 1, 3, 4, 5, 7, 8], dtype=int)
GRID = np.array([0.00, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40])
HIGH_CALL_RATE = 0.65
NEAREST_MAIN = 60
NEAREST_HIGH = 40
FLOAT_FORMAT = "%.10g"
SEED = 731_947


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates-per-model", type=int, default=196)
    parser.add_argument("--bootstrap-replicates", type=int, default=200)
    parser.add_argument("--prepare-only", action="store_true")
    return parser.parse_args()


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_model_module():
    spec = importlib.util.spec_from_file_location("founder_model_159", MODEL_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load model script: {MODEL_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_control_matrix(matched):
    wanted = set(KEYS + ["sample", "GT", "DP", "GQ"])
    calls = pd.read_csv(
        CONTROL_INPUT, sep="\t", usecols=lambda x: x in wanted, low_memory=False
    )
    missing = wanted - set(calls.columns)
    if missing:
        raise RuntimeError(f"Control input missing columns: {sorted(missing)}")

    calls = calls.loc[
        ~calls["sample"].isin(FORBIDDEN)
        & ~calls["chrom"].isin(SEX_CHROMS)
    ].copy()

    if set(calls["sample"].astype(str)) & FORBIDDEN:
        raise RuntimeError("A forbidden T2/T3 control survived exclusion.")

    site_order = matched[KEYS].copy()
    site_order["site_index"] = np.arange(len(site_order))

    calls = calls.merge(
        site_order,
        on=KEYS,
        how="inner",
        validate="many_to_one",
    )

    calls["DP"] = pd.to_numeric(calls["DP"], errors="coerce")
    calls["GQ"] = pd.to_numeric(calls["GQ"], errors="coerce")

    gt = (
        calls["GT"]
        .astype(str)
        .str.replace("|", "/", regex=False)
        .str.replace("1/0", "0/1", regex=False)
    )

    calls["dosage"] = gt.map({
        "0/0": 0.0,
        "0/1": 1.0,
        "1/1": 2.0,
    })

    calls.loc[
        (calls["DP"] < 2) | (calls["GQ"] < 20),
        "dosage",
    ] = np.nan

    matrix = calls.pivot(
        index="site_index",
        columns="sample",
        values="dosage",
    )
    matrix = matrix.reindex(np.arange(len(matched)))

    if matrix.shape[1] != 47:
        raise RuntimeError(
            f"Expected 47 control columns; found {matrix.shape[1]}."
        )

    donor = pd.DataFrame({
        "sample": matrix.columns,
        "callable_sites": matrix.notna().sum(axis=0).to_numpy(),
        "heterozygous_sites": matrix.eq(1.0).sum(axis=0).to_numpy(),
    })

    donor["total_sites"] = len(matrix)
    donor["call_rate"] = (
        donor["callable_sites"]
        / donor["total_sites"]
    )
    donor["heterozygosity"] = (
        donor["heterozygous_sites"]
        / donor["callable_sites"]
    )
    donor["high_call_sensitivity_set"] = (
        donor["call_rate"] >= HIGH_CALL_RATE
    )

    included = donor.sort_values(
        ["heterozygosity", "sample"]
    ).reset_index(drop=True)

    if len(included) != 47:
        raise RuntimeError(
            f"Expected all 47 controls; found {len(included)}."
        )

    included["heterozygosity_stratum"] = np.minimum(
        2,
        np.floor(
            np.arange(len(included)) * 3 / len(included)
        ).astype(int),
    )

    donor = donor.merge(
        included[["sample", "heterozygosity_stratum"]],
        on="sample",
        how="left",
    )

    dosage = matrix.to_numpy(dtype=float)
    called = np.isfinite(dosage)
    alternate = np.nansum(dosage, axis=1)
    chromosomes = 2.0 * called.sum(axis=1)

    allele_frequency = (
        alternate + 0.5
    ) / (
        chromosomes + 1.0
    )
    allele_frequency = np.clip(
        allele_frequency,
        1e-4,
        1 - 1e-4,
    )

    return (
        matrix,
        donor.sort_values("sample"),
        included,
        allele_frequency,
    )


def marginal_evidence(likelihood, epsilon=None):
    if epsilon is None:
        return logsumexp(
            likelihood[:, COMPATIBLE]
            - math.log(len(COMPATIBLE)),
            axis=1,
        )

    epsilon = float(
        np.clip(epsilon, 1e-8, 0.5)
    )

    log_prior = np.empty(9, dtype=float)

    log_prior[COMPATIBLE] = (
        math.log1p(-epsilon)
        - math.log(len(COMPATIBLE))
    )
    log_prior[OPPOSING] = (
        math.log(epsilon)
        - math.log(len(OPPOSING))
    )

    return logsumexp(
        likelihood + log_prior[None, :],
        axis=1,
    )


def fit_null(model, arrays, emission, starts, weights=None):
    if weights is None:
        weights = np.ones(
            len(arrays["left_alt"]),
            dtype=float,
        )

    def objective(parameters):
        likelihood = model.state_log_likelihoods(
            arrays,
            float(parameters[0]),
            float(parameters[1]),
            emission,
        )
        return -float(
            np.dot(
                weights,
                marginal_evidence(likelihood),
            )
        )

    candidates = []
    successful = False

    for start_values in starts:
        initial = np.asarray(
            start_values,
            dtype=float,
        )
        initial_value = objective(initial)

        if np.isfinite(initial_value):
            candidates.append(
                (
                    float(initial_value),
                    initial.copy(),
                )
            )

        result = minimize(
            objective,
            initial,
            method="L-BFGS-B",
            bounds=(
                (0.0, 0.45),
                (0.0, 0.45),
            ),
            options={
                "maxiter": 200,
                "ftol": 1e-9,
            },
        )

        if (
            np.isfinite(result.fun)
            and np.all(np.isfinite(result.x))
        ):
            candidates.append(
                (
                    float(result.fun),
                    np.asarray(
                        result.x,
                        dtype=float,
                    ),
                )
            )
            successful = (
                successful
                or bool(result.success)
            )

    if not candidates:
        raise RuntimeError(
            "Every three-founder objective "
            "evaluation was non-finite."
        )

    if not successful:
        fallback_start = min(
            candidates,
            key=lambda x: x[0],
        )[1]

        fallback = minimize(
            objective,
            fallback_start,
            method="Powell",
            bounds=(
                (0.0, 0.45),
                (0.0, 0.45),
            ),
            options={
                "maxiter": 500,
                "xtol": 1e-6,
                "ftol": 1e-8,
            },
        )

        if (
            np.isfinite(fallback.fun)
            and np.all(np.isfinite(fallback.x))
        ):
            candidates.append(
                (
                    float(fallback.fun),
                    np.asarray(
                        fallback.x,
                        dtype=float,
                    ),
                )
            )

        print(
            "Three-founder L-BFGS-B fallback used.",
            flush=True,
        )

    best_fun, best_x = min(
        candidates,
        key=lambda x: x[0],
    )

    return {
        "c_left": float(best_x[0]),
        "c_right": float(best_x[1]),
        "epsilon": 0.0,
        "log_evidence": -float(best_fun),
    }


def fit_alternative(
    model,
    arrays,
    emission,
    starts,
    weights=None,
):
    if weights is None:
        weights = np.ones(
            len(arrays["left_alt"]),
            dtype=float,
        )

    def objective(parameters):
        epsilon = float(
            expit(parameters[2]) * 0.5
        )

        likelihood = model.state_log_likelihoods(
            arrays,
            float(parameters[0]),
            float(parameters[1]),
            emission,
        )

        return -float(
            np.dot(
                weights,
                marginal_evidence(
                    likelihood,
                    epsilon,
                ),
            )
        )

    candidates = []
    successful = False

    for c_left, c_right, epsilon in starts:
        z = float(
            logit(
                np.clip(
                    epsilon / 0.5,
                    1e-7,
                    1 - 1e-7,
                )
            )
        )

        initial = np.array(
            [c_left, c_right, z],
            dtype=float,
        )
        initial_value = objective(initial)

        if np.isfinite(initial_value):
            candidates.append(
                (
                    float(initial_value),
                    initial.copy(),
                )
            )

        result = minimize(
            objective,
            initial,
            method="L-BFGS-B",
            bounds=(
                (0.0, 0.45),
                (0.0, 0.45),
                (-16.0, 16.0),
            ),
            options={
                "maxiter": 250,
                "ftol": 1e-9,
            },
        )

        if (
            np.isfinite(result.fun)
            and np.all(np.isfinite(result.x))
        ):
            candidates.append(
                (
                    float(result.fun),
                    np.asarray(
                        result.x,
                        dtype=float,
                    ),
                )
            )
            successful = (
                successful
                or bool(result.success)
            )

    if not candidates:
        raise RuntimeError(
            "Every four-founder objective "
            "evaluation was non-finite."
        )

    if not successful:
        fallback_start = min(
            candidates,
            key=lambda x: x[0],
        )[1]

        fallback = minimize(
            objective,
            fallback_start,
            method="Powell",
            bounds=(
                (0.0, 0.45),
                (0.0, 0.45),
                (-16.0, 16.0),
            ),
            options={
                "maxiter": 750,
                "xtol": 1e-6,
                "ftol": 1e-8,
            },
        )

        if (
            np.isfinite(fallback.fun)
            and np.all(np.isfinite(fallback.x))
        ):
            candidates.append(
                (
                    float(fallback.fun),
                    np.asarray(
                        fallback.x,
                        dtype=float,
                    ),
                )
            )

        print(
            "Four-founder alternative "
            "L-BFGS-B fallback used.",
            flush=True,
        )

    best_fun, best_x = min(
        candidates,
        key=lambda x: x[0],
    )

    return {
        "c_left": float(best_x[0]),
        "c_right": float(best_x[1]),
        "epsilon": float(
            expit(best_x[2]) * 0.5
        ),
        "log_evidence": -float(best_fun),
    }


def fitted_pair(
    model,
    arrays,
    emission,
    null_starts,
    alternative_starts,
    weights=None,
):
    null = fit_null(
        model,
        arrays,
        emission,
        null_starts,
        weights,
    )

    alternative_starts = list(
        alternative_starts
    ) + [
        (
            null["c_left"],
            null["c_right"],
            1e-7,
        )
    ]

    alternative = fit_alternative(
        model,
        arrays,
        emission,
        alternative_starts,
        weights,
    )

    statistic = max(
        0.0,
        2.0 * (
            alternative["log_evidence"]
            - null["log_evidence"]
        ),
    )

    return null, alternative, statistic


def imputed_dosage(
    matrix,
    sample,
    allele_frequency,
    rng,
):
    values = matrix[sample].to_numpy(
        dtype=float
    ).copy()

    missing = ~np.isfinite(values)

    values[missing] = rng.binomial(
        2,
        allele_frequency[missing],
    )

    return values.astype(int)


def gamete(dosage, rng):
    result = np.zeros(
        len(dosage),
        dtype=int,
    )

    result[dosage == 2] = 1
    heterozygous = dosage == 1

    result[heterozygous] = rng.integers(
        0,
        2,
        heterozygous.sum(),
    )

    return result


def choose_donors(
    included,
    replicate,
    rng,
):
    samples = included[
        "sample"
    ].astype(str).tolist()

    first_index = (
        replicate - 1
    ) % len(samples)

    first = samples[first_index]

    offset = 1 + (
        (replicate * 17)
        % (len(samples) - 1)
    )

    second = samples[
        (first_index + offset)
        % len(samples)
    ]

    if second == first:
        raise RuntimeError(
            "Donor pairing produced the same control twice."
        )

    return first, second


def true_genotypes(
    founder_model,
    matrix,
    included,
    allele_frequency,
    replicate,
    rng,
):
    donor_a, donor_b = choose_donors(
        included,
        replicate,
        rng,
    )

    dosage_a = imputed_dosage(
        matrix,
        donor_a,
        allele_frequency,
        rng,
    )
    dosage_b = imputed_dosage(
        matrix,
        donor_b,
        allele_frequency,
        rng,
    )

    if founder_model == "four_founder":
        left = dosage_a
        right = dosage_b
        construction = (
            "two_real_nonmixed_diploid_controls"
        )
    else:
        shared = gamete(
            dosage_a,
            rng,
        )
        private_b = gamete(
            dosage_b,
            rng,
        )
        constructed = shared + private_b

        if replicate % 2 == 0:
            left = dosage_a
            right = constructed
            construction = (
                "real_left_plus_constructed_"
                "shared_haplotype_right"
            )
        else:
            left = constructed
            right = dosage_a
            construction = (
                "constructed_shared_haplotype_"
                "left_plus_real_right"
            )

    return (
        left,
        right,
        donor_a,
        donor_b,
        construction,
    )


def simulate_arrays(
    model,
    true_left,
    true_right,
    c_left,
    c_right,
    depths,
    emission,
    rng,
):
    biological_left = (
        (1 - c_left) * (true_left / 2)
        + c_left * (true_right / 2)
    )
    biological_right = (
        c_right * (true_left / 2)
        + (1 - c_right) * (true_right / 2)
    )

    mean_left = model.map_fraction(
        biological_left,
        emission["error_rate"],
        emission["bias_log_odds"],
    )
    mean_right = model.map_fraction(
        biological_right,
        emission["error_rate"],
        emission["bias_log_odds"],
    )

    concentration = emission["concentration"]

    latent_left = rng.beta(
        mean_left * concentration,
        (1 - mean_left) * concentration,
    )
    latent_right = rng.beta(
        mean_right * concentration,
        (1 - mean_right) * concentration,
    )

    left_alt = rng.binomial(
        depths["left_total"].astype(int),
        latent_left,
    )
    right_alt = rng.binomial(
        depths["right_total"].astype(int),
        latent_right,
    )

    return {
        "left_alt": left_alt.astype(float),
        "left_total": depths["left_total"],
        "right_alt": right_alt.astype(float),
        "right_total": depths["right_total"],
    }


def simulation_schedule(number):
    combinations = list(
        itertools.product(GRID, GRID)
    )

    return [
        combinations[
            index % len(combinations)
        ]
        for index in range(number)
    ]


def simulation_preflight(
    included,
    donor,
    number,
):
    schedule = simulation_schedule(number)
    donor_lookup = donor.set_index("sample")

    donor_a = []
    donor_b = []
    high_call_pairs = 0

    for replicate in range(1, number + 1):
        first, second = choose_donors(
            included,
            replicate,
            None,
        )

        donor_a.append(first)
        donor_b.append(second)

        high_call_pairs += int(
            donor_lookup.loc[
                first,
                "high_call_sensitivity_set",
            ]
            and donor_lookup.loc[
                second,
                "high_call_sensitivity_set",
            ]
        )

    grid_counts = pd.Series(
        schedule
    ).value_counts()

    result = {
        "simulations_per_founder_model": number,
        "mixing_grid_pairs": len(set(schedule)),
        "minimum_replicates_per_grid_pair": int(
            grid_counts.min()
        ),
        "maximum_replicates_per_grid_pair": int(
            grid_counts.max()
        ),
        "distinct_donor_a_controls": len(
            set(donor_a)
        ),
        "distinct_donor_b_controls": len(
            set(donor_b)
        ),
        "high_call_control_count": int(
            donor[
                "high_call_sensitivity_set"
            ].sum()
        ),
        "both_donors_high_call_simulations_per_model": (
            high_call_pairs
        ),
    }

    if result["distinct_donor_a_controls"] != 47:
        raise RuntimeError(
            "Simulation schedule does not use "
            "every control as donor A."
        )

    if high_call_pairs < 10:
        raise RuntimeError(
            "Simulation schedule has fewer than "
            "10 high-call donor pairs per model; "
            "increase --replicates-per-model or "
            "review control callability."
        )

    return result


def calibration_table(
    simulations,
    observed_null,
    observed_alternative,
):
    rows = []

    for founder_model, fit_prefix, observed_fit in [
        (
            "three_founder",
            "null",
            observed_null,
        ),
        (
            "four_founder",
            "alternative",
            observed_alternative,
        ),
    ]:
        subset = simulations.loc[
            simulations[
                "true_founder_model"
            ] == founder_model
        ]

        for side in ["left", "right"]:
            true = subset[
                f"true_c_{side}"
            ].to_numpy(float)

            estimated = subset[
                f"{fit_prefix}_c_{side}"
            ].to_numpy(float)

            slope, intercept = np.polyfit(
                true,
                estimated,
                1,
            )

            if slope <= 0:
                raise RuntimeError(
                    "Non-positive mixing-calibration "
                    f"slope for {founder_model} "
                    f"{side}: {slope}"
                )

            predicted = intercept + slope * true
            residual = estimated - predicted

            calibrated = (
                observed_fit[f"c_{side}"]
                - intercept
            ) / slope
            calibrated = float(
                np.clip(calibrated, 0, 0.45)
            )

            residual_low, residual_high = (
                np.quantile(
                    residual,
                    [0.025, 0.975],
                )
            )

            ci_low = (
                observed_fit[f"c_{side}"]
                - intercept
                - residual_high
            ) / slope
            ci_high = (
                observed_fit[f"c_{side}"]
                - intercept
                - residual_low
            ) / slope

            rows.append({
                "founder_model": founder_model,
                "side": side,
                "observed_raw_estimate": (
                    observed_fit[f"c_{side}"]
                ),
                "calibration_intercept": intercept,
                "calibration_slope": slope,
                "simulation_bias": float(
                    np.mean(estimated - true)
                ),
                "simulation_RMSE": float(
                    np.sqrt(
                        np.mean(
                            (estimated - true) ** 2
                        )
                    )
                ),
                "calibrated_estimate": calibrated,
                "simulation_calibrated_95_low": (
                    float(
                        np.clip(
                            ci_low,
                            0,
                            0.45,
                        )
                    )
                ),
                "simulation_calibrated_95_high": (
                    float(
                        np.clip(
                            ci_high,
                            0,
                            0.45,
                        )
                    )
                ),
            })

    return pd.DataFrame(rows)


def plot_results(
    simulations,
    bootstrap,
    mixing,
    donor,
    observed_statistic,
    threshold,
):
    colors = {
        "three_founder": "#0072B2",
        "four_founder": "#D55E00",
    }

    with PdfPages(FIGURE_OUT) as pdf:
        figure, axes = plt.subplots(
            2,
            2,
            figsize=(10.5, 8.0),
        )

        ax = axes[0, 0]

        values = [
            simulations.loc[
                simulations[
                    "true_founder_model"
                ] == model,
                "founder_LR_statistic",
            ]
            for model in [
                "three_founder",
                "four_founder",
            ]
        ]

        parts = ax.violinplot(
            values,
            positions=[1, 2],
            showmedians=True,
            showextrema=False,
        )

        for body, model_name in zip(
            parts["bodies"],
            [
                "three_founder",
                "four_founder",
            ],
        ):
            body.set_facecolor(
                colors[model_name]
            )
            body.set_alpha(0.65)

        ax.axhline(
            observed_statistic,
            color="#000000",
            linewidth=1.5,
            label="Observed",
        )
        ax.axhline(
            threshold,
            color="#666666",
            linestyle="--",
            linewidth=1.1,
            label=(
                "95% three-founder threshold"
            ),
        )

        ax.set_xticks(
            [1, 2],
            [
                "Three-founder\nsimulations",
                "Four-founder\nsimulations",
            ],
        )
        ax.set_ylabel(
            "Nested founder likelihood-ratio statistic"
        )
        ax.set_title(
            "Founder-model calibration"
        )
        ax.legend(
            frameon=False,
            fontsize=7,
        )

        ax = axes[0, 1]

        for founder_model, marker in [
            ("three_founder", "o"),
            ("four_founder", "s"),
        ]:
            subset = simulations.loc[
                simulations[
                    "true_founder_model"
                ] == founder_model
            ]

            estimate_column = (
                "null_c_left"
                if founder_model == "three_founder"
                else "alternative_c_left"
            )

            ax.scatter(
                subset["true_c_left"],
                subset[estimate_column],
                s=10,
                alpha=0.35,
                marker=marker,
                color=colors[founder_model],
            )

        ax.plot(
            [0, 0.42],
            [0, 0.42],
            color="#333333",
            linewidth=1,
        )
        ax.set_xlim(-0.01, 0.42)
        ax.set_ylim(-0.01, 0.45)
        ax.set_xlabel(
            "Simulated mixing into Left"
        )
        ax.set_ylabel(
            "Estimated mixing into Left"
        )
        ax.set_title(
            "Mixing recovery: Left"
        )

        ax = axes[1, 0]

        for founder_model, marker in [
            ("three_founder", "o"),
            ("four_founder", "s"),
        ]:
            subset = simulations.loc[
                simulations[
                    "true_founder_model"
                ] == founder_model
            ]

            estimate_column = (
                "null_c_right"
                if founder_model == "three_founder"
                else "alternative_c_right"
            )

            ax.scatter(
                subset["true_c_right"],
                subset[estimate_column],
                s=10,
                alpha=0.35,
                marker=marker,
                color=colors[founder_model],
            )

        ax.plot(
            [0, 0.42],
            [0, 0.42],
            color="#333333",
            linewidth=1,
        )
        ax.set_xlim(-0.01, 0.42)
        ax.set_ylim(-0.01, 0.45)
        ax.set_xlabel(
            "Simulated mixing into Right"
        )
        ax.set_ylabel(
            "Estimated mixing into Right"
        )
        ax.set_title(
            "Mixing recovery: Right"
        )

        ax = axes[1, 1]

        for model_name, prefix in [
            ("three_founder", "null"),
            ("four_founder", "alternative"),
        ]:
            ax.scatter(
                bootstrap[
                    f"{prefix}_c_left_calibrated"
                ],
                bootstrap[
                    f"{prefix}_c_right_calibrated"
                ],
                s=15,
                alpha=0.35,
                color=colors[model_name],
                label=model_name.replace(
                    "_",
                    " ",
                ),
            )

        ax.plot(
            [0, 0.45],
            [0, 0.45],
            color="#555555",
            linewidth=0.8,
            linestyle="--",
        )
        ax.set_xlim(0, 0.45)
        ax.set_ylim(0, 0.45)
        ax.set_xlabel(
            "Calibrated bootstrap mixing into Left"
        )
        ax.set_ylabel(
            "Calibrated bootstrap mixing into Right"
        )
        ax.set_title(
            "Genomic-block uncertainty"
        )
        ax.legend(
            frameon=False,
            fontsize=8,
        )

        figure.tight_layout()

        pdf.savefig(
            figure,
            bbox_inches="tight",
        )
        figure.savefig(
            FIGURE_PNG,
            dpi=240,
            bbox_inches="tight",
        )

        plt.close(figure)

        figure, axes = plt.subplots(
            1,
            2,
            figsize=(11.2, 4.5),
        )

        ax = axes[0]

        point_colors = np.where(
            donor["high_call_sensitivity_set"],
            "#009E73",
            "#999999",
        )

        ax.scatter(
            donor["call_rate"],
            donor["heterozygosity"],
            c=point_colors,
            s=30,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.4,
        )

        ax.axvline(
            HIGH_CALL_RATE,
            color="#555555",
            linestyle="--",
            linewidth=1,
        )

        annotated = pd.concat([
            donor.nsmallest(
                2,
                "call_rate",
            ),
            donor.nlargest(
                2,
                "heterozygosity",
            ),
        ]).drop_duplicates("sample")

        for _, row in annotated.iterrows():
            ax.annotate(
                row["sample"],
                (
                    row["call_rate"],
                    row["heterozygosity"],
                ),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=6,
            )

        ax.set_xlabel(
            "Callable fraction of 8,152 sites"
        )
        ax.set_ylabel(
            "Heterozygosity among callable sites"
        )
        ax.set_title(
            "All 47 simulation controls"
        )

        ax = axes[1]

        for founder_model in [
            "three_founder",
            "four_founder",
        ]:
            subset = simulations.loc[
                simulations[
                    "true_founder_model"
                ] == founder_model
            ]

            ax.scatter(
                subset[
                    "true_opposing_homozygote_fraction"
                ],
                subset[
                    "founder_LR_statistic"
                ],
                s=14,
                alpha=0.45,
                color=colors[founder_model],
                label=founder_model.replace(
                    "_",
                    " ",
                ),
            )

        ax.axhline(
            observed_statistic,
            color="#000000",
            linewidth=1.3,
            label="observed statistic",
        )
        ax.set_xlabel(
            "True opposing-homozygote fraction"
        )
        ax.set_ylabel(
            "Nested founder likelihood-ratio statistic"
        )
        ax.set_title(
            "Signal defining the founder comparison"
        )
        ax.legend(
            frameon=False,
            fontsize=7,
        )

        figure.tight_layout()
        pdf.savefig(
            figure,
            bbox_inches="tight",
        )
        plt.close(figure)

        figure, axes = plt.subplots(
            1,
            2,
            figsize=(11.2, 4.5),
        )

        for ax, founder_model, prefix in [
            (
                axes[0],
                "three_founder",
                "null",
            ),
            (
                axes[1],
                "four_founder",
                "alternative",
            ),
        ]:
            subset = simulations.loc[
                simulations[
                    "true_founder_model"
                ] == founder_model
            ].copy()

            subset["true_mean"] = (
                subset["true_c_left"]
                + subset["true_c_right"]
            ) / 2

            subset["estimated_mean"] = (
                subset[f"{prefix}_c_left"]
                + subset[f"{prefix}_c_right"]
            ) / 2

            subset["bias"] = (
                subset["estimated_mean"]
                - subset["true_mean"]
            )

            grid = subset.pivot_table(
                index="true_c_right",
                columns="true_c_left",
                values="bias",
                aggfunc="mean",
            ).reindex(
                index=GRID,
                columns=GRID,
            )

            image = ax.imshow(
                grid.to_numpy(),
                origin="lower",
                aspect="auto",
                extent=[
                    GRID.min(),
                    GRID.max(),
                    GRID.min(),
                    GRID.max(),
                ],
                cmap="BrBG",
                vmin=-0.08,
                vmax=0.08,
            )

            ax.set_xlabel(
                "True mixing into Left"
            )
            ax.set_ylabel(
                "True mixing into Right"
            )
            ax.set_title(
                founder_model.replace(
                    "_",
                    " ",
                ).title()
            )

        colorbar = figure.colorbar(
            image,
            ax=axes.ravel().tolist(),
            shrink=0.85,
        )
        colorbar.set_label(
            "Mean estimated-minus-true mixing"
        )

        figure.suptitle(
            "Mixing-estimation bias across "
            "the simulation grid",
            y=0.99,
        )
        figure.subplots_adjust(
            left=0.08,
            right=0.88,
            bottom=0.12,
            top=0.88,
            wspace=0.25,
        )

        pdf.savefig(
            figure,
            bbox_inches="tight",
        )
        plt.close(figure)


def main():
    args = arguments()

    if (
        args.replicates_per_model < 47
        or args.bootstrap_replicates < 10
    ):
        raise RuntimeError(
            "Use at least 47 simulations/model "
            "so every control is represented, "
            "and at least 10 bootstraps."
        )

    for path in [
        MODEL_SCRIPT,
        MATCHED_INPUT,
        FITS_INPUT,
        CONTROL_INPUT,
        AUDIT_INPUT,
    ]:
        if not path.is_file():
            raise RuntimeError(
                f"Missing required input: {path}"
            )

    model = load_model_module()
    emission = model.load_emission()

    matched = pd.read_csv(
        MATCHED_INPUT,
        sep="\t",
        low_memory=False,
    )

    required = set(KEYS + [
        "left_alt_depth",
        "left_AD_total",
        "right_alt_depth",
        "right_AD_total",
        "block_id",
    ])

    missing = required - set(matched.columns)

    if missing:
        raise RuntimeError(
            "Matched input missing columns: "
            f"{sorted(missing)}"
        )

    if matched["chrom"].isin(
        SEX_CHROMS
    ).any():
        raise RuntimeError(
            "A sex-chromosome row survived "
            "into matched input."
        )

    arrays = {
        "left_alt": matched[
            "left_alt_depth"
        ].to_numpy(float),
        "left_total": matched[
            "left_AD_total"
        ].to_numpy(float),
        "right_alt": matched[
            "right_alt_depth"
        ].to_numpy(float),
        "right_total": matched[
            "right_AD_total"
        ].to_numpy(float),
    }

    (
        matrix,
        donor,
        included,
        allele_frequency,
    ) = build_control_matrix(matched)

    donor.to_csv(
        DONOR_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
    )

    preflight = simulation_preflight(
        included,
        donor,
        args.replicates_per_model,
    )

    print(
        f"Controls included in simulation: "
        f"{len(included)} of 47"
    )
    print(
        included.groupby(
            "heterozygosity_stratum"
        ).size().to_string()
    )
    print(
        "\nSIMULATION SCHEDULE PREFLIGHT"
    )
    print(
        pd.Series(preflight).to_string()
    )

    if args.prepare_only:
        print(
            f"Prepared donor audit: {DONOR_OUT}"
        )
        return

    initial_fits = pd.read_csv(
        FITS_INPUT,
        sep="\t",
    ).set_index("model")

    c3 = initial_fits.loc[
        "three_founder"
    ]
    c4 = initial_fits.loc[
        "four_founder"
    ]

    null_starts = [
        (
            float(c3["c_left_from_right"]),
            float(c3["c_right_from_left"]),
        ),
        (0.10, 0.10),
        (0.25, 0.25),
    ]

    alternative_starts = [
        (
            float(c4["c_left_from_right"]),
            float(c4["c_right_from_left"]),
            0.01,
        ),
        (0.15, 0.15, 0.05),
        (0.25, 0.25, 0.15),
    ]

    (
        observed_null,
        observed_alt,
        observed_statistic,
    ) = fitted_pair(
        model,
        arrays,
        emission,
        null_starts,
        alternative_starts,
    )

    observed_table = pd.DataFrame([
        {
            "fit": "three_founder_null",
            **observed_null,
            "LR_statistic_vs_null": 0.0,
        },
        {
            "fit": (
                "opposing_state_weight_alternative"
            ),
            **observed_alt,
            "LR_statistic_vs_null": (
                observed_statistic
            ),
        },
    ])

    observed_table.to_csv(
        OBSERVED_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
    )

    print("\nOBSERVED NESTED FIT")
    print(
        observed_table.to_string(index=False)
    )

    rng = np.random.default_rng(SEED)
    donor_lookup = donor.set_index("sample")
    simulation_rows = []

    schedule = simulation_schedule(
        args.replicates_per_model
    )

    for founder_model in [
        "three_founder",
        "four_founder",
    ]:
        for replicate, (
            true_c_left,
            true_c_right,
        ) in enumerate(
            schedule,
            start=1,
        ):
            (
                true_left,
                true_right,
                donor_a,
                donor_b,
                construction,
            ) = true_genotypes(
                founder_model,
                matrix,
                included,
                allele_frequency,
                replicate,
                rng,
            )

            simulated = simulate_arrays(
                model,
                true_left,
                true_right,
                true_c_left,
                true_c_right,
                arrays,
                emission,
                rng,
            )

            local_null_starts = [
                (0.20, 0.20)
            ]
            local_alt_starts = [
                (0.20, 0.20, 0.03)
            ]

            (
                null,
                alternative,
                statistic,
            ) = fitted_pair(
                model,
                simulated,
                emission,
                local_null_starts,
                local_alt_starts,
            )

            simulation_rows.append({
                "true_founder_model": (
                    founder_model
                ),
                "replicate": replicate,
                "true_c_left": true_c_left,
                "true_c_right": true_c_right,
                "donor_a": donor_a,
                "donor_b": donor_b,
                "donor_a_heterozygosity": (
                    donor_lookup.loc[
                        donor_a,
                        "heterozygosity",
                    ]
                ),
                "donor_b_heterozygosity": (
                    donor_lookup.loc[
                        donor_b,
                        "heterozygosity",
                    ]
                ),
                "donor_a_call_rate": (
                    donor_lookup.loc[
                        donor_a,
                        "call_rate",
                    ]
                ),
                "donor_b_call_rate": (
                    donor_lookup.loc[
                        donor_b,
                        "call_rate",
                    ]
                ),
                "both_donors_high_call": bool(
                    donor_lookup.loc[
                        donor_a,
                        "high_call_sensitivity_set",
                    ]
                    and donor_lookup.loc[
                        donor_b,
                        "high_call_sensitivity_set",
                    ]
                ),
                "construction": construction,
                "true_opposing_homozygote_fraction": (
                    float(
                        np.mean(
                            (
                                (true_left == 0)
                                & (true_right == 2)
                            )
                            | (
                                (true_left == 2)
                                & (true_right == 0)
                            )
                        )
                    )
                ),
                "null_c_left": null["c_left"],
                "null_c_right": null["c_right"],
                "alternative_c_left": (
                    alternative["c_left"]
                ),
                "alternative_c_right": (
                    alternative["c_right"]
                ),
                "estimated_opposing_state_weight": (
                    alternative["epsilon"]
                ),
                "founder_LR_statistic": statistic,
            })

            if replicate % 10 == 0:
                print(
                    f"{founder_model}: completed "
                    f"{replicate} simulations",
                    flush=True,
                )

    simulations = pd.DataFrame(
        simulation_rows
    )

    for founder_model in [
        "three_founder",
        "four_founder",
    ]:
        used = set(
            simulations.loc[
                simulations[
                    "true_founder_model"
                ] == founder_model,
                "donor_a",
            ]
        )

        missing_donors = sorted(
            set(
                included[
                    "sample"
                ].astype(str)
            ) - used
        )

        if missing_donors:
            raise RuntimeError(
                "Not all 47 controls were used "
                f"under {founder_model}: "
                f"{missing_donors}"
            )

    distances = np.empty(
        len(simulations)
    )

    for index, row in simulations.iterrows():
        reference = (
            observed_null
            if row[
                "true_founder_model"
            ] == "three_founder"
            else observed_alt
        )

        distances[index] = math.hypot(
            row["true_c_left"]
            - reference["c_left"],
            row["true_c_right"]
            - reference["c_right"],
        )

    simulations[
        "distance_from_corresponding_"
        "observed_mixing"
    ] = distances

    simulations.to_csv(
        SIM_OUT,
        sep="\t",
        index=False,
        compression="gzip",
        float_format=FLOAT_FORMAT,
    )

    nearest_parts = []

    for founder_model in [
        "three_founder",
        "four_founder",
    ]:
        nearest_parts.append(
            simulations.loc[
                simulations[
                    "true_founder_model"
                ] == founder_model
            ].nsmallest(
                min(
                    NEAREST_MAIN,
                    args.replicates_per_model,
                ),
                (
                    "distance_from_corresponding_"
                    "observed_mixing"
                ),
            )
        )

    nearest = pd.concat(
        nearest_parts,
        ignore_index=True,
    )

    nearest_m3 = nearest.loc[
        nearest[
            "true_founder_model"
        ] == "three_founder"
    ]
    nearest_m4 = nearest.loc[
        nearest[
            "true_founder_model"
        ] == "four_founder"
    ]

    m3_threshold = float(
        nearest_m3[
            "founder_LR_statistic"
        ].quantile(0.95)
    )

    empirical_p = float(
        (
            1
            + (
                nearest_m3[
                    "founder_LR_statistic"
                ] >= observed_statistic
            ).sum()
        )
        / (1 + len(nearest_m3))
    )

    m4_detection = float(
        (
            nearest_m4[
                "founder_LR_statistic"
            ] > m3_threshold
        ).mean()
    )

    high_call = simulations.loc[
        simulations[
            "both_donors_high_call"
        ]
    ].copy()

    high_call_nearest_parts = []

    for founder_model in [
        "three_founder",
        "four_founder",
    ]:
        part = high_call.loc[
            high_call[
                "true_founder_model"
            ] == founder_model
        ]

        high_call_nearest_parts.append(
            part.nsmallest(
                min(
                    NEAREST_HIGH,
                    len(part),
                ),
                (
                    "distance_from_corresponding_"
                    "observed_mixing"
                ),
            )
        )

    high_call_nearest = pd.concat(
        high_call_nearest_parts,
        ignore_index=True,
    )

    high_m3 = high_call_nearest.loc[
        high_call_nearest[
            "true_founder_model"
        ] == "three_founder"
    ]
    high_m4 = high_call_nearest.loc[
        high_call_nearest[
            "true_founder_model"
        ] == "four_founder"
    ]

    if (
        len(high_m3) < 10
        or len(high_m4) < 10
    ):
        raise RuntimeError(
            "Insufficient high-call "
            "sensitivity simulations."
        )

    high_threshold = float(
        high_m3[
            "founder_LR_statistic"
        ].quantile(0.95)
    )

    high_empirical_p = float(
        (
            1
            + (
                high_m3[
                    "founder_LR_statistic"
                ] >= observed_statistic
            ).sum()
        )
        / (1 + len(high_m3))
    )

    high_m4_detection = float(
        (
            high_m4[
                "founder_LR_statistic"
            ] > high_threshold
        ).mean()
    )

    summary_rows = []

    for founder_model in [
        "three_founder",
        "four_founder",
    ]:
        subset = simulations.loc[
            simulations[
                "true_founder_model"
            ] == founder_model
        ]

        close = nearest.loc[
            nearest[
                "true_founder_model"
            ] == founder_model
        ]

        high_close = high_call_nearest.loc[
            high_call_nearest[
                "true_founder_model"
            ] == founder_model
        ]

        for scope, values in [
            (
                "all_mixing_values",
                subset,
            ),
            (
                (
                    f"nearest_{NEAREST_MAIN}_"
                    "all_47_controls"
                ),
                close,
            ),
            (
                "nearest_high_call_sensitivity",
                high_close,
            ),
        ]:
            summary_rows.append({
                "true_founder_model": (
                    founder_model
                ),
                "scope": scope,
                "replicates": len(values),
                "LR_median": values[
                    "founder_LR_statistic"
                ].median(),
                "LR_95th_percentile": values[
                    "founder_LR_statistic"
                ].quantile(0.95),
                "epsilon_median": values[
                    "estimated_opposing_state_weight"
                ].median(),
                "true_opposing_fraction_median": (
                    values[
                        "true_opposing_"
                        "homozygote_fraction"
                    ].median()
                ),
            })

    pd.DataFrame(
        summary_rows
    ).to_csv(
        SIM_SUMMARY_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
    )

    mixing = calibration_table(
        simulations,
        observed_null,
        observed_alt,
    )

    mixing.to_csv(
        MIXING_OUT,
        sep="\t",
        index=False,
        float_format=FLOAT_FORMAT,
    )

    block_ids = matched[
        "block_id"
    ].astype(str).to_numpy()

    unique_blocks = np.unique(
        block_ids
    )

    bootstrap_rows = []

    for replicate in range(
        1,
        args.bootstrap_replicates + 1,
    ):
        sampled = rng.choice(
            unique_blocks,
            size=len(unique_blocks),
            replace=True,
        )

        block_counts = pd.Series(
            sampled
        ).value_counts().to_dict()

        weights = np.array([
            block_counts.get(block, 0)
            for block in block_ids
        ], dtype=float)

        (
            null,
            alternative,
            statistic,
        ) = fitted_pair(
            model,
            arrays,
            emission,
            [
                (
                    observed_null["c_left"],
                    observed_null["c_right"],
                )
            ],
            [
                (
                    observed_alt["c_left"],
                    observed_alt["c_right"],
                    observed_alt["epsilon"],
                )
            ],
            weights,
        )

        bootstrap_rows.append({
            "replicate": replicate,
            "null_c_left": null["c_left"],
            "null_c_right": null["c_right"],
            "alternative_c_left": (
                alternative["c_left"]
            ),
            "alternative_c_right": (
                alternative["c_right"]
            ),
            "alternative_epsilon": (
                alternative["epsilon"]
            ),
            "founder_LR_statistic": statistic,
        })

        if replicate % 10 == 0:
            print(
                "block bootstrap: completed "
                f"{replicate}",
                flush=True,
            )

    bootstrap = pd.DataFrame(
        bootstrap_rows
    )

    for _, calibration in mixing.iterrows():
        founder_model = calibration[
            "founder_model"
        ]
        side = calibration["side"]

        prefix = (
            "null"
            if founder_model == "three_founder"
            else "alternative"
        )

        raw_column = (
            f"{prefix}_c_{side}"
        )
        calibrated_column = (
            f"{raw_column}_calibrated"
        )

        bootstrap[
            calibrated_column
        ] = np.clip(
            (
                bootstrap[raw_column]
                - calibration[
                    "calibration_intercept"
                ]
            )
            / calibration[
                "calibration_slope"
            ],
            0,
            0.45,
        )

    bootstrap.to_csv(
        BOOTSTRAP_OUT,
        sep="\t",
        index=False,
        compression="gzip",
        float_format=FLOAT_FORMAT,
    )

    conclusion = (
        (
            "reject_three_founder_in_favor_"
            "of_four_founder_signal"
        )
        if observed_statistic > m3_threshold
        else "three_founder_not_rejected"
    )

    final_rows = [
        (
            "observed_nested_LR_statistic",
            observed_statistic,
        ),
        (
            "observed_opposing_state_weight",
            observed_alt["epsilon"],
        ),
        (
            "nearest_three_founder_simulations",
            len(nearest_m3),
        ),
        (
            (
                "three_founder_simulated_"
                "95_percent_threshold"
            ),
            m3_threshold,
        ),
        (
            "three_founder_empirical_p_value",
            empirical_p,
        ),
        (
            (
                "four_founder_detection_"
                "fraction_at_threshold"
            ),
            m4_detection,
        ),
        (
            (
                "high_call_three_founder_"
                "simulated_95_percent_threshold"
            ),
            high_threshold,
        ),
        (
            (
                "high_call_three_founder_"
                "empirical_p_value"
            ),
            high_empirical_p,
        ),
        (
            (
                "high_call_four_founder_"
                "detection_fraction_at_threshold"
            ),
            high_m4_detection,
        ),
        (
            "calibrated_founder_conclusion",
            conclusion,
        ),
        (
            "founder_support_interpretation",
            (
                "simulation_calibrated_nested_LR_"
                "test_not_posterior_model_probability"
            ),
        ),
        (
            "inference_scope",
            "Hmr18_optix_matched_sites_only",
        ),
        (
            "observed_three_founder_c_left",
            observed_null["c_left"],
        ),
        (
            "observed_three_founder_c_right",
            observed_null["c_right"],
        ),
        (
            "observed_four_founder_c_left",
            observed_alt["c_left"],
        ),
        (
            "observed_four_founder_c_right",
            observed_alt["c_right"],
        ),
        (
            "block_bootstrap_replicates",
            len(bootstrap),
        ),
        (
            "simulation_replicates_per_model",
            args.replicates_per_model,
        ),
        (
            "simulation_donors",
            len(included),
        ),
        (
            "all_47_controls_used_under_each_model",
            True,
        ),
        (
            "T2_T3_included",
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
    ]

    for _, row in mixing.iterrows():
        label = (
            f"{row['founder_model']}_"
            f"c_{row['side']}"
        )

        final_rows += [
            (
                f"calibrated_{label}",
                row["calibrated_estimate"],
            ),
            (
                (
                    f"simulation_calibrated_{label}_"
                    "2.5_percent"
                ),
                row[
                    "simulation_calibrated_95_low"
                ],
            ),
            (
                (
                    f"simulation_calibrated_{label}_"
                    "97.5_percent"
                ),
                row[
                    "simulation_calibrated_95_high"
                ],
            ),
        ]

    for column in [
        "null_c_left",
        "null_c_right",
        "alternative_c_left",
        "alternative_c_right",
        "alternative_epsilon",
        "null_c_left_calibrated",
        "null_c_right_calibrated",
        "alternative_c_left_calibrated",
        "alternative_c_right_calibrated",
    ]:
        low, median, high = bootstrap[
            column
        ].quantile([
            0.025,
            0.5,
            0.975,
        ])

        final_rows += [
            (
                f"bootstrap_{column}_2.5_percent",
                low,
            ),
            (
                f"bootstrap_{column}_median",
                median,
            ),
            (
                f"bootstrap_{column}_97.5_percent",
                high,
            ),
        ]

    final = pd.DataFrame(
        final_rows,
        columns=["statistic", "value"],
    )

    final.to_csv(
        FINAL_OUT,
        sep="\t",
        index=False,
    )

    plot_results(
        simulations,
        bootstrap,
        mixing,
        donor,
        observed_statistic,
        m3_threshold,
    )

    outputs = [
        OBSERVED_OUT,
        SIM_OUT,
        SIM_SUMMARY_OUT,
        MIXING_OUT,
        BOOTSTRAP_OUT,
        FINAL_OUT,
        DONOR_OUT,
        FIGURE_OUT,
        FIGURE_PNG,
    ]

    CHECKSUM_OUT.write_text(
        "\n".join(
            (
                f"{sha256(path)}  "
                f"{path.relative_to(PROJECT)}"
            )
            for path in outputs
        ) + "\n",
        encoding="utf-8",
    )

    print(
        "\nFINAL CALIBRATED SUMMARY"
    )
    print(
        final.to_string(index=False)
    )

    print(
        "\nMIXING CALIBRATION"
    )
    print(
        mixing.to_string(index=False)
    )

    print(f"\nFigure: {FIGURE_OUT}")
    print(
        f"Runtime: Python "
        f"{platform.python_version()}, "
        f"NumPy {np.__version__}, "
        f"pandas {pd.__version__}, "
        f"SciPy {scipy.__version__}"
    )


if __name__ == "__main__":
    main()
