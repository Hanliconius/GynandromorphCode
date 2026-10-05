#!/usr/bin/env python3

import argparse
import csv
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("analysis/genomewide_mixing_corrected_v1")
BASE = ROOT / "expanded_calibration"
OUT = ROOT / "empirical_full_sib_calibration_v1"
TASKS = OUT / "tasks"
FINAL = OUT / "final"
REPRESENTATIVES = OUT / "representative_datasets"
FIGURES = OUT / "figures"

OBSERVED_JSON = BASE / "observed_expanded.json"
RESULTS = FINAL / "all_empirical_sibling_simulations.tsv.gz"
SUMMARY = FINAL / "empirical_sibling_simulation_summary.tsv"
REPRESENTATIVE_TABLE = (
    FINAL
    / "representative_simulations.tsv"
)
WINDOW_EVIDENCE = (
    FINAL
    / "observed_and_simulated_1Mb_evidence.tsv"
)


def load_module(path, name):
    specification = importlib.util.spec_from_file_location(
        name,
        path,
    )

    module = importlib.util.module_from_spec(
        specification
    )

    specification.loader.exec_module(module)

    return module


stage201 = load_module(
    "scripts/201_simulate_empirical_full_sib_mixing.py",
    "stage201_empirical_siblings",
)

pilot = stage201.pilot

BLUE = "#0072B2"
LIGHT_BLUE = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#D55E00"
GREY = "#666666"
EDGES = np.linspace(0, 1, 61)


def observed_fit():
    document = json.loads(
        OBSERVED_JSON.read_text()
    )

    fit = document["full"]

    for key in (
        "null",
        "alternative",
        "LR",
    ):
        if key not in fit:
            raise RuntimeError(
                f"Observed fit lacks {key}"
            )

    return fit


def collect():
    FINAL.mkdir(
        parents=True,
        exist_ok=True,
    )

    pieces = []
    missing = []

    for task in range(660):
        result = (
            TASKS
            / f"task_{task:03d}.tsv"
        )

        completion = (
            TASKS
            / f"task_{task:03d}.complete.tsv"
        )

        if (
            not result.is_file()
            or result.stat().st_size == 0
            or not completion.is_file()
            or completion.stat().st_size == 0
        ):
            missing.append(task)
            continue

        table = pd.read_csv(
            result,
            sep="\t",
        )

        if len(table) != 2:
            raise RuntimeError(
                f"Task {task}: expected two rows; "
                f"found {len(table)}"
            )

        expected_models = {
            "three_founder",
            "empirical_full_sib",
        }

        if set(table["generating_model"]) != expected_models:
            raise RuntimeError(
                f"Task {task}: generating-model rows "
                "are incorrect"
            )

        if set(table["task"]) != {task}:
            raise RuntimeError(
                f"Task {task}: task identifier mismatch"
            )

        pieces.append(table)

    if missing:
        raise RuntimeError(
            "Missing or incomplete simulation tasks: "
            + ",".join(
                str(task)
                for task in missing
            )
        )

    results = pd.concat(
        pieces,
        ignore_index=True,
    )

    if len(results) != 1320:
        raise RuntimeError(
            f"Expected 1,320 rows; "
            f"found {len(results)}"
        )

    numeric_columns = [
        "true_c_left",
        "true_c_right",
        "genotype_contrast",
        "IBD0_fraction",
        "IBD1_fraction",
        "IBD2_fraction",
        "estimated_c_left_null",
        "estimated_c_right_null",
        "estimated_c_left_alt",
        "estimated_c_right_alt",
        "epsilon",
        "LR",
        "sites",
    ]

    if not np.isfinite(
        results[numeric_columns].to_numpy(
            dtype=float
        )
    ).all():
        raise RuntimeError(
            "Non-finite simulation result detected"
        )

    if set(results["sites"]) != {1_000_000}:
        raise RuntimeError(
            "Not every simulation used exactly "
            "1,000,000 sites"
        )

    expected_counts = (
        results.groupby(
            [
                "scenario",
                "generating_model",
            ]
        )
        .size()
    )

    for scenario in (
        "exact_20",
        "observed_matched",
    ):
        for model_name in (
            "three_founder",
            "empirical_full_sib",
        ):
            count = int(
                expected_counts.loc[
                    scenario,
                    model_name,
                ]
            )

            if count != 330:
                raise RuntimeError(
                    f"{scenario}/{model_name}: "
                    f"expected 330 rows; found {count}"
                )

    results.to_csv(
        RESULTS,
        sep="\t",
        index=False,
        compression="gzip",
    )

    observed_lr = float(
        observed_fit()["LR"]
    )

    summary_rows = []

    for (
        scenario,
        generating_model,
    ), group in results.groupby(
        [
            "scenario",
            "generating_model",
        ],
        sort=True,
    ):
        values = group["LR"].to_numpy()

        summary_rows.append(
            {
                "scenario": scenario,
                "generating_model":
                    generating_model,
                "simulations": len(group),
                "sites_per_simulation":
                    int(group["sites"].iloc[0]),
                "mean_LR":
                    float(np.mean(values)),
                "median_LR":
                    float(np.median(values)),
                "LR_2.5_percentile":
                    float(np.quantile(values, 0.025)),
                "LR_97.5_percentile":
                    float(np.quantile(values, 0.975)),
                "minimum_LR":
                    float(np.min(values)),
                "maximum_LR":
                    float(np.max(values)),
                "observed_LR": observed_lr,
                "fraction_simulated_LR_below_observed":
                    float(np.mean(values < observed_lr)),
                "fraction_simulated_LR_at_or_above_observed":
                    float(np.mean(values >= observed_lr)),
            }
        )

    pd.DataFrame(summary_rows).to_csv(
        SUMMARY,
        sep="\t",
        index=False,
    )

    audit = {
        "status": "PASS",
        "simulation_tasks": 660,
        "simulation_rows": len(results),
        "generating_models": [
            "three_founder",
            "empirical_full_sib",
        ],
        "scenarios": [
            "exact_20",
            "observed_matched",
        ],
        "rows_per_scenario_and_model": 330,
        "sites_per_simulation": 1_000_000,
        "sibling_pairs": 66,
        "replicates_per_pair_and_scenario": 5,
        "observed_LR": observed_lr,
        "likelihood_changed": False,
        "emission_model_changed": False,
    }

    (
        FINAL
        / "collection.audit.json"
    ).write_text(
        json.dumps(
            audit,
            indent=2,
        )
        + "\n"
    )

    (
        FINAL
        / "collection.complete.tsv"
    ).write_text(
        "status\ttasks\trows\tsites_per_simulation\n"
        "PASS\t660\t1320\t1000000\n"
    )

    print(
        "EMPIRICAL-SIBLING SIMULATION "
        "COLLECTION COMPLETE"
    )

    print(
        pd.DataFrame(summary_rows).to_string(
            index=False
        )
    )


def reconstruct(row, data):
    schedule = (
        stage201.simulation_schedule()
    )

    task = int(row["task"])
    configuration = schedule[task]

    left_index = int(
        configuration[
            "offspring_a_index"
        ]
    )

    right_index = int(
        configuration[
            "offspring_b_index"
        ]
    )

    left_sample = configuration["offspring_a"]
    right_sample = configuration["offspring_b"]

    if int(
        configuration["repetition"]
    ) % 2 == 1:
        left_index, right_index = (
            right_index,
            left_index,
        )

        left_sample, right_sample = (
            right_sample,
            left_sample,
        )

    controls = data["controls"]
    called = controls >= 0

    allele_frequency = (
        np.where(
            called,
            controls,
            0,
        ).sum(axis=1)
        + 0.5
    ) / (
        2 * called.sum(axis=1)
        + 1
    )

    generating_model = row["generating_model"]

    founder_index = (
        0
        if generating_model
        == "three_founder"
        else 1
    )

    rng = np.random.default_rng(
        201000000
        + task * 2
        + founder_index
    )

    sibling_left = stage201.impute(
        data["siblings"][:, left_index],
        allele_frequency,
        rng,
    )

    sibling_right = stage201.impute(
        data["siblings"][:, right_index],
        allele_frequency,
        rng,
    )

    if generating_model == "three_founder":
        constructed = (
            stage201.gamete(
                sibling_left,
                rng,
            )
            + stage201.gamete(
                sibling_right,
                rng,
            )
        )

        if task % 2 == 0:
            true_left = sibling_left
            true_right = constructed
        else:
            true_left = constructed
            true_right = sibling_left
    else:
        true_left = sibling_left
        true_right = sibling_right

    synthetic = pilot.simulate_counts(
        true_left,
        true_right,
        float(row["true_c_left"]),
        float(row["true_c_right"]),
        data,
        pilot.load_emission(),
        rng,
    )

    return (
        synthetic,
        true_left,
        true_right,
        left_sample,
        right_sample,
    )


def select_representatives(results):
    rows = []

    for scenario in (
        "exact_20",
        "observed_matched",
    ):
        for generating_model in (
            "three_founder",
            "empirical_full_sib",
        ):
            subset = results.loc[
                (
                    results["scenario"]
                    == scenario
                )
                & (
                    results["generating_model"]
                    == generating_model
                )
            ].copy()

            median_lr = float(
                subset["LR"].median()
            )

            row = subset.loc[
                (
                    subset["LR"]
                    - median_lr
                ).abs().idxmin()
            ].copy()

            row["selection_label"] = (
                f"{scenario}_"
                f"{generating_model}_median"
            )

            rows.append(row)

    full_sibling = results.loc[
        (
            results["scenario"]
            == "exact_20"
        )
        & (
            results["generating_model"]
            == "empirical_full_sib"
        )
    ].copy()

    target = float(
        observed_fit()["LR"]
    )

    closest = full_sibling.loc[
        (
            full_sibling["LR"]
            - target
        ).abs().idxmin()
    ].copy()

    closest["selection_label"] = (
        "exact_20_empirical_full_sib_"
        "closest_to_observed"
    )

    rows.append(closest)

    return pd.DataFrame(rows)


def site_contributions(
    data,
    null_fit,
    alternative_fit,
):
    emission = pilot.load_emission()

    null_states = (
        pilot.model.state_log_likelihoods(
            data,
            float(null_fit["c_left"]),
            float(null_fit["c_right"]),
            emission,
        )
    )

    alternative_states = (
        pilot.model.state_log_likelihoods(
            data,
            float(
                alternative_fit["c_left"]
            ),
            float(
                alternative_fit["c_right"]
            ),
            emission,
        )
    )

    null_evidence = pilot.model_evidence(
        null_states,
        0.0,
    )

    alternative_evidence = (
        pilot.model_evidence(
            alternative_states,
            float(
                alternative_fit["epsilon"]
            ),
        )
    )

    return 2.0 * (
        alternative_evidence
        - null_evidence
    )


def window_rows(
    label,
    source_type,
    data,
    contributions,
):
    rows = []

    for chromosome in range(1, 21):
        chromosome_mask = (
            data["chrom"]
            == chromosome
        )

        positions = (
            data["pos"][chromosome_mask]
        )

        values = (
            contributions[chromosome_mask]
        )

        starts = (
            (positions - 1)
            // 1_000_000
            * 1_000_000
        )

        for start in np.unique(starts):
            mask = starts == start
            sites = int(mask.sum())
            contribution = float(
                values[mask].sum()
            )

            rows.append(
                {
                    "source_label": label,
                    "source_type": source_type,
                    "chromosome": chromosome,
                    "window_start":
                        int(start),
                    "window_end":
                        int(start + 1_000_000),
                    "sites": sites,
                    "LR_contribution":
                        contribution,
                    "LR_contribution_per_site":
                        contribution / sites,
                }
            )

    return rows


def evidence():
    FINAL.mkdir(
        parents=True,
        exist_ok=True,
    )

    REPRESENTATIVES.mkdir(
        parents=True,
        exist_ok=True,
    )

    results = pd.read_csv(
        RESULTS,
        sep="\t",
    )

    representatives = (
        select_representatives(results)
    )

    representatives.to_csv(
        REPRESENTATIVE_TABLE,
        sep="\t",
        index=False,
    )

    data = stage201.load_input()
    all_window_rows = []

    observed = observed_fit()

    observed_contributions = site_contributions(
        data,
        observed["null"],
        observed["alternative"],
    )

    all_window_rows.extend(
        window_rows(
            "observed_LR",
            "observed",
            data,
            observed_contributions,
        )
    )

    observed_sum = float(
        observed_contributions.sum()
    )

    if not np.isclose(
        observed_sum,
        float(observed["LR"]),
        rtol=1e-8,
        atol=1e-5,
    ):
        raise RuntimeError(
            "Observed site contributions do not "
            "sum to the observed LR: "
            f"{observed_sum} versus "
            f"{observed['LR']}"
        )

    reconstruction_rows = []

    for _, row in representatives.iterrows():
        (
            synthetic,
            true_left,
            true_right,
            left_sample,
            right_sample,
        ) = reconstruct(
            row,
            data,
        )

        label = row["selection_label"]

        null_fit = {
            "c_left":
                float(
                    row[
                        "estimated_c_left_null"
                    ]
                ),
            "c_right":
                float(
                    row[
                        "estimated_c_right_null"
                    ]
                ),
        }

        alternative_fit = {
            "c_left":
                float(
                    row[
                        "estimated_c_left_alt"
                    ]
                ),
            "c_right":
                float(
                    row[
                        "estimated_c_right_alt"
                    ]
                ),
            "epsilon":
                float(row["epsilon"]),
        }

        contributions = site_contributions(
            synthetic,
            null_fit,
            alternative_fit,
        )

        contribution_sum = float(
            contributions.sum()
        )

        if not np.isclose(
            contribution_sum,
            float(row["LR"]),
            rtol=1e-8,
            atol=1e-5,
        ):
            raise RuntimeError(
                f"{label}: reconstructed site "
                f"contributions sum to "
                f"{contribution_sum}, expected "
                f"{row['LR']}"
            )

        all_window_rows.extend(
            window_rows(
                label,
                "simulation",
                data,
                contributions,
            )
        )

        dataset_path = (
            REPRESENTATIVES
            / f"{label}.npz"
        )

        np.savez_compressed(
            dataset_path,
            chrom=data["chrom"],
            pos=data["pos"],
            left_alt=synthetic["left_alt"],
            left_total=synthetic["left_total"],
            right_alt=synthetic["right_alt"],
            right_total=synthetic["right_total"],
            true_left=true_left,
            true_right=true_right,
        )

        reconstruction_rows.append(
            {
                "selection_label": label,
                "task": int(row["task"]),
                "generating_model":
                    row["generating_model"],
                "scenario": row["scenario"],
                "donor_left": left_sample,
                "donor_right": right_sample,
                "true_c_left":
                    float(row["true_c_left"]),
                "true_c_right":
                    float(row["true_c_right"]),
                "LR": float(row["LR"]),
                "site_contribution_sum":
                    contribution_sum,
                "dataset": str(dataset_path),
            }
        )

    pd.DataFrame(
        all_window_rows
    ).to_csv(
        WINDOW_EVIDENCE,
        sep="\t",
        index=False,
    )

    pd.DataFrame(
        reconstruction_rows
    ).to_csv(
        FINAL
        / "representative_reconstruction.audit.tsv",
        sep="\t",
        index=False,
    )

    audit = {
        "status": "PASS",
        "observed_sites": 1_000_000,
        "representative_simulations":
            len(representatives),
        "windows_per_dataset": 284,
        "window_rows":
            len(all_window_rows),
        "window_size": 1_000_000,
        "window_statistic":
            "fixed-global-parameter LR contribution",
        "independent_window_tests": False,
        "observed_LR":
            float(observed["LR"]),
        "observed_window_sum":
            observed_sum,
    }

    (
        FINAL
        / "window_evidence.audit.json"
    ).write_text(
        json.dumps(
            audit,
            indent=2,
        )
        + "\n"
    )

    (
        FINAL
        / "window_evidence.complete.tsv"
    ).write_text(
        "status\trepresentatives\t"
        "windows_per_dataset\n"
        f"PASS\t{len(representatives)}\t284\n"
    )

    print(
        "OBSERVED AND EMPIRICAL-SIBLING "
        "WINDOW EVIDENCE COMPLETE"
    )

    print(
        pd.DataFrame(
            reconstruction_rows
        ).to_string(
            index=False
        )
    )


def add_centroids(
    axis,
    c_left,
    c_right,
):
    emission = pilot.load_emission()

    expected_left = (
        (1.0 - c_left)
        * pilot.model.GL
        + c_left
        * pilot.model.GR
    )

    expected_right = (
        c_right
        * pilot.model.GL
        + (1.0 - c_right)
        * pilot.model.GR
    )

    x = pilot.model.map_fraction(
        expected_left,
        emission["error_rate"],
        emission["bias_log_odds"],
    )

    y = pilot.model.map_fraction(
        expected_right,
        emission["error_rate"],
        emission["bias_log_odds"],
    )

    axis.scatter(
        x[pilot.COMPATIBLE],
        y[pilot.COMPATIBLE],
        marker="s",
        s=31,
        facecolor="white",
        edgecolor=BLUE,
        linewidth=1.0,
        zorder=5,
    )

    axis.scatter(
        x[pilot.OPPOSING],
        y[pilot.OPPOSING],
        marker="D",
        s=34,
        facecolor=ORANGE,
        edgecolor="white",
        linewidth=0.7,
        zorder=6,
    )


def density(left_alt, left_total, right_alt, right_total):
    valid = (
        (left_total > 0)
        & (right_total > 0)
    )

    x = (
        left_alt[valid]
        / left_total[valid]
    )

    y = (
        right_alt[valid]
        / right_total[valid]
    )

    return np.histogram2d(
        x,
        y,
        bins=[
            EDGES,
            EDGES,
        ],
    )[0].T


def style_axis(axis):
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def make_density_page(
    pdf,
    observed,
    representatives,
):
    from matplotlib.colors import LogNorm

    selected_labels = [
        "exact_20_three_founder_median",
        "exact_20_empirical_full_sib_median",
        (
            "exact_20_empirical_full_sib_"
            "closest_to_observed"
        ),
    ]

    datasets = [
        (
            "Observed L/R",
            observed,
            float(
                observed_fit()[
                    "alternative"
                ]["c_left"]
            ),
            float(
                observed_fit()[
                    "alternative"
                ]["c_right"]
            ),
        )
    ]

    for label in selected_labels:
        row = representatives.loc[
            representatives[
                "selection_label"
            ]
            == label
        ].iloc[0]

        with np.load(
            REPRESENTATIVES
            / f"{label}.npz",
            allow_pickle=False,
        ) as source:
            dataset = {
                key: source[key]
                for key in source.files
            }

        title = (
            "Three-founder simulation"
            if row["generating_model"]
            == "three_founder"
            else (
                "Full-sibling simulation"
                if label.endswith("_median")
                else "Full sibling nearest observed LR"
            )
        )

        datasets.append(
            (
                title,
                dataset,
                float(row["true_c_left"]),
                float(row["true_c_right"]),
            )
        )

    densities = [
        density(
            dataset["left_alt"],
            dataset["left_total"],
            dataset["right_alt"],
            dataset["right_total"],
        )
        for _, dataset, _, _ in datasets
    ]

    maximum = max(
        float(values.max())
        for values in densities
    )

    figure, axes = (
        __import__(
            "matplotlib.pyplot",
            fromlist=["plt"],
        )
        .subplots(
            1,
            4,
            figsize=(14.2, 3.9),
            constrained_layout=True,
        )
    )

    artist = None

    for axis, (
        title,
        dataset,
        c_left,
        c_right,
    ), values in zip(
        axes,
        datasets,
        densities,
    ):
        artist = axis.pcolormesh(
            EDGES,
            EDGES,
            np.ma.masked_where(
                values == 0,
                values,
            ),
            cmap="Greys",
            norm=LogNorm(
                vmin=1,
                vmax=maximum,
            ),
            shading="flat",
        )

        add_centroids(
            axis,
            c_left,
            c_right,
        )

        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.set_aspect("equal")
        axis.set_xlabel(
            "Left ALT-read fraction"
        )
        axis.set_title(
            title,
            loc="left",
            weight="bold",
            fontsize=10,
        )

        style_axis(axis)

    axes[0].set_ylabel(
        "Right ALT-read fraction"
    )

    figure.colorbar(
        artist,
        ax=axes,
        label="Sites per hex-equivalent bin",
        shrink=0.78,
    )

    figure.suptitle(
        "Observed and equal-scope allele-fraction patterns",
        x=0.01,
        ha="left",
        weight="bold",
        fontsize=14,
    )

    pdf.savefig(
        figure,
        bbox_inches="tight",
    )

    summary_png = (
        FIGURES
        / "empirical_full_sib_AF_summary.png"
    )

    figure.savefig(
        summary_png,
        dpi=220,
        bbox_inches="tight",
    )

    __import__(
        "matplotlib.pyplot",
        fromlist=["plt"],
    ).close(figure)


def plot():
    import matplotlib

    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    FIGURES.mkdir(
        parents=True,
        exist_ok=True,
    )

    results = pd.read_csv(
        RESULTS,
        sep="\t",
    )

    representatives = pd.read_csv(
        REPRESENTATIVE_TABLE,
        sep="\t",
    )

    windows = pd.read_csv(
        WINDOW_EVIDENCE,
        sep="\t",
    )

    data = stage201.load_input()

    observed = {
        "left_alt": data["left_alt"],
        "left_total": data["left_total"],
        "right_alt": data["right_alt"],
        "right_total": data["right_total"],
    }

    observed_lr = float(
        observed_fit()["LR"]
    )

    pdf_path = (
        FIGURES
        / "empirical_full_sib_mixing_narrative.pdf"
    )

    with PdfPages(pdf_path) as pdf:
        make_density_page(
            pdf,
            observed,
            representatives,
        )

        figure, axes = plt.subplots(
            1,
            2,
            figsize=(11.5, 4.4),
            constrained_layout=True,
        )

        for axis, scenario in zip(
            axes,
            (
                "exact_20",
                "observed_matched",
            ),
        ):
            subsets = [
                results.loc[
                    (
                        results["scenario"]
                        == scenario
                    )
                    & (
                        results[
                            "generating_model"
                        ]
                        == model_name
                    ),
                    "LR",
                ].to_numpy()
                for model_name in (
                    "three_founder",
                    "empirical_full_sib",
                )
            ]

            violin = axis.violinplot(
                subsets,
                positions=[1, 2],
                showmedians=True,
                showextrema=False,
            )

            for body, color in zip(
                violin["bodies"],
                (
                    BLUE,
                    ORANGE,
                ),
            ):
                body.set_facecolor(color)
                body.set_edgecolor("white")
                body.set_alpha(0.72)

            violin["cmedians"].set_color(
                "black"
            )

            axis.axhline(
                observed_lr,
                color=GREEN,
                linewidth=1.8,
                linestyle="--",
                label="Observed LR",
            )

            axis.set_xticks(
                [1, 2],
                [
                    "Three founder",
                    "Empirical full sib",
                ],
            )

            axis.set_ylabel(
                "Likelihood-ratio statistic"
            )

            axis.set_title(
                scenario.replace("_", " "),
                loc="left",
                weight="bold",
            )

            axis.legend(
                frameon=False,
            )

            style_axis(axis)

        figure.suptitle(
            "Does realistic full-sibling genotype space reproduce the observed evidence?",
            x=0.01,
            ha="left",
            weight="bold",
            fontsize=14,
        )

        pdf.savefig(
            figure,
            bbox_inches="tight",
        )

        plt.close(figure)

        figure, axes = plt.subplots(
            1,
            2,
            figsize=(11.5, 4.4),
            constrained_layout=True,
        )

        colors = {
            "three_founder": BLUE,
            "empirical_full_sib": ORANGE,
        }

        for generating_model, group in results.groupby(
            "generating_model"
        ):
            axes[0].scatter(
                group["true_c_left"],
                group["estimated_c_left_alt"],
                s=12,
                alpha=0.28,
                color=colors[generating_model],
                label=generating_model.replace("_", " "),
            )

            axes[1].scatter(
                group["true_c_right"],
                group["estimated_c_right_alt"],
                s=12,
                alpha=0.28,
                color=colors[generating_model],
                label=generating_model.replace("_", " "),
            )

        for axis, side in zip(
            axes,
            ("left", "right"),
        ):
            axis.plot(
                [0.15, 0.25],
                [0.15, 0.25],
                color="black",
                linewidth=1,
                linestyle="--",
            )

            axis.set_xlim(0.15, 0.25)
            axis.set_ylim(0.15, 0.25)
            axis.set_xlabel(
                f"True {side} mixing fraction"
            )
            axis.set_ylabel(
                f"Recovered {side} mixing fraction"
            )
            axis.set_title(
                side.capitalize(),
                loc="left",
                weight="bold",
            )

            style_axis(axis)

        axes[0].legend(
            frameon=False,
        )

        figure.suptitle(
            "Equal-scope mixing-fraction recovery",
            x=0.01,
            ha="left",
            weight="bold",
            fontsize=14,
        )

        pdf.savefig(
            figure,
            bbox_inches="tight",
        )

        plt.close(figure)

        full_sib = results.loc[
            (
                results["scenario"]
                == "exact_20"
            )
            & (
                results["generating_model"]
                == "empirical_full_sib"
            )
        ].copy()

        pair_summary = (
            full_sib.groupby(
                "pair_index",
                as_index=False,
            )
            .agg(
                {
                    "genotype_contrast": "mean",
                    "IBD0_fraction": "first",
                    "IBD1_fraction": "first",
                    "IBD2_fraction": "first",
                    "LR": "mean",
                }
            )
        )

        figure, axes = plt.subplots(
            2,
            2,
            figsize=(10.5, 8.0),
            constrained_layout=True,
        )

        variables = [
            (
                "genotype_contrast",
                "Genotype contrast",
            ),
            (
                "IBD0_fraction",
                "IBD0 fraction",
            ),
            (
                "IBD1_fraction",
                "IBD1 fraction",
            ),
            (
                "IBD2_fraction",
                "IBD2 fraction",
            ),
        ]

        for axis, (
            variable,
            label,
        ) in zip(
            axes.flat,
            variables,
        ):
            x = pair_summary[
                variable
            ].to_numpy()

            y = pair_summary[
                "LR"
            ].to_numpy()

            axis.scatter(
                x,
                y,
                color=ORANGE,
                edgecolor="white",
                linewidth=0.35,
                s=38,
                alpha=0.82,
            )

            if np.unique(x).size > 1:
                coefficients = np.polyfit(
                    x,
                    y,
                    1,
                )

                ordered = np.sort(x)

                axis.plot(
                    ordered,
                    np.polyval(
                        coefficients,
                        ordered,
                    ),
                    color=BLUE,
                    linewidth=1.6,
                )

            axis.axhline(
                observed_lr,
                color=GREEN,
                linestyle="--",
                linewidth=1.3,
            )

            axis.set_xlabel(label)
            axis.set_ylabel("Mean LR")
            style_axis(axis)

        figure.suptitle(
            "How sibling genotype structure changes four-founder evidence",
            x=0.01,
            ha="left",
            weight="bold",
            fontsize=14,
        )

        pdf.savefig(
            figure,
            bbox_inches="tight",
        )

        plt.close(figure)

        manhattan_labels = [
            "observed_LR",
            "exact_20_three_founder_median",
            (
                "exact_20_empirical_full_sib_"
                "closest_to_observed"
            ),
        ]

        chromosome_lengths = (
            windows.groupby("chromosome")[
                "window_end"
            ]
            .max()
            .to_dict()
        )

        offsets = {}
        running = 0

        for chromosome in range(1, 21):
            offsets[chromosome] = running

            running += int(
                chromosome_lengths[
                    chromosome
                ]
            )

        figure, axes = plt.subplots(
            3,
            1,
            figsize=(13, 8.5),
            sharex=True,
            constrained_layout=True,
        )

        titles = [
            "Observed L/R",
            "Median three-founder simulation",
            "Full-sibling simulation nearest observed LR",
        ]

        for axis, label, title in zip(
            axes,
            manhattan_labels,
            titles,
        ):
            subset = windows.loc[
                windows["source_label"]
                == label
            ].copy()

            subset["genome_position"] = [
                offsets[int(chromosome)]
                + float(start)
                + 500_000
                for chromosome, start in zip(
                    subset["chromosome"],
                    subset["window_start"],
                )
            ]

            for chromosome, group in subset.groupby(
                "chromosome"
            ):
                axis.scatter(
                    group["genome_position"],
                    group[
                        "LR_contribution_per_site"
                    ],
                    s=15,
                    color=(
                        BLUE
                        if int(chromosome) % 2
                        else LIGHT_BLUE
                    ),
                )

            axis.axhline(
                0,
                color=GREY,
                linewidth=0.8,
            )

            axis.set_ylabel(
                "LR contribution/site"
            )

            axis.set_title(
                title,
                loc="left",
                weight="bold",
            )

            style_axis(axis)

        chromosome_centres = [
            offsets[chromosome]
            + chromosome_lengths[
                chromosome
            ] / 2
            for chromosome in range(1, 21)
        ]

        axes[-1].set_xticks(
            chromosome_centres,
            [
                str(chromosome)
                for chromosome in range(1, 21)
            ],
        )

        axes[-1].set_xlabel(
            "Autosome"
        )

        figure.suptitle(
            "Chromosome-scale distribution of founder-model evidence",
            x=0.01,
            ha="left",
            weight="bold",
            fontsize=14,
        )

        pdf.savefig(
            figure,
            bbox_inches="tight",
        )

        plt.close(figure)

    summary = pd.read_csv(
        SUMMARY,
        sep="\t",
    )

    audit = {
        "status": "PASS",
        "pages": 5,
        "simulation_rows": len(results),
        "sites_per_simulation": 1_000_000,
        "sibling_pairs": 66,
        "representative_datasets":
            len(representatives),
        "window_size": 1_000_000,
        "window_statistic":
            "fixed-global-parameter LR contribution",
        "observed_LR": observed_lr,
        "posterior_probability_assigned": False,
        "likelihood_changed": False,
        "pdf": str(pdf_path),
        "summary_png": str(
            FIGURES
            / "empirical_full_sib_AF_summary.png"
        ),
        "simulation_summary":
            summary.to_dict(
                orient="records"
            ),
    }

    (
        FIGURES
        / "empirical_full_sib_narrative.audit.json"
    ).write_text(
        json.dumps(
            audit,
            indent=2,
        )
        + "\n"
    )

    print(
        "EMPIRICAL FULL-SIBLING "
        "NARRATIVE COMPLETE"
    )

    print(
        summary.to_string(
            index=False
        )
    )

    print(f"NARRATIVE_PDF\t{pdf_path}")

    print(
        "SUMMARY_PNG\t"
        + str(
            FIGURES
            / "empirical_full_sib_AF_summary.png"
        )
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "stage",
        choices=[
            "collect",
            "evidence",
            "plot",
        ],
    )

    arguments = parser.parse_args()

    if arguments.stage == "collect":
        collect()
    elif arguments.stage == "evidence":
        evidence()
    else:
        plot()


if __name__ == "__main__":
    main()
