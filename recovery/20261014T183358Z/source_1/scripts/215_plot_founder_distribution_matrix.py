#!/usr/bin/env python3

import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(
    "analysis/genomewide_mixing_corrected_v1"
)

STAGE214_SCRIPT = Path(
    "scripts/214_plot_founder_distribution_ladder.py"
)

OUTPUT_DIRECTORY = (
    ROOT
    / "sequential_founder_narrative"
    / "distribution_matrix"
)

OUTPUT_PDF = (
    OUTPUT_DIRECTORY
    / "sequential_founder_distribution_matrix.pdf"
)

OUTPUT_PNG = (
    OUTPUT_DIRECTORY
    / "sequential_founder_distribution_matrix.png"
)

OUTPUT_AUDIT = (
    OUTPUT_DIRECTORY
    / "sequential_founder_distribution_matrix.audit.json"
)


BLUE = "#2C9BCB"
ORANGE = "#F28E2B"
BLACK = "#202020"

BLUE_BACKGROUND = "#EDF7FB"
ORANGE_BACKGROUND = "#FFF4E8"


def load_stage214():
    if (
        not STAGE214_SCRIPT.is_file()
        or STAGE214_SCRIPT.stat().st_size == 0
    ):
        raise RuntimeError(
            f"Missing Stage 214 script: "
            f"{STAGE214_SCRIPT}"
        )

    specification = (
        importlib.util.spec_from_file_location(
            "stage214_distribution_source",
            STAGE214_SCRIPT,
        )
    )

    module = (
        importlib.util.module_from_spec(
            specification
        )
    )

    specification.loader.exec_module(
        module
    )

    return module


def percentile_statement(q):
    percentage = 100.0 * q

    if q >= 0.5:
        if percentage > 99.9:
            value = ">99.9%"
        elif percentage >= 99.0:
            value = f"{percentage:.1f}%"
        else:
            value = f"{percentage:.0f}%"

        return (
            f"Observed above {value} "
            "of simulations"
        )

    above = 100.0 * (1.0 - q)

    if above > 99.9:
        value = ">99.9%"
    elif above >= 99.0:
        value = f"{above:.1f}%"
    else:
        value = f"{above:.0f}%"

    return (
        f"Observed below {value} "
        "of simulations"
    )


def distribution_limits(
    values,
    observed,
    show_observed,
):
    values = np.asarray(
        values,
        dtype=float,
    )

    q005, q995 = np.quantile(
        values,
        [
            0.005,
            0.995,
        ],
    )

    if show_observed:
        lower = min(
            q005,
            observed,
        )

        upper = max(
            q995,
            observed,
        )

        span = upper - lower

        if span <= 0:
            span = max(
                abs(upper),
                1e-8,
            )

        lower = max(
            0.0,
            lower - 0.06 * span,
        )

        upper = (
            upper
            + 0.10 * span
        )

    else:
        lower = max(
            0.0,
            q005,
        )

        upper = max(
            q995,
            float(np.max(values)),
            1e-10,
        )

        span = upper - lower

        if span <= 0:
            span = max(
                abs(upper),
                1e-8,
            )

        lower = max(
            0.0,
            lower - 0.04 * span,
        )

        upper = (
            upper
            + 0.10 * span
        )

    if upper <= lower:
        upper = lower + 1e-8

    return lower, upper


def draw_distribution_cell(
    axis,
    values,
    observed,
    colour,
    background,
    q,
    rng,
    show_observed,
    model_label,
):
    values = np.asarray(
        values,
        dtype=float,
    )

    values = values[
        np.isfinite(values)
    ]

    if values.size == 0:
        raise RuntimeError(
            "Cannot draw an empty distribution"
        )

    axis.set_facecolor(
        background
    )

    if values.size > 220:
        selected = rng.choice(
            values,
            size=220,
            replace=False,
        )
    else:
        selected = values

    jitter = rng.uniform(
        -0.055,
        0.055,
        size=selected.size,
    )

    axis.scatter(
        selected,
        jitter,
        s=11,
        color=colour,
        alpha=0.18,
        edgecolor="none",
        rasterized=True,
        zorder=1,
    )

    q025, q25, median, q75, q975 = (
        np.quantile(
            values,
            [
                0.025,
                0.25,
                0.50,
                0.75,
                0.975,
            ],
        )
    )

    axis.hlines(
        0.0,
        q025,
        q975,
        color=colour,
        linewidth=2.2,
        alpha=0.82,
        zorder=2,
    )

    axis.hlines(
        0.0,
        q25,
        q75,
        color=colour,
        linewidth=8.0,
        alpha=0.95,
        zorder=3,
    )

    axis.scatter(
        [median],
        [0.0],
        s=72,
        color=colour,
        edgecolor="white",
        linewidth=1.2,
        zorder=4,
    )

    lower, upper = distribution_limits(
        values,
        observed,
        show_observed,
    )

    axis.set_xlim(
        lower,
        upper,
    )

    axis.set_ylim(
        -0.22,
        0.34,
    )

    if show_observed:
        axis.axvline(
            observed,
            color=BLACK,
            linewidth=2.0,
            zorder=5,
        )

        axis.text(
            observed,
            0.205,
            "Observed",
            ha="center",
            va="bottom",
            fontsize=8.2,
            fontweight="bold",
            color=BLACK,
        )

    else:
        axis.annotate(
            f"Observed {observed:.2f} is off-scale",
            xy=(
                0.985,
                0.72,
            ),
            xycoords="axes fraction",
            xytext=(
                0.48,
                0.72,
            ),
            textcoords="axes fraction",
            arrowprops={
                "arrowstyle": "-|>",
                "color": BLACK,
                "linewidth": 1.2,
            },
            ha="center",
            va="center",
            fontsize=8.0,
            color=BLACK,
        )

    axis.text(
        0.018,
        0.86,
        model_label,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8.7,
        color=colour,
        fontweight="bold",
    )

    axis.text(
        0.985,
        0.86,
        percentile_statement(q),
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=8.1,
        color="#444444",
    )

    axis.text(
        0.018,
        0.10,
        (
            f"median {median:.3g}; "
            f"95% interval "
            f"{q025:.3g}-{q975:.3g}"
        ),
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.6,
        color="#555555",
    )

    axis.set_yticks([])

    axis.tick_params(
        axis="x",
        labelsize=8,
        length=3,
    )

    axis.grid(
        axis="x",
        color="white",
        linewidth=1.0,
        alpha=0.95,
        zorder=0,
    )

    for spine in [
        "top",
        "right",
        "left",
    ]:
        axis.spines[
            spine
        ].set_visible(False)

    axis.spines[
        "bottom"
    ].set_color("#888888")


def main(validate_only=False):
    stage214 = load_stage214()

    (
        results,
        distributions,
        observed_lr,
        observed_sites,
        observed_metric,
    ) = stage214.assemble_results()

    expected_keys = {
        (
            stage["order"],
            role,
        )
        for stage in stage214.STAGES
        for role in [
            "three",
            "four",
        ]
    }

    missing = (
        expected_keys
        - set(distributions)
    )

    if missing:
        raise RuntimeError(
            "Missing distributions: "
            f"{sorted(missing)}"
        )

    print(
        f"observed_LR={observed_lr}"
    )

    print(
        f"observed_sites={observed_sites}"
    )

    print(
        "observed_LR_per_1000_sites="
        f"{observed_metric}"
    )

    print(
        results[
            [
                "stage",
                "model_role",
                "generating_model",
                "n",
                "median",
                "q025",
                "q975",
                "observed_percentile",
            ]
        ].to_string(index=False)
    )

    if validate_only:
        print("VALIDATION_PASS")
        return

    figure, axes = plt.subplots(
        nrows=len(stage214.STAGES),
        ncols=2,
        figsize=(13.2, 9.0),
        constrained_layout=False,
    )

    figure.subplots_adjust(
        left=0.285,
        right=0.975,
        top=0.83,
        bottom=0.105,
        hspace=0.42,
        wspace=0.16,
    )

    rng = np.random.default_rng(
        2150917
    )

    for row_index, stage in enumerate(
        stage214.STAGES
    ):
        three_axis = axes[
            row_index,
            0,
        ]

        four_axis = axes[
            row_index,
            1,
        ]

        three_values = distributions[
            (
                stage["order"],
                "three",
            )
        ]

        four_values = distributions[
            (
                stage["order"],
                "four",
            )
        ]

        three_result = results.loc[
            results["stage"].eq(
                stage["order"]
            )
            & results[
                "model_role"
            ].eq("three")
        ].iloc[0]

        four_result = results.loc[
            results["stage"].eq(
                stage["order"]
            )
            & results[
                "model_role"
            ].eq("four")
        ].iloc[0]

        draw_distribution_cell(
            three_axis,
            three_values,
            observed_metric,
            BLUE,
            BLUE_BACKGROUND,
            float(
                three_result[
                    "observed_percentile"
                ]
            ),
            rng,
            False,
            "Three-founder",
        )

        four_name = (
            "Four-founder"
            if stage["source"]
            == "expanded"
            else
            "Empirical full-sibling"
        )

        draw_distribution_cell(
            four_axis,
            four_values,
            observed_metric,
            ORANGE,
            ORANGE_BACKGROUND,
            float(
                four_result[
                    "observed_percentile"
                ]
            ),
            rng,
            True,
            four_name,
        )

        three_axis.set_ylabel(
            stage["label"],
            rotation=0,
            ha="right",
            va="center",
            fontsize=10.2,
            labelpad=18,
        )

        three_axis.yaxis.set_label_coords(
            -0.06,
            0.50,
        )

    axes[-1, 0].set_xlabel(
        "LR evidence per 1,000 sites\n"
        "(independent three-founder zoom)",
        fontsize=10,
    )

    axes[-1, 1].set_xlabel(
        "LR evidence per 1,000 sites\n"
        "(independent four-founder/sibling scale)",
        fontsize=10,
    )

    figure.text(
        0.055,
        0.955,
        (
            "How successive calibration layers change "
            "three- and four-founder evidence"
        ),
        fontsize=18,
        fontweight="bold",
        ha="left",
        va="top",
        color=BLACK,
    )

    figure.text(
        0.055,
        0.914,
        (
            "Each cell has its own labelled x-scale so that "
            "both narrow and broad simulation distributions "
            "remain visible."
        ),
        fontsize=10.4,
        ha="left",
        va="top",
        color="#555555",
    )

    figure.text(
        0.46,
        0.865,
        "THREE-FOUNDER SIMULATIONS",
        fontsize=11,
        fontweight="bold",
        color=BLUE,
        ha="center",
        va="center",
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": BLUE_BACKGROUND,
            "edgecolor": BLUE,
            "linewidth": 1.2,
        },
    )

    figure.text(
        0.785,
        0.865,
        (
            "FOUR-FOUNDER / "
            "EMPIRICAL FULL-SIB SIMULATIONS"
        ),
        fontsize=11,
        fontweight="bold",
        color=ORANGE,
        ha="center",
        va="center",
        bbox={
            "boxstyle": "round,pad=0.35",
            "facecolor": ORANGE_BACKGROUND,
            "edgecolor": ORANGE,
            "linewidth": 1.2,
        },
    )

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=BLUE,
            marker="o",
            linewidth=6,
            label=(
                "Simulation median and "
                "interquartile range"
            ),
        ),
        Line2D(
            [0],
            [0],
            color=BLACK,
            linewidth=2,
            label=(
                "Observed LR when on-scale"
            ),
        ),
    ]

    figure.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(
            0.63,
            0.027,
        ),
        frameon=False,
        ncol=2,
        fontsize=9,
    )

    figure.text(
        0.055,
        0.020,
        (
            "Faint points are simulation replicates; thin lines "
            "are 95% intervals; thick lines are interquartile "
            "ranges; circles are medians. Scales differ among "
            "cells and must be read from their numeric axes."
        ),
        fontsize=8.6,
        color="#666666",
        ha="left",
        va="bottom",
    )

    OUTPUT_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure.savefig(
        OUTPUT_PDF,
        bbox_inches="tight",
    )

    figure.savefig(
        OUTPUT_PNG,
        dpi=240,
        bbox_inches="tight",
    )

    plt.close(figure)

    audit = {
        "observed_LR": observed_lr,
        "observed_sites": observed_sites,
        "observed_LR_per_1000_sites": (
            observed_metric
        ),
        "rows": len(stage214.STAGES),
        "columns": 2,
        "left_column": (
            "three-founder simulations with "
            "independent local scales"
        ),
        "right_column": (
            "four-founder or empirical full-sib "
            "simulations with independent local scales"
        ),
        "observed_display": (
            "black vertical line when on-scale; "
            "off-scale arrow otherwise"
        ),
        "posterior_probability_assigned": False,
        "pdf": str(OUTPUT_PDF),
        "png": str(OUTPUT_PNG),
    }

    OUTPUT_AUDIT.write_text(
        json.dumps(
            audit,
            indent=2,
        )
        + "\n"
    )

    print("STAGE_215_COMPLETE")

    print(
        "PDF\t"
        f"{OUTPUT_PDF.resolve()}"
    )

    print(
        "PNG\t"
        f"{OUTPUT_PNG.resolve()}"
    )

    print(
        "AUDIT\t"
        f"{OUTPUT_AUDIT.resolve()}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--validate-only",
        action="store_true",
    )

    arguments = parser.parse_args()

    main(
        validate_only=arguments.validate_only
    )
