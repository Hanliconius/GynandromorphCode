#!/usr/bin/env python3

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


ROOT = Path(
    "analysis/genomewide_mixing_corrected_v1"
)

EXPANDED_TABLE = (
    ROOT
    / "expanded_calibration_v2"
    / "simulations.tsv.gz"
)

SIBLING_TABLE = (
    ROOT
    / "empirical_full_sib_calibration_v1"
    / "final"
    / "all_empirical_sibling_simulations.tsv.gz"
)

OBSERVED_AUDIT = (
    ROOT
    / "expanded_calibration_v2"
    / "autosomal_1Mb_evidence.audit.json"
)

OUTPUT_DIRECTORY = (
    ROOT
    / "sequential_founder_narrative"
    / "distribution_ladder"
)

OUTPUT_PDF = (
    OUTPUT_DIRECTORY
    / "sequential_founder_distribution_ladder.pdf"
)

OUTPUT_PNG = (
    OUTPUT_DIRECTORY
    / "sequential_founder_distribution_ladder.png"
)

OUTPUT_TABLE = (
    OUTPUT_DIRECTORY
    / "sequential_founder_distribution_ladder.tsv"
)

OUTPUT_AUDIT = (
    OUTPUT_DIRECTORY
    / "sequential_founder_distribution_ladder.audit.json"
)


BLUE = "#2C9BCB"
ORANGE = "#F28E2B"
BLACK = "#202020"
GREY = "#777777"


STAGES = [
    {
        "order": 1,
        "label": (
            "1  Broad parameter space\n"
            "mixing 0-45%"
        ),
        "source": "expanded",
        "subset": "global",
        "three": "three",
        "four": "four",
    },
    {
        "order": 2,
        "label": (
            "2  Focus near observed mixing\n"
            "14-26%"
        ),
        "source": "expanded",
        "subset": "focused",
        "three": "three",
        "four": "four",
    },
    {
        "order": 3,
        "label": (
            "3  Exact 20% mixing\n"
            "all donor pairs"
        ),
        "source": "expanded",
        "subset": "repeat_20",
        "three": "three",
        "four": "four",
    },
    {
        "order": 4,
        "label": (
            "4  Add empirical full-sib structure\n"
            "exact 20% mixing"
        ),
        "source": "sibling",
        "subset": "exact_20",
        "three": "three_founder",
        "four": "empirical_full_sib",
    },
    {
        "order": 5,
        "label": (
            "5  Full siblings + fitted mixing\n"
            "19.7% Left / 20.24% Right"
        ),
        "source": "sibling",
        "subset": "observed_matched",
        "three": "three_founder",
        "four": "empirical_full_sib",
    },
]


def require_file(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty input: {path}"
        )


def empirical_statistics(values, observed):
    values = np.asarray(
        values,
        dtype=float,
    )

    values = values[
        np.isfinite(values)
    ]

    if values.size == 0:
        raise RuntimeError(
            "Empty simulation distribution"
        )

    less = int(
        np.sum(values < observed)
    )

    equal = int(
        np.sum(
            np.isclose(
                values,
                observed,
                rtol=1e-12,
                atol=1e-12,
            )
        )
    )

    q = (
        less
        + 0.5 * equal
        + 0.5
    ) / (
        values.size + 1.0
    )

    compatibility = 2.0 * min(
        q,
        1.0 - q,
    )

    quantiles = np.quantile(
        values,
        [
            0.025,
            0.25,
            0.50,
            0.75,
            0.975,
        ],
    )

    return {
        "n": int(values.size),
        "minimum": float(np.min(values)),
        "q025": float(quantiles[0]),
        "q25": float(quantiles[1]),
        "median": float(quantiles[2]),
        "q75": float(quantiles[3]),
        "q975": float(quantiles[4]),
        "maximum": float(np.max(values)),
        "observed_percentile": float(q),
        "compatibility": float(
            compatibility
        ),
    }


def load_inputs():
    require_file(EXPANDED_TABLE)
    require_file(SIBLING_TABLE)
    require_file(OBSERVED_AUDIT)

    expanded = pd.read_csv(
        EXPANDED_TABLE,
        sep="\t",
        compression="gzip",
    )

    siblings = pd.read_csv(
        SIBLING_TABLE,
        sep="\t",
        compression="gzip",
    )

    expanded_required = {
        "family",
        "generating_model",
        "LR",
        "sites",
    }

    sibling_required = {
        "scenario",
        "generating_model",
        "LR",
        "sites",
    }

    missing_expanded = (
        expanded_required
        - set(expanded.columns)
    )

    missing_siblings = (
        sibling_required
        - set(siblings.columns)
    )

    if missing_expanded:
        raise RuntimeError(
            "Expanded table missing columns: "
            f"{sorted(missing_expanded)}"
        )

    if missing_siblings:
        raise RuntimeError(
            "Sibling table missing columns: "
            f"{sorted(missing_siblings)}"
        )

    observed_audit = json.loads(
        OBSERVED_AUDIT.read_text()
    )

    observed_lr = float(
        observed_audit["global_LR"]
    )

    observed_sites = int(
        observed_audit["sites"]
    )

    observed_metric = (
        1000.0
        * observed_lr
        / observed_sites
    )

    return (
        expanded,
        siblings,
        observed_lr,
        observed_sites,
        observed_metric,
    )


def select_distribution(
    stage,
    model_role,
    expanded,
    siblings,
):
    model_name = stage[model_role]

    if stage["source"] == "expanded":
        subset = expanded.loc[
            expanded["family"].eq(
                stage["subset"]
            )
            & expanded[
                "generating_model"
            ].eq(model_name)
        ].copy()

    elif stage["source"] == "sibling":
        subset = siblings.loc[
            siblings["scenario"].eq(
                stage["subset"]
            )
            & siblings[
                "generating_model"
            ].eq(model_name)
        ].copy()

    else:
        raise RuntimeError(
            "Unknown stage source: "
            f"{stage['source']}"
        )

    lr = pd.to_numeric(
        subset["LR"],
        errors="coerce",
    )

    sites = pd.to_numeric(
        subset["sites"],
        errors="coerce",
    )

    valid = (
        lr.notna()
        & sites.notna()
        & np.isfinite(lr)
        & np.isfinite(sites)
        & sites.gt(0)
    )

    lr = lr.loc[valid].to_numpy(
        dtype=float
    )

    sites = sites.loc[valid].to_numpy(
        dtype=float
    )

    if lr.size == 0:
        raise RuntimeError(
            f"No valid rows for stage "
            f"{stage['order']} / {model_role}"
        )

    if np.any(lr < -1e-6):
        raise RuntimeError(
            f"Materially negative LR in stage "
            f"{stage['order']} / {model_role}"
        )

    lr = np.maximum(
        lr,
        0.0,
    )

    metric = (
        1000.0
        * lr
        / sites
    )

    return metric, sorted(
        np.unique(sites.astype(int)).tolist()
    )


def assemble_results():
    (
        expanded,
        siblings,
        observed_lr,
        observed_sites,
        observed_metric,
    ) = load_inputs()

    distributions = {}
    records = []

    for stage in STAGES:
        for role in [
            "three",
            "four",
        ]:
            values, site_counts = (
                select_distribution(
                    stage,
                    role,
                    expanded,
                    siblings,
                )
            )

            statistics = empirical_statistics(
                values,
                observed_metric,
            )

            key = (
                int(stage["order"]),
                role,
            )

            distributions[key] = values

            record = {
                "stage": int(stage["order"]),
                "stage_label": stage["label"],
                "source": stage["source"],
                "subset": stage["subset"],
                "model_role": role,
                "generating_model": stage[role],
                "site_counts": ",".join(
                    str(value)
                    for value in site_counts
                ),
                "metric": "LR_per_1000_sites",
                "observed_metric": (
                    observed_metric
                ),
                **statistics,
            }

            records.append(record)

    results = pd.DataFrame(records)

    return (
        results,
        distributions,
        observed_lr,
        observed_sites,
        observed_metric,
    )


def axis_limits(
    distributions,
    observed_metric,
):
    three_values = np.concatenate(
        [
            distributions[
                (stage["order"], "three")
            ]
            for stage in STAGES
        ]
    )

    four_values = np.concatenate(
        [
            distributions[
                (stage["order"], "four")
            ]
            for stage in STAGES
        ]
    )

    left_high = float(
        np.quantile(
            three_values,
            0.995,
        )
    )

    left_high = max(
        left_high * 1.18,
        float(np.max(three_values)) * 1.03,
        1e-8,
    )

    right_low = min(
        float(
            np.quantile(
                four_values,
                0.005,
            )
        ),
        observed_metric,
    )

    right_low = max(
        right_low * 0.90,
        left_high * 2.5,
    )

    right_high = max(
        float(
            np.quantile(
                four_values,
                0.995,
            )
        ),
        observed_metric,
    )

    right_high = max(
        right_high * 1.08,
        right_low * 1.2,
    )

    if right_low >= observed_metric:
        right_low = max(
            left_high * 2.5,
            observed_metric * 0.80,
        )

    if right_low <= left_high:
        raise RuntimeError(
            "No clean axis break exists: "
            f"left_high={left_high}, "
            f"right_low={right_low}"
        )

    return (
        left_high,
        right_low,
        right_high,
    )


def draw_distribution(
    axis,
    values,
    y,
    colour,
    rng,
):
    values = np.asarray(
        values,
        dtype=float,
    )

    if values.size > 180:
        selected = rng.choice(
            values,
            size=180,
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
        y + jitter,
        s=9,
        color=colour,
        alpha=0.16,
        linewidth=0,
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
        y,
        q025,
        q975,
        color=colour,
        linewidth=2.0,
        alpha=0.80,
        zorder=2,
    )

    axis.hlines(
        y,
        q25,
        q75,
        color=colour,
        linewidth=7.0,
        alpha=0.92,
        zorder=3,
    )

    axis.scatter(
        [median],
        [y],
        s=64,
        color=colour,
        edgecolor="white",
        linewidth=1.1,
        zorder=4,
    )


def percentile_label(value):
    percentile = 100.0 * value

    if percentile < 0.1:
        return "<0.1"

    if percentile > 99.9:
        return ">99.9"

    if percentile < 10.0:
        return f"{percentile:.1f}"

    if percentile > 90.0:
        return f"{percentile:.1f}"

    return f"{percentile:.0f}"


def make_plot(
    results,
    distributions,
    observed_metric,
    limits,
):
    left_high, right_low, right_high = (
        limits
    )

    figure = plt.figure(
        figsize=(13.2, 7.2)
    )

    grid = figure.add_gridspec(
        nrows=1,
        ncols=2,
        width_ratios=[
            1.15,
            3.85,
        ],
        left=0.30,
        right=0.88,
        top=0.78,
        bottom=0.17,
        wspace=0.045,
    )

    left_axis = figure.add_subplot(
        grid[0, 0]
    )

    right_axis = figure.add_subplot(
        grid[0, 1],
        sharey=left_axis,
    )

    rng = np.random.default_rng(
        2140917
    )

    stage_y = {
        stage["order"]:
            len(STAGES) - stage["order"]
        for stage in STAGES
    }

    offset = 0.13

    for stage in STAGES:
        centre = stage_y[
            stage["order"]
        ]

        three_y = centre + offset
        four_y = centre - offset

        draw_distribution(
            left_axis,
            distributions[
                (stage["order"], "three")
            ],
            three_y,
            BLUE,
            rng,
        )

        draw_distribution(
            right_axis,
            distributions[
                (stage["order"], "four")
            ],
            four_y,
            ORANGE,
            rng,
        )

        left_axis.text(
            left_high * 0.97,
            three_y,
            "3F",
            ha="right",
            va="center",
            color=BLUE,
            fontsize=8.5,
            fontweight="bold",
        )

        right_axis.text(
            right_low
            + 0.015
            * (right_high - right_low),
            four_y,
            "4F",
            ha="left",
            va="center",
            color=ORANGE,
            fontsize=8.5,
            fontweight="bold",
        )

        three_row = results.loc[
            results["stage"].eq(
                stage["order"]
            )
            & results["model_role"].eq(
                "three"
            )
        ].iloc[0]

        four_row = results.loc[
            results["stage"].eq(
                stage["order"]
            )
            & results["model_role"].eq(
                "four"
            )
        ].iloc[0]

        annotation = (
            "Observed percentile: "
            f"3F {percentile_label(three_row['observed_percentile'])}%"
            "  |  "
            f"4F {percentile_label(four_row['observed_percentile'])}%"
        )

        right_axis.text(
            1.025,
            centre,
            annotation,
            transform=(
                right_axis
                .get_yaxis_transform()
            ),
            ha="left",
            va="center",
            fontsize=8.4,
            color="#4A4A4A",
            clip_on=False,
        )

    right_axis.axvline(
        observed_metric,
        color=BLACK,
        linewidth=2.1,
        zorder=0,
    )

    right_axis.text(
        observed_metric,
        len(STAGES) - 0.40,
        (
            "Observed\n"
            f"{observed_metric:.3g}"
        ),
        ha="center",
        va="bottom",
        color=BLACK,
        fontsize=9,
        fontweight="bold",
    )

    left_axis.annotate(
        "Observed lies beyond this zoom",
        xy=(
            left_high,
            len(STAGES) - 0.55,
        ),
        xytext=(
            left_high * 0.50,
            len(STAGES) - 0.55,
        ),
        arrowprops={
            "arrowstyle": "-|>",
            "color": BLACK,
            "linewidth": 1.1,
        },
        ha="center",
        va="center",
        fontsize=8.2,
        color=BLACK,
    )

    labels = [
        stage["label"]
        for stage in STAGES
    ]

    y_ticks = [
        stage_y[stage["order"]]
        for stage in STAGES
    ]

    left_axis.set_yticks(
        y_ticks
    )

    left_axis.set_yticklabels(
        labels,
        fontsize=10.2,
    )

    left_axis.set_xlim(
        0.0,
        left_high,
    )

    right_axis.set_xlim(
        right_low,
        right_high,
    )

    left_axis.set_ylim(
        -0.55,
        len(STAGES) - 0.35,
    )

    left_axis.set_xlabel(
        "Near-zero region",
        fontsize=10,
    )

    right_axis.set_xlabel(
        "Likelihood-ratio evidence per 1,000 sites",
        fontsize=11,
    )

    left_axis.ticklabel_format(
        axis="x",
        style="sci",
        scilimits=(-2, 2),
    )

    for axis in [
        left_axis,
        right_axis,
    ]:
        axis.grid(
            axis="x",
            color="#E8E8E8",
            linewidth=0.8,
            zorder=0,
        )

        axis.spines["top"].set_visible(
            False
        )

        axis.spines["right"].set_visible(
            False
        )

        axis.tick_params(
            axis="y",
            length=0,
        )

    right_axis.spines[
        "left"
    ].set_visible(False)

    right_axis.tick_params(
        axis="y",
        left=False,
        labelleft=False,
    )

    left_axis.spines[
        "right"
    ].set_visible(False)

    break_size = 0.012

    break_arguments = {
        "color": BLACK,
        "clip_on": False,
        "linewidth": 1.0,
    }

    left_axis.plot(
        (
            1.0 - break_size,
            1.0 + break_size,
        ),
        (
            -break_size,
            break_size,
        ),
        transform=left_axis.transAxes,
        **break_arguments,
    )

    left_axis.plot(
        (
            1.0 - break_size,
            1.0 + break_size,
        ),
        (
            1.0 - break_size,
            1.0 + break_size,
        ),
        transform=left_axis.transAxes,
        **break_arguments,
    )

    right_axis.plot(
        (
            -break_size,
            break_size,
        ),
        (
            -break_size,
            break_size,
        ),
        transform=right_axis.transAxes,
        **break_arguments,
    )

    right_axis.plot(
        (
            -break_size,
            break_size,
        ),
        (
            1.0 - break_size,
            1.0 + break_size,
        ),
        transform=right_axis.transAxes,
        **break_arguments,
    )

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=BLUE,
            marker="o",
            linewidth=5,
            label="Three-founder simulations",
        ),
        Line2D(
            [0],
            [0],
            color=ORANGE,
            marker="o",
            linewidth=5,
            label=(
                "Four-founder / empirical "
                "full-sib simulations"
            ),
        ),
        Line2D(
            [0],
            [0],
            color=BLACK,
            linewidth=2,
            label="Observed LR",
        ),
    ]

    figure.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(
            0.30,
            0.835,
        ),
        frameon=False,
        ncol=3,
        fontsize=9.5,
    )

    figure.text(
        0.055,
        0.945,
        (
            "Sequential calibration retains the "
            "three- versus four-founder comparison"
        ),
        fontsize=18,
        fontweight="bold",
        ha="left",
        va="top",
        color=BLACK,
    )

    figure.text(
        0.055,
        0.895,
        (
            "Points are simulation replicates; thin lines show "
            "95% intervals, thick lines show interquartile "
            "ranges, and circles show medians."
        ),
        fontsize=10.5,
        ha="left",
        va="top",
        color="#555555",
    )

    figure.text(
        0.055,
        0.055,
        (
            "The broken axis magnifies the near-zero "
            "three-founder distributions. LR is divided by "
            "the number of analysed sites before scaling to "
            "1,000 sites. Percentiles report the observed "
            "position within each simulated distribution."
        ),
        fontsize=9.0,
        ha="left",
        va="bottom",
        color="#666666",
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


def main(validate_only=False):
    (
        results,
        distributions,
        observed_lr,
        observed_sites,
        observed_metric,
    ) = assemble_results()

    limits = axis_limits(
        distributions,
        observed_metric,
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
        "axis_limits="
        f"{limits}"
    )

    print(
        results[
            [
                "stage",
                "subset",
                "model_role",
                "generating_model",
                "n",
                "median",
                "q025",
                "q975",
                "observed_percentile",
                "compatibility",
            ]
        ].to_string(index=False)
    )

    if validate_only:
        print("VALIDATION_PASS")
        return

    OUTPUT_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    results.to_csv(
        OUTPUT_TABLE,
        sep="\t",
        index=False,
    )

    make_plot(
        results,
        distributions,
        observed_metric,
        limits,
    )

    audit = {
        "observed_LR": observed_lr,
        "observed_sites": observed_sites,
        "observed_LR_per_1000_sites": (
            observed_metric
        ),
        "stages": len(STAGES),
        "paired_three_vs_four_at_every_stage": (
            True
        ),
        "distribution_display": (
            "replicate points, 95% interval, IQR, median"
        ),
        "axis": (
            "broken linear axis with separate near-zero zoom"
        ),
        "site_scope_normalization": (
            "LR divided by sites and multiplied by 1000"
        ),
        "posterior_probability_assigned": False,
        "axis_limits": {
            "near_zero_maximum": limits[0],
            "right_minimum": limits[1],
            "right_maximum": limits[2],
        },
        "pdf": str(OUTPUT_PDF),
        "png": str(OUTPUT_PNG),
        "table": str(OUTPUT_TABLE),
    }

    OUTPUT_AUDIT.write_text(
        json.dumps(
            audit,
            indent=2,
        )
        + "\n"
    )

    print("STAGE_214_COMPLETE")

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
