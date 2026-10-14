#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.special import logsumexp


PROJECT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT / "analysis/autosomal_three_vs_four_founder"
OUTDIR = ANALYSIS / "narrative_figures"
MODEL159 = PROJECT / "scripts/159_fit_autosomal_three_vs_four_founder.py"
MODEL160 = PROJECT / "scripts/160_simulate_calibrate_autosomal_founder_models.py"

MATCHED_FILE = ANALYSIS / "matched_autosomal_focal_counts.tsv.gz"
OBSERVED_FILE = ANALYSIS / "observed_nested_founder_fit.tsv"
SIM_FILE = ANALYSIS / "founder_simulation_replicates.tsv.gz"
MIXING_FILE = ANALYSIS / "mixing_calibration.tsv"
BOOTSTRAP_FILE = ANALYSIS / "observed_10kb_block_bootstrap.tsv.gz"
FINAL_FILE = ANALYSIS / "calibrated_founder_and_mixing_summary.tsv"
DONOR_FILE = ANALYSIS / "simulation_all_47_controls.tsv"

SITE_OUT = ANALYSIS / "observed_nested_LR_site_contributions.tsv.gz"
BLOCK_OUT = ANALYSIS / "observed_nested_LR_10kb_block_contributions.tsv"
NARRATIVE_PDF = OUTDIR / "three_vs_four_founder_mixing_narrative.pdf"
SUMMARY_PDF = OUTDIR / "three_vs_four_founder_mixing_summary_figure.pdf"
SUMMARY_PNG = OUTDIR / "three_vs_four_founder_mixing_summary_figure.png"
AUDIT_OUT = OUTDIR / "three_vs_four_founder_mixing_narrative.audit.tsv"

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
BLACK = "#222222"
GREY = "#777777"
OPTIX_START = 1_088_602
OPTIX_END = 1_097_846
KEYS = ["chrom", "pos", "ref", "alt"]


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def numeric_summary(path):
    table = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    if set(table.columns) != {"statistic", "value"}:
        raise RuntimeError(
            f"Unexpected summary columns in {path}"
        )
    return dict(zip(table["statistic"], table["value"]))


def as_float(summary, key):
    if key not in summary:
        raise RuntimeError(
            f"Missing summary statistic: {key}"
        )
    return float(summary[key])


def set_style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "axes.linewidth": 0.8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
    })


def clean_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(
        axis="y",
        color="#ECECEC",
        linewidth=0.6,
        zorder=0,
    )


def add_page_heading(
    figure,
    step,
    title,
    subtitle,
):
    figure.text(
        0.055,
        0.965,
        f"STEP {step}",
        color=GREY,
        fontsize=10,
        weight="bold",
        va="top",
    )
    figure.text(
        0.055,
        0.925,
        title,
        color=BLACK,
        fontsize=18,
        weight="bold",
        va="top",
    )
    figure.text(
        0.055,
        0.885,
        subtitle,
        color="#4A4A4A",
        fontsize=9.5,
        va="top",
    )


def model_schematic(ax, compact=False):
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    box = dict(
        boxstyle="round,pad=0.35",
        linewidth=1.2,
    )

    if compact:
        y3 = 0.70
        y4 = 0.27

        ax.text(
            0.03,
            y3,
            "Three founders",
            weight="bold",
            color=BLUE,
            va="center",
        )
        ax.text(
            0.34,
            y3,
            "Left: A + B",
            bbox={
                **box,
                "ec": BLUE,
                "fc": "#EAF4FA",
            },
            va="center",
        )
        ax.text(
            0.70,
            y3,
            "Right: A + C",
            bbox={
                **box,
                "ec": BLUE,
                "fc": "#EAF4FA",
            },
            va="center",
        )
        ax.text(
            0.50,
            y3 - 0.17,
            "Shared A; opposing homozygotes impossible",
            color=GREEN,
            ha="center",
            fontsize=8,
        )

        ax.text(
            0.03,
            y4,
            "Four founders",
            weight="bold",
            color=ORANGE,
            va="center",
        )
        ax.text(
            0.34,
            y4,
            "Left: A + B",
            bbox={
                **box,
                "ec": ORANGE,
                "fc": "#FCEFE8",
            },
            va="center",
        )
        ax.text(
            0.70,
            y4,
            "Right: C + D",
            bbox={
                **box,
                "ec": ORANGE,
                "fc": "#FCEFE8",
            },
            va="center",
        )
        ax.text(
            0.50,
            y4 - 0.17,
            "No required shared gamete; all genotype pairs possible",
            color=ORANGE,
            ha="center",
            fontsize=8,
        )
        return

    ax.text(
        0.25,
        0.92,
        "Three-founder hypothesis",
        ha="center",
        weight="bold",
        fontsize=15,
        color=BLUE,
    )
    ax.text(
        0.75,
        0.92,
        "Four-founder hypothesis",
        ha="center",
        weight="bold",
        fontsize=15,
        color=ORANGE,
    )
    ax.text(
        0.25,
        0.78,
        "Left = A + B",
        ha="center",
        bbox={
            **box,
            "ec": BLUE,
            "fc": "#EAF4FA",
        },
        fontsize=12,
    )
    ax.text(
        0.25,
        0.61,
        "Right = A + C",
        ha="center",
        bbox={
            **box,
            "ec": BLUE,
            "fc": "#EAF4FA",
        },
        fontsize=12,
    )
    ax.annotate(
        "shared gamete A",
        xy=(0.25, 0.69),
        xytext=(0.43, 0.69),
        arrowprops=dict(
            arrowstyle="-[,widthB=2.6",
            color=GREEN,
            lw=1.5,
        ),
        color=GREEN,
        va="center",
        fontsize=10,
    )
    ax.text(
        0.25,
        0.43,
        "0/0 | 1/1 and 1/1 | 0/0\nare impossible",
        ha="center",
        color=BLUE,
        fontsize=11,
        weight="bold",
    )
    ax.text(
        0.75,
        0.78,
        "Left = A + B",
        ha="center",
        bbox={
            **box,
            "ec": ORANGE,
            "fc": "#FCEFE8",
        },
        fontsize=12,
    )
    ax.text(
        0.75,
        0.61,
        "Right = C + D",
        ha="center",
        bbox={
            **box,
            "ec": ORANGE,
            "fc": "#FCEFE8",
        },
        fontsize=12,
    )
    ax.text(
        0.75,
        0.43,
        "No obligatory shared gamete\n"
        "All nine genotype pairs are possible",
        ha="center",
        color=ORANGE,
        fontsize=11,
        weight="bold",
    )
    ax.text(
        0.50,
        0.23,
        "Sequencing mixture is fitted separately from founder state",
        ha="center",
        weight="bold",
        fontsize=12,
        color=BLACK,
    )
    ax.text(
        0.50,
        0.12,
        "Observed Left = (1 - cL) x true Left + "
        "cL x true Right\n"
        "Observed Right = cR x true Left + "
        "(1 - cR) x true Right",
        ha="center",
        va="center",
        family="monospace",
        fontsize=10,
        bbox=dict(
            boxstyle="round,pad=0.5",
            ec="#AAAAAA",
            fc="#F7F7F7",
        ),
    )


def calculate_site_evidence(
    matched,
    observed,
    model159,
    model160,
):
    emission = model159.load_emission()

    arrays = {
        "left_alt":
            matched["left_alt_depth"].to_numpy(float),
        "left_total":
            matched["left_AD_total"].to_numpy(float),
        "right_alt":
            matched["right_alt_depth"].to_numpy(float),
        "right_total":
            matched["right_AD_total"].to_numpy(float),
    }

    null = observed.loc[
        observed["fit"] == "three_founder_null"
    ].iloc[0]

    alt = observed.loc[
        observed["fit"]
        == "opposing_state_weight_alternative"
    ].iloc[0]

    ll_null = model159.state_log_likelihoods(
        arrays,
        float(null["c_left"]),
        float(null["c_right"]),
        emission,
    )
    ll_alt = model159.state_log_likelihoods(
        arrays,
        float(alt["c_left"]),
        float(alt["c_right"]),
        emission,
    )

    null_site = model160.marginal_evidence(ll_null)
    alt_site = model160.marginal_evidence(
        ll_alt,
        float(alt["epsilon"]),
    )

    epsilon = float(alt["epsilon"])

    log_prior = np.empty(9, dtype=float)
    log_prior[model160.COMPATIBLE] = (
        math.log1p(-epsilon)
        - math.log(len(model160.COMPATIBLE))
    )
    log_prior[model160.OPPOSING] = (
        math.log(epsilon)
        - math.log(len(model160.OPPOSING))
    )

    log_joint = (
        ll_alt
        + log_prior[None, :]
    )
    posterior = np.exp(
        log_joint
        - logsumexp(
            log_joint,
            axis=1,
        )[:, None]
    )

    site = matched[
        KEYS + ["block_id"]
    ].copy()

    site["left_AF"] = (
        arrays["left_alt"]
        / arrays["left_total"]
    )
    site["right_AF"] = (
        arrays["right_alt"]
        / arrays["right_total"]
    )
    site["nested_LR_contribution"] = (
        2
        * (
            alt_site
            - null_site
        )
    )
    site["opposing_state_posterior"] = (
        posterior[
            :,
            model160.OPPOSING,
        ].sum(axis=1)
    )
    site["block_start"] = (
        site["pos"].astype(int)
        // 10_000
    ) * 10_000

    block = site.groupby(
        [
            "chrom",
            "block_start",
            "block_id",
        ],
        as_index=False,
    ).agg(
        sites=(
            "pos",
            "size",
        ),
        nested_LR_contribution=(
            "nested_LR_contribution",
            "sum",
        ),
        opposing_posterior_sum=(
            "opposing_state_posterior",
            "sum",
        ),
        opposing_posterior_mean=(
            "opposing_state_posterior",
            "mean",
        ),
    )

    return site, block, emission


def plot_af_panel(
    ax,
    site,
    observed,
    emission,
    model159,
    show_legend=True,
):
    hb = ax.hexbin(
        site["left_AF"],
        site["right_AF"],
        gridsize=43,
        mincnt=1,
        bins="log",
        cmap="Greys",
        linewidths=0,
    )

    alt = observed.loc[
        observed["fit"]
        == "opposing_state_weight_alternative"
    ].iloc[0]

    c_left = float(alt["c_left"])
    c_right = float(alt["c_right"])

    gl = model159.GL
    gr = model159.GR

    expected_left = model159.map_fraction(
        (1 - c_left) * gl
        + c_left * gr,
        emission["error_rate"],
        emission["bias_log_odds"],
    )
    expected_right = model159.map_fraction(
        c_right * gl
        + (1 - c_right) * gr,
        emission["error_rate"],
        emission["bias_log_odds"],
    )

    opposing = np.array(
        [2, 6]
    )
    compatible = np.array([
        index
        for index in range(9)
        if index not in opposing
    ])

    ax.scatter(
        expected_left[compatible],
        expected_right[compatible],
        marker="s",
        s=48,
        facecolor="white",
        edgecolor=BLUE,
        linewidth=1.2,
        label=(
            "Three-founder-compatible centroids"
        ),
        zorder=4,
    )
    ax.scatter(
        expected_left[opposing],
        expected_right[opposing],
        marker="D",
        s=58,
        facecolor=ORANGE,
        edgecolor="white",
        linewidth=0.8,
        label=(
            "Four-founder-only centroids"
        ),
        zorder=5,
    )

    for index in opposing:
        ax.annotate(
            model159.STATE_LABELS[index],
            (
                expected_left[index],
                expected_right[index],
            ),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
            color=ORANGE,
            weight="bold",
        )

    ax.set(
        xlim=(-0.02, 1.02),
        ylim=(-0.02, 1.02),
        xlabel=(
            "Observed alternate-read fraction: Left"
        ),
        ylabel=(
            "Observed alternate-read fraction: Right"
        ),
    )
    ax.set_aspect(
        "equal",
        adjustable="box",
    )
    ax.set_title(
        "Observed read fractions and "
        "fitted genotype-pair centroids"
    )

    if show_legend:
        ax.legend(
            frameon=False,
            loc="upper left",
            fontsize=7,
        )

    return hb


def plot_calibration(
    ax,
    simulations,
    mixing,
    founder,
    side,
):
    if founder == "three_founder":
        prefix = "null"
        color = BLUE
    else:
        prefix = "alternative"
        color = ORANGE

    subset = simulations.loc[
        simulations["true_founder_model"]
        == founder
    ].copy()

    xcol = f"true_c_{side}"
    ycol = f"{prefix}_c_{side}"

    grouped = subset.groupby(xcol)[ycol]
    x = np.array(
        sorted(
            subset[xcol].unique()
        )
    )
    median = (
        grouped.median()
        .reindex(x)
        .to_numpy()
    )
    low = (
        grouped.quantile(0.025)
        .reindex(x)
        .to_numpy()
    )
    high = (
        grouped.quantile(0.975)
        .reindex(x)
        .to_numpy()
    )

    ax.scatter(
        subset[xcol],
        subset[ycol],
        s=7,
        alpha=0.12,
        color=color,
        rasterized=True,
    )
    ax.fill_between(
        x,
        low,
        high,
        color=color,
        alpha=0.16,
        linewidth=0,
    )
    ax.plot(
        x,
        median,
        color=color,
        linewidth=2,
        label="simulation median",
    )
    ax.plot(
        [0, 0.42],
        [0, 0.42],
        color=GREY,
        linestyle="--",
        linewidth=1,
        label="perfect recovery",
    )

    row = mixing.loc[
        (
            mixing["founder_model"]
            == founder
        )
        & (
            mixing["side"]
            == side
        )
    ].iloc[0]

    ax.scatter(
        [row["calibrated_estimate"]],
        [row["observed_raw_estimate"]],
        marker="*",
        s=120,
        color=BLACK,
        zorder=6,
        label="observed projection",
    )
    ax.set(
        xlim=(-0.01, 0.42),
        ylim=(-0.01, 0.46),
        xlabel=(
            f"True mixing into "
            f"{side.title()}"
        ),
        ylabel=(
            f"Recovered mixing into "
            f"{side.title()}"
        ),
    )
    ax.set_title(
        f"{founder.replace('_', ' ').title()}: "
        f"{side.title()}"
    )
    clean_axis(ax)


def nearest_simulations(
    simulations,
    founder,
    number,
    high_call=False,
):
    subset = simulations.loc[
        simulations["true_founder_model"]
        == founder
    ].copy()

    if high_call:
        subset = subset.loc[
            subset[
                "both_donors_high_call"
            ].astype(bool)
        ]

    return subset.nsmallest(
        min(
            number,
            len(subset),
        ),
        "distance_from_corresponding_observed_mixing",
    )


def plot_founder_panel(
    ax,
    simulations,
    observed_lr,
    threshold,
    high_call=False,
):
    number = (
        40
        if high_call
        else 60
    )

    m3 = nearest_simulations(
        simulations,
        "three_founder",
        number,
        high_call,
    )
    m4 = nearest_simulations(
        simulations,
        "four_founder",
        number,
        high_call,
    )

    def transform(values):
        return np.log10(
            1
            + np.asarray(
                values,
                dtype=float,
            )
        )

    all_x = np.concatenate([
        transform(
            m3["founder_LR_statistic"]
        ),
        transform(
            m4["founder_LR_statistic"]
        ),
        transform(
            [observed_lr]
        ),
    ])

    bins = np.linspace(
        all_x.min(),
        all_x.max(),
        24,
    )

    ax.hist(
        transform(
            m3["founder_LR_statistic"]
        ),
        bins=bins,
        density=True,
        histtype="stepfilled",
        alpha=0.42,
        color=BLUE,
        label=(
            f"Three-founder (n={len(m3)})"
        ),
    )
    ax.hist(
        transform(
            m4["founder_LR_statistic"]
        ),
        bins=bins,
        density=True,
        histtype="step",
        linewidth=2,
        color=ORANGE,
        label=(
            f"Four-founder (n={len(m4)})"
        ),
    )
    ax.axvline(
        transform([threshold])[0],
        color=GREY,
        linestyle="--",
        linewidth=1.3,
        label=(
            "three-founder 95% threshold"
        ),
    )
    ax.axvline(
        transform([observed_lr])[0],
        color=BLACK,
        linewidth=2,
        label=(
            f"observed LR = "
            f"{observed_lr:.1f}"
        ),
    )

    tick_values = np.array([
        0,
        0.1,
        1,
        10,
        100,
        1000,
    ])
    visible = tick_values[
        np.log10(
            1
            + tick_values
        )
        <= all_x.max() + 0.1
    ]

    ax.set_xticks(
        np.log10(
            1
            + visible
        ),
        [
            f"{value:g}"
            for value in visible
        ],
    )
    ax.set_xlabel(
        "Nested LR statistic "
        "(log-scaled spacing)"
    )
    ax.set_ylabel(
        "Simulation density"
    )
    ax.set_title(
        "High-call donor sensitivity"
        if high_call
        else "All 47 controls"
    )

    clean_axis(ax)

    ax.legend(
        frameon=False,
        fontsize=7,
    )

    exceedances = int(
        (
            m3["founder_LR_statistic"]
            >= observed_lr
        ).sum()
    )
    mc_p = (
        1
        + exceedances
    ) / (
        1
        + len(m3)
    )
    power = float(
        (
            m4["founder_LR_statistic"]
            > threshold
        ).mean()
    )

    ax.text(
        0.98,
        0.63,
        f"Null exceedances: "
        f"{exceedances}/{len(m3)}\n"
        f"Monte Carlo p = {mc_p:.4f}\n"
        f"Detection fraction = {power:.2f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7.5,
        bbox=dict(
            boxstyle="round,pad=0.3",
            fc="white",
            ec="#BBBBBB",
            alpha=0.9,
        ),
    )


def interval_rows(mixing, bootstrap):
    rows = []

    combinations = [
        (
            "three_founder",
            "null",
        ),
        (
            "four_founder",
            "alternative",
        ),
    ]

    for founder, prefix in combinations:
        for side in ["left", "right"]:
            row = mixing.loc[
                (
                    mixing["founder_model"]
                    == founder
                )
                & (
                    mixing["side"]
                    == side
                )
            ].iloc[0]

            boot = bootstrap[
                f"{prefix}_c_{side}_calibrated"
            ]

            rows.append({
                "founder":
                    founder,
                "side":
                    side,
                "estimate":
                    row["calibrated_estimate"],
                "simulation_low":
                    row[
                        "simulation_calibrated_95_low"
                    ],
                "simulation_high":
                    row[
                        "simulation_calibrated_95_high"
                    ],
                "bootstrap_low":
                    boot.quantile(0.025),
                "bootstrap_high":
                    boot.quantile(0.975),
            })

    return pd.DataFrame(rows)


def plot_interval_panel(ax, intervals):
    order = [
        (
            "four_founder",
            "right",
        ),
        (
            "four_founder",
            "left",
        ),
        (
            "three_founder",
            "right",
        ),
        (
            "three_founder",
            "left",
        ),
    ]

    labels = {
        (
            "four_founder",
            "right",
        ):
            "Four-founder: Left -> Right",
        (
            "four_founder",
            "left",
        ):
            "Four-founder: Right -> Left",
        (
            "three_founder",
            "right",
        ):
            "Three-founder: Left -> Right",
        (
            "three_founder",
            "left",
        ):
            "Three-founder: Right -> Left",
    }

    for y, key in enumerate(order):
        row = intervals.loc[
            (
                intervals["founder"]
                == key[0]
            )
            & (
                intervals["side"]
                == key[1]
            )
        ].iloc[0]

        color = (
            ORANGE
            if key[0] == "four_founder"
            else BLUE
        )

        ax.plot(
            [
                100
                * row["simulation_low"],
                100
                * row["simulation_high"],
            ],
            [y, y],
            color=color,
            alpha=0.35,
            linewidth=3,
            solid_capstyle="round",
        )
        ax.plot(
            [
                100
                * row["bootstrap_low"],
                100
                * row["bootstrap_high"],
            ],
            [y, y],
            color=color,
            linewidth=7,
            solid_capstyle="round",
        )
        ax.scatter(
            100
            * row["estimate"],
            y,
            s=55,
            color=color,
            edgecolor="white",
            linewidth=0.8,
            zorder=5,
        )
        ax.text(
            100
            * row["estimate"]
            + 0.6,
            y,
            f"{100 * row['estimate']:.1f}%",
            va="center",
            fontsize=8,
        )

    ax.set_yticks(
        range(4),
        [
            labels[key]
            for key in order
        ],
    )
    ax.set_xlabel(
        "Estimated reciprocal mixing (%)"
    )
    ax.set_title(
        "Bias-corrected estimates and "
        "two uncertainty sources"
    )

    clean_axis(ax)

    ax.legend(
        handles=[
            Patch(
                facecolor=GREY,
                alpha=0.35,
                label=(
                    "simulation-recovery "
                    "95% interval"
                ),
            ),
            Patch(
                facecolor=GREY,
                label=(
                    "calibrated block-bootstrap "
                    "95% interval"
                ),
            ),
        ],
        frameon=False,
        fontsize=7,
        loc="lower right",
    )


def plot_bootstrap_panel(ax, bootstrap):
    combinations = [
        (
            "Three-founder",
            "null",
            BLUE,
        ),
        (
            "Four-founder",
            "alternative",
            ORANGE,
        ),
    ]

    for founder, prefix, color in combinations:
        x = (
            bootstrap[
                f"{prefix}_c_left_calibrated"
            ]
            * 100
        )
        y = (
            bootstrap[
                f"{prefix}_c_right_calibrated"
            ]
            * 100
        )

        ax.scatter(
            x,
            y,
            s=12,
            alpha=0.22,
            color=color,
            label=founder,
            rasterized=True,
        )
        ax.scatter(
            x.median(),
            y.median(),
            marker="*",
            s=110,
            color=color,
            edgecolor="white",
            linewidth=0.8,
            zorder=5,
        )

    limits = [
        min(
            ax.get_xlim()[0],
            ax.get_ylim()[0],
        ),
        max(
            ax.get_xlim()[1],
            ax.get_ylim()[1],
        ),
    ]

    ax.plot(
        limits,
        limits,
        color=GREY,
        linestyle="--",
        linewidth=1,
    )
    ax.set(
        xlim=limits,
        ylim=limits,
        xlabel="Right material in Left (%)",
        ylabel="Left material in Right (%)",
        title=(
            "Paired 10-kb "
            "block-bootstrap estimates"
        ),
    )

    clean_axis(ax)

    ax.legend(
        frameon=False,
    )


def plot_block_panel(
    ax,
    block,
    posterior=False,
):
    x = (
        block["block_start"]
        / 1e6
    )
    width = 0.009

    if posterior:
        ax.bar(
            x,
            block[
                "opposing_posterior_sum"
            ],
            width=width,
            color=ORANGE,
            alpha=0.75,
            edgecolor="none",
        )
        ax.set_ylabel(
            "Summed opposing-state posterior"
        )
    else:
        colors = np.where(
            block[
                "nested_LR_contribution"
            ] >= 0,
            ORANGE,
            BLUE,
        )
        ax.bar(
            x,
            block[
                "nested_LR_contribution"
            ],
            width=width,
            color=colors,
            alpha=0.82,
            edgecolor="none",
        )
        ax.axhline(
            0,
            color=BLACK,
            linewidth=0.8,
        )
        ax.set_ylabel(
            "Contribution to nested LR"
        )

    ax.axvspan(
        OPTIX_START / 1e6,
        OPTIX_END / 1e6,
        color=GREEN,
        alpha=0.16,
        label="optix",
    )
    ax.set_xlabel(
        "Hmr18 position (Mb)"
    )

    clean_axis(ax)


def create_narrative(
    site,
    block,
    observed,
    simulations,
    mixing,
    bootstrap,
    donor,
    summary,
    emission,
    model159,
):
    observed_lr = as_float(
        summary,
        "observed_nested_LR_statistic",
    )
    threshold = as_float(
        summary,
        "three_founder_simulated_95_percent_threshold",
    )
    intervals = interval_rows(
        mixing,
        bootstrap,
    )

    with PdfPages(NARRATIVE_PDF) as pdf:
        figure = plt.figure(
            figsize=(11.7, 8.3)
        )
        add_page_heading(
            figure,
            1,
            "Define the biological alternatives",
            "Founder sharing and sequencing mixture "
            "are distinct parts of the model.",
        )
        axis = figure.add_axes([
            0.05,
            0.08,
            0.90,
            0.74,
        ])
        model_schematic(axis)
        pdf.savefig(
            figure,
            bbox_inches="tight",
        )
        plt.close(figure)

        figure, axis = plt.subplots(
            figsize=(8.6, 8.0)
        )
        figure.subplots_adjust(
            top=0.78,
            left=0.12,
            right=0.88,
            bottom=0.11,
        )
        add_page_heading(
            figure,
            2,
            "Fit the observed Left and Right read fractions",
            "No focal GT calls are used: the model "
            "works directly from allele depths.",
        )
        plot_af_panel(
            axis,
            site,
            observed,
            emission,
            model159,
        )
        axis.set_title("")
        axis.legend(
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.10),
            ncol=2,
            fontsize=7,
        )
        pdf.savefig(
            figure,
            bbox_inches="tight",
        )
        plt.close(figure)

        figure, axes = plt.subplots(
            2,
            2,
            figsize=(11.7, 8.3),
        )
        figure.subplots_adjust(
            top=0.78,
            left=0.08,
            right=0.97,
            bottom=0.09,
            hspace=0.38,
            wspace=0.28,
        )
        add_page_heading(
            figure,
            3,
            "Calibrate how accurately mixing is recovered",
            "Each point is a simulated control pair; "
            "bands span the recovery distribution.",
        )

        founders = [
            "three_founder",
            "three_founder",
            "four_founder",
            "four_founder",
        ]
        sides = [
            "left",
            "right",
            "left",
            "right",
        ]

        for axis, founder, side in zip(
            axes.ravel(),
            founders,
            sides,
        ):
            plot_calibration(
                axis,
                simulations,
                mixing,
                founder,
                side,
            )

        axes[0, 0].legend(
            frameon=False,
            fontsize=7,
        )
        pdf.savefig(
            figure,
            bbox_inches="tight",
        )
        plt.close(figure)

        figure, axes = plt.subplots(
            1,
            2,
            figsize=(11.7, 6.3),
        )
        figure.subplots_adjust(
            top=0.72,
            left=0.08,
            right=0.97,
            bottom=0.14,
            wspace=0.25,
        )
        add_page_heading(
            figure,
            4,
            "Calibrate evidence for four founders",
            "The observed statistic is compared only "
            "with simulations near the observed mixing.",
        )
        plot_founder_panel(
            axes[0],
            simulations,
            observed_lr,
            threshold,
            False,
        )
        high_threshold = as_float(
            summary,
            "high_call_three_founder_"
            "simulated_95_percent_threshold",
        )
        plot_founder_panel(
            axes[1],
            simulations,
            observed_lr,
            high_threshold,
            True,
        )
        pdf.savefig(
            figure,
            bbox_inches="tight",
        )
        plt.close(figure)

        figure, axes = plt.subplots(
            1,
            2,
            figsize=(11.7, 6.3),
        )
        figure.subplots_adjust(
            top=0.72,
            left=0.15,
            right=0.97,
            bottom=0.14,
            wspace=0.32,
        )
        add_page_heading(
            figure,
            5,
            "Report reciprocal mixing with uncertainty",
            "Thick intervals resample genomic blocks; "
            "thin intervals include simulation recovery.",
        )
        plot_interval_panel(
            axes[0],
            intervals,
        )
        plot_bootstrap_panel(
            axes[1],
            bootstrap,
        )
        pdf.savefig(
            figure,
            bbox_inches="tight",
        )
        plt.close(figure)

        figure, axes = plt.subplots(
            2,
            1,
            figsize=(11.7, 7.5),
            sharex=True,
        )
        figure.subplots_adjust(
            top=0.76,
            left=0.10,
            right=0.97,
            bottom=0.10,
            hspace=0.16,
        )
        add_page_heading(
            figure,
            6,
            "Locate the founder-model signal along Hmr18",
            "Coherent multi-block signal is more persuasive "
            "than isolated site-level deviations.",
        )
        plot_block_panel(
            axes[0],
            block,
            False,
        )
        axes[0].set_xlabel("")
        axes[0].set_title(
            "Fixed-fit block contributions; "
            f"total LR = "
            f"{block['nested_LR_contribution'].sum():.2f}"
        )
        plot_block_panel(
            axes[1],
            block,
            True,
        )
        axes[1].legend(
            frameon=False,
            loc="upper right",
        )
        pdf.savefig(
            figure,
            bbox_inches="tight",
        )
        plt.close(figure)


def create_summary(
    site,
    block,
    observed,
    simulations,
    mixing,
    bootstrap,
    summary,
    emission,
    model159,
):
    observed_lr = as_float(
        summary,
        "observed_nested_LR_statistic",
    )
    threshold = as_float(
        summary,
        "three_founder_simulated_95_percent_threshold",
    )
    intervals = interval_rows(
        mixing,
        bootstrap,
    )

    figure = plt.figure(
        figsize=(11.7, 9.0)
    )
    grid = figure.add_gridspec(
        2,
        2,
        left=0.08,
        right=0.97,
        bottom=0.08,
        top=0.87,
        wspace=0.30,
        hspace=0.34,
    )
    axes = [
        figure.add_subplot(
            grid[row, column]
        )
        for row in range(2)
        for column in range(2)
    ]

    figure.suptitle(
        "Reciprocal mixing and founder-model "
        "calibration at Hmr18/optix",
        fontsize=16,
        weight="bold",
        y=0.975,
    )

    model_schematic(
        axes[0],
        compact=True,
    )
    axes[0].set_title(
        "A  Biological hypotheses",
        loc="left",
        weight="bold",
    )

    plot_founder_panel(
        axes[1],
        simulations,
        observed_lr,
        threshold,
        False,
    )
    axes[1].set_title(
        "",
        loc="center",
    )
    axes[1].set_title(
        "B  Simulation-calibrated founder test",
        loc="left",
        weight="bold",
    )

    plot_interval_panel(
        axes[2],
        intervals,
    )
    axes[2].set_title(
        "",
        loc="center",
    )
    axes[2].set_title(
        "C  Reciprocal mixing estimates",
        loc="left",
        weight="bold",
    )

    plot_block_panel(
        axes[3],
        block,
        False,
    )
    axes[3].set_title(
        "D  Genomic distribution of founder signal",
        loc="left",
        weight="bold",
    )

    figure.savefig(
        SUMMARY_PDF,
        bbox_inches="tight",
    )
    figure.savefig(
        SUMMARY_PNG,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)


def main():
    required = [
        MODEL159,
        MODEL160,
        MATCHED_FILE,
        OBSERVED_FILE,
        SIM_FILE,
        MIXING_FILE,
        BOOTSTRAP_FILE,
        FINAL_FILE,
        DONOR_FILE,
    ]

    missing = [
        str(path)
        for path in required
        if not path.is_file()
    ]

    if missing:
        raise RuntimeError(
            "Missing required inputs:\n"
            + "\n".join(missing)
        )

    OUTDIR.mkdir(
        parents=True,
        exist_ok=True,
    )
    set_style()

    model159 = load_module(
        MODEL159,
        "model159_for_narrative",
    )
    model160 = load_module(
        MODEL160,
        "model160_for_narrative",
    )

    matched = pd.read_csv(
        MATCHED_FILE,
        sep="\t",
        low_memory=False,
    )
    observed = pd.read_csv(
        OBSERVED_FILE,
        sep="\t",
    )
    simulations = pd.read_csv(
        SIM_FILE,
        sep="\t",
        low_memory=False,
    )
    mixing = pd.read_csv(
        MIXING_FILE,
        sep="\t",
    )
    bootstrap = pd.read_csv(
        BOOTSTRAP_FILE,
        sep="\t",
    )
    donor = pd.read_csv(
        DONOR_FILE,
        sep="\t",
    )
    summary = numeric_summary(
        FINAL_FILE
    )

    site, block, emission = calculate_site_evidence(
        matched,
        observed,
        model159,
        model160,
    )

    site.to_csv(
        SITE_OUT,
        sep="\t",
        index=False,
        compression="gzip",
    )
    block.to_csv(
        BLOCK_OUT,
        sep="\t",
        index=False,
    )

    observed_lr = as_float(
        summary,
        "observed_nested_LR_statistic",
    )
    contribution_sum = float(
        block[
            "nested_LR_contribution"
        ].sum()
    )

    if not np.isclose(
        contribution_sum,
        observed_lr,
        rtol=1e-6,
        atol=1e-5,
    ):
        raise RuntimeError(
            "Block LR contributions sum to "
            f"{contribution_sum}, "
            f"expected {observed_lr}"
        )

    create_narrative(
        site,
        block,
        observed,
        simulations,
        mixing,
        bootstrap,
        donor,
        summary,
        emission,
        model159,
    )
    create_summary(
        site,
        block,
        observed,
        simulations,
        mixing,
        bootstrap,
        summary,
        emission,
        model159,
    )

    audit = pd.DataFrame(
        [
            (
                "matched_sites",
                len(site),
            ),
            (
                "ten_kb_blocks",
                len(block),
            ),
            (
                "observed_nested_LR",
                observed_lr,
            ),
            (
                "block_contribution_sum",
                contribution_sum,
            ),
            (
                "controls",
                donor["sample"].nunique(),
            ),
            (
                "three_founder_simulations",
                (
                    simulations[
                        "true_founder_model"
                    ]
                    == "three_founder"
                ).sum(),
            ),
            (
                "four_founder_simulations",
                (
                    simulations[
                        "true_founder_model"
                    ]
                    == "four_founder"
                ).sum(),
            ),
            (
                "bootstrap_replicates",
                len(bootstrap),
            ),
            (
                "focal_GT_used",
                False,
            ),
            (
                "narrative_pdf",
                NARRATIVE_PDF,
            ),
            (
                "summary_pdf",
                SUMMARY_PDF,
            ),
        ],
        columns=[
            "item",
            "value",
        ],
    )

    audit.to_csv(
        AUDIT_OUT,
        sep="\t",
        index=False,
    )

    print(
        "NARRATIVE VISUALIZATION COMPLETE"
    )
    print(
        audit.to_string(index=False)
    )

    outputs = [
        NARRATIVE_PDF,
        SUMMARY_PDF,
        SUMMARY_PNG,
        SITE_OUT,
        BLOCK_OUT,
        AUDIT_OUT,
    ]

    for path in outputs:
        if (
            not path.is_file()
            or path.stat().st_size == 0
        ):
            raise RuntimeError(
                f"Missing or empty output: {path}"
            )

        print(
            f"PRESENT\t{path}\t"
            f"{path.stat().st_size} bytes"
        )


if __name__ == "__main__":
    main()
