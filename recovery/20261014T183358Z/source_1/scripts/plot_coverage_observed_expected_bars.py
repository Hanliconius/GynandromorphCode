#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# ============================================================
# INPUT
# ============================================================

RUN = Path(
    "analysis/chromosome_coverage_bars/"
    "run_20260923T162015902039Z"
)

INPUT = RUN / "chromosome_coverage.tsv"

if not INPUT.is_file():
    raise FileNotFoundError(INPUT)

d = pd.read_csv(
    INPUT,
    sep="\t",
)

required = {
    "label",
    "autosome",
    "acgt_bases",
    "left_t2_relative_depth",
    "left_t3_relative_depth",
    "right_t2_relative_depth",
    "right_t3_relative_depth",
}

missing = required - set(d.columns)

if missing:
    raise RuntimeError(
        "chromosome_coverage.tsv is missing columns: "
        + ", ".join(sorted(missing))
    )


# ============================================================
# SAMPLE DEFINITIONS
# ============================================================

SAMPLES = [
    "left_t2",
    "left_t3",
    "right_t2",
    "right_t3",
]

DISPLAY = {
    "left_t2": "T2L",
    "left_t3": "T3L",
    "right_t2": "T2R",
    "right_t3": "T3R",
}

# Keep Left related blues and Right related oranges,
# with T2/T3 distinguished by shade.
COLOURS = {
    "left_t2": "#2478A7",
    "left_t3": "#62B0D2",
    "right_t2": "#E6842A",
    "right_t3": "#F2B56B",
}


# ============================================================
# EXPECTED RELATIVE COVERAGE
#
# Relative to diploid autosomal coverage:
#
#   Left  = ZW
#   Right = ZZ
#
# Autosomes: 2 copies on both sides -> 1.0
# W:         1 copy Left, 0 Right   -> 0.5 / 0
# Z:         1 copy Left, 2 Right   -> 0.5 / 1
# ============================================================

EXPECTED = {
    "Autosomes": {
        "left_t2": 1.0,
        "left_t3": 1.0,
        "right_t2": 1.0,
        "right_t3": 1.0,
    },
    "W": {
        "left_t2": 0.5,
        "left_t3": 0.5,
        "right_t2": 0.0,
        "right_t3": 0.0,
    },
    "Z": {
        "left_t2": 0.5,
        "left_t3": 0.5,
        "right_t2": 1.0,
        "right_t3": 1.0,
    },
}


# ============================================================
# OBSERVED VALUES
# ============================================================

autosomes = d.loc[
    d["autosome"].astype(str).str.lower().eq("true")
].copy()

if autosomes.empty:
    raise RuntimeError(
        "No autosomes found in chromosome_coverage.tsv"
    )

w_rows = d.loc[
    d["label"].astype(str).str.upper().eq("W")
]

z_rows = d.loc[
    d["label"].astype(str).str.upper().eq("Z")
]

if len(w_rows) != 1:
    raise RuntimeError(
        f"Expected exactly one W row; found {len(w_rows)}"
    )

if len(z_rows) != 1:
    raise RuntimeError(
        f"Expected exactly one Z row; found {len(z_rows)}"
    )

w_row = w_rows.iloc[0]
z_row = z_rows.iloc[0]


def weighted_autosomal_relative_depth(sample):
    column = f"{sample}_relative_depth"

    values = pd.to_numeric(
        autosomes[column],
        errors="raise",
    ).to_numpy(dtype=float)

    weights = pd.to_numeric(
        autosomes["acgt_bases"],
        errors="raise",
    ).to_numpy(dtype=float)

    if (
        not np.isfinite(values).all()
        or not np.isfinite(weights).all()
        or np.any(weights <= 0)
    ):
        raise RuntimeError(
            f"Invalid autosomal values for {sample}"
        )

    return float(
        np.average(
            values,
            weights=weights,
        )
    )


OBSERVED = {
    "Autosomes": {},
    "W": {},
    "Z": {},
}

for sample in SAMPLES:

    column = f"{sample}_relative_depth"

    OBSERVED["Autosomes"][sample] = (
        weighted_autosomal_relative_depth(sample)
    )

    OBSERVED["W"][sample] = float(
        w_row[column]
    )

    OBSERVED["Z"][sample] = float(
        z_row[column]
    )


# ============================================================
# WRITE EXACT VALUES USED IN FIGURE
# ============================================================

rows = []

for chromosome_class in [
    "Autosomes",
    "W",
    "Z",
]:
    for sample in SAMPLES:

        observed = OBSERVED[
            chromosome_class
        ][sample]

        expected = EXPECTED[
            chromosome_class
        ][sample]

        rows.append(
            {
                "chromosome_class": chromosome_class,
                "sample": DISPLAY[sample],
                "sample_key": sample,
                "observed_relative_coverage": observed,
                "expected_relative_coverage": expected,
                "deviation": observed - expected,
            }
        )

summary = pd.DataFrame(rows)

summary_file = (
    RUN
    / "coverage_observed_expected_summary.tsv"
)

summary.to_csv(
    summary_file,
    sep="\t",
    index=False,
)


# ============================================================
# FIGURE
# ============================================================

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.linewidth": 0.7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

fig, axes = plt.subplots(
    1,
    3,
    figsize=(10.6, 4.0),
    sharey=True,
)

fig.subplots_adjust(
    left=0.09,
    right=0.985,
    bottom=0.22,
    top=0.83,
    wspace=0.17,
)

classes = [
    "Autosomes",
    "W",
    "Z",
]

x = np.arange(
    len(SAMPLES),
)

observed_max = max(
    value
    for chromosome_class
    in OBSERVED.values()
    for value
    in chromosome_class.values()
)

expected_max = max(
    value
    for chromosome_class
    in EXPECTED.values()
    for value
    in chromosome_class.values()
)

ymax = max(
    1.15,
    observed_max + 0.14,
    expected_max + 0.14,
)


for ax, chromosome_class in zip(
    axes,
    classes,
):

    observed = np.array(
        [
            OBSERVED[
                chromosome_class
            ][sample]
            for sample in SAMPLES
        ],
        dtype=float,
    )

    expected = np.array(
        [
            EXPECTED[
                chromosome_class
            ][sample]
            for sample in SAMPLES
        ],
        dtype=float,
    )

    colours = [
        COLOURS[sample]
        for sample in SAMPLES
    ]

    # --------------------------------------------------------
    # OBSERVED COVERAGE
    # --------------------------------------------------------

    bars = ax.bar(
        x,
        observed,
        width=0.66,
        color=colours,
        edgecolor="white",
        linewidth=0.6,
        zorder=3,
    )

    # --------------------------------------------------------
    # EXPECTED COVERAGE
    # --------------------------------------------------------

    ax.plot(
        x,
        expected,
        color="#222222",
        linestyle=(0, (4, 2)),
        linewidth=1.25,
        zorder=6,
    )

    ax.scatter(
        x,
        expected,
        s=30,
        facecolor="white",
        edgecolor="#222222",
        linewidth=1.1,
        zorder=7,
    )

    # --------------------------------------------------------
    # DEVIATION FROM EXPECTATION
    # --------------------------------------------------------

    for xi, obs, exp in zip(
        x,
        observed,
        expected,
    ):

        deviation = obs - exp

        ax.plot(
            [xi, xi],
            [obs, exp],
            color="#444444",
            linewidth=1.0,
            zorder=5,
        )

        # Place delta label away from both observed and expected.
        if deviation >= 0:
            text_y = max(obs, exp) + 0.035
            va = "bottom"
        else:
            text_y = min(obs, exp) - 0.035
            va = "top"

        # Don't put negative-position labels below the axis.
        if text_y < 0.015:
            text_y = max(obs, exp) + 0.035
            va = "bottom"

        ax.text(
            xi,
            text_y,
            f"Δ {deviation:+.2f}",
            ha="center",
            va=va,
            fontsize=7.5,
            color="#444444",
            zorder=8,
        )

    # --------------------------------------------------------
    # PANEL FORMATTING
    # --------------------------------------------------------

    ax.set_title(
        chromosome_class,
        fontsize=12,
        fontweight="bold",
        pad=9,
    )

    ax.set_xticks(
        x,
    )

    ax.set_xticklabels(
        [
            DISPLAY[sample]
            for sample in SAMPLES
        ],
        fontsize=9,
    )

    ax.set_ylim(
        0,
        ymax,
    )

    ax.grid(
        axis="y",
        color="#E5E5E5",
        linewidth=0.65,
        zorder=0,
    )

    ax.set_axisbelow(
        True
    )

    ax.spines[
        "top"
    ].set_visible(False)

    ax.spines[
        "right"
    ].set_visible(False)

    ax.tick_params(
        axis="both",
        labelsize=9,
    )


axes[0].set_ylabel(
    "Relative coverage",
    fontsize=10.5,
)


# ============================================================
# TITLE / LEGEND
# ============================================================

fig.suptitle(
    "Observed chromosome coverage relative to dosage expectation",
    fontsize=14,
    fontweight="bold",
    y=0.965,
)

sample_handles = [
    Patch(
        facecolor=COLOURS[sample],
        edgecolor="none",
        label=DISPLAY[sample],
    )
    for sample in SAMPLES
]

expected_handle = Line2D(
    [0],
    [0],
    color="#222222",
    linestyle=(0, (4, 2)),
    linewidth=1.25,
    marker="o",
    markerfacecolor="white",
    markeredgecolor="#222222",
    markersize=5,
    label="Expected dosage",
)

fig.legend(
    handles=(
        sample_handles
        + [expected_handle]
    ),
    loc="lower center",
    bbox_to_anchor=(
        0.5,
        0.035,
    ),
    ncol=5,
    frameon=False,
    fontsize=9,
)


# ============================================================
# SAVE
# ============================================================

pdf_file = (
    RUN
    / "coverage_observed_vs_expected_bars.pdf"
)

png_file = (
    RUN
    / "coverage_observed_vs_expected_bars.png"
)

fig.savefig(
    pdf_file,
)

fig.savefig(
    png_file,
    dpi=300,
)

plt.close(fig)


# ============================================================
# REPORT
# ============================================================

print()
print(
    "VALUES USED"
)

print(
    summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.6f}",
    )
)

print()
print(
    "WROTE",
    summary_file.resolve(),
    sep="\t",
)

print(
    "WROTE",
    pdf_file.resolve(),
    sep="\t",
)

print(
    "WROTE",
    png_file.resolve(),
    sep="\t",
)
