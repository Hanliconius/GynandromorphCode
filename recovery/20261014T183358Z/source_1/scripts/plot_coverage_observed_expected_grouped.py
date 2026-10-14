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

d = pd.read_csv(INPUT, sep="\t")

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

COLOURS = {
    "left_t2": "#2478A7",
    "left_t3": "#62B0D2",
    "right_t2": "#E6842A",
    "right_t3": "#F2B56B",
}

CLASSES = [
    "Autosomes",
    "W",
    "Z",
]


# ============================================================
# EXPECTED RELATIVE COVERAGE
#
# Left = ZW ; Right = ZZ
# Autosomes relative expectation is 1 for all samples.
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
    raise RuntimeError("No autosomes found.")

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

    return float(np.average(values, weights=weights))


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
    OBSERVED["W"][sample] = float(w_row[column])
    OBSERVED["Z"][sample] = float(z_row[column])


# ============================================================
# SAVE THE EXACT NUMBERS USED
# ============================================================

rows = []

for chromosome_class in CLASSES:
    for sample in SAMPLES:
        observed = OBSERVED[chromosome_class][sample]
        expected = EXPECTED[chromosome_class][sample]

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

summary_file = RUN / "coverage_observed_expected_grouped_summary.tsv"
summary.to_csv(summary_file, sep="\t", index=False)


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

fig = plt.figure(figsize=(8.6, 4.8))
ax = fig.add_axes([0.10, 0.20, 0.86, 0.64])

group_centres = np.arange(len(CLASSES), dtype=float)
bar_width = 0.16
offsets = np.array([-1.5, -0.5, 0.5, 1.5]) * bar_width

observed_max = max(
    value
    for chromosome_class in OBSERVED.values()
    for value in chromosome_class.values()
)

expected_max = max(
    value
    for chromosome_class in EXPECTED.values()
    for value in chromosome_class.values()
)

ymax = max(1.15, observed_max + 0.16, expected_max + 0.16)

for sample, offset in zip(SAMPLES, offsets):

    x = group_centres + offset

    observed = np.array(
        [OBSERVED[c][sample] for c in CLASSES],
        dtype=float,
    )

    expected = np.array(
        [EXPECTED[c][sample] for c in CLASSES],
        dtype=float,
    )

    # Observed bars.
    ax.bar(
        x,
        observed,
        width=bar_width * 0.92,
        color=COLOURS[sample],
        edgecolor="white",
        linewidth=0.6,
        zorder=3,
        label=DISPLAY[sample],
    )

    # Expected markers only.
    ax.scatter(
        x,
        expected,
        s=26,
        facecolor="white",
        edgecolor="#222222",
        linewidth=1.0,
        zorder=6,
    )

    # Observed-expected deviation segments.
    for xi, obs, exp in zip(x, observed, expected):
        ax.plot(
            [xi, xi],
            [obs, exp],
            color="#444444",
            linewidth=0.95,
            zorder=4,
        )


ax.set_xlim(
    group_centres[0] - 0.45,
    group_centres[-1] + 0.45,
)

ax.set_ylim(0, ymax)

ax.set_xticks(group_centres)
ax.set_xticklabels(CLASSES, fontsize=10.5)

ax.set_ylabel("Relative coverage", fontsize=10.5)

ax.grid(
    axis="y",
    color="#E5E5E5",
    linewidth=0.55,
    zorder=0,
)

ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="both", labelsize=9)

fig.suptitle(
    "Observed and expected chromosome coverage",
    fontsize=14,
    fontweight="bold",
    y=0.95,
)

fig.text(
    0.10,
    0.885,
    "Bars show observed relative coverage; open circles show expected dosage.",
    fontsize=9.5,
    color="#555555",
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
    [0], [0],
    color="#222222",
    linestyle="None",
    linewidth=0,
    marker="o",
    markerfacecolor="white",
    markeredgecolor="#222222",
    markersize=5,
    label="Expected dosage",
)

fig.legend(
    handles=sample_handles + [expected_handle],
    loc="lower center",
    bbox_to_anchor=(0.5, 0.03),
    ncol=5,
    frameon=False,
    fontsize=9,
)

pdf_file = RUN / "coverage_observed_vs_expected_grouped.pdf"
png_file = RUN / "coverage_observed_vs_expected_grouped.png"

fig.savefig(pdf_file)
fig.savefig(png_file, dpi=300)

plt.close(fig)

print()
print("VALUES USED")
print(
    summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.6f}",
    )
)
print()
print("WROTE", summary_file.resolve(), sep="\t")
print("WROTE", pdf_file.resolve(), sep="\t")
print("WROTE", png_file.resolve(), sep="\t")
