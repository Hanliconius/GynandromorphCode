#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter


ROOT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

OUTDIR = ROOT / "analysis/smudgeplot_4legs_combined"
OUTDIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# Inputs
# ============================================================

PANELS = [
    {
        "title": "Diploid control",
        "subtitle": "ERR1051012  •  L=4  •  1n=8.01×  •  AB",
        "path": ROOT / "analysis/smudgeplot_arrangements/ERR1051012_krista/kmerpairs_L4.smu",
    },
    {
        "title": "Pseudo-tetraploid",
        "subtitle": "ERR1307111 + ERR1307112  •  1n=24.85×",
        "path": ROOT / "analysis/smudgeplot_arrangements/pseudoTetraploid_ERR1307111_ERR1307112/kmerpairs.smu",
    },
    {
        "title": "Left",
        "subtitle": "TN Left: T2 + T3",
        "path": ROOT / "analysis/smudgeplot_arrangements/TN_Left/kmerpairs.smu",
    },
    {
        "title": "Right",
        "subtitle": "TN Right: T2 + T3",
        "path": ROOT / "analysis/smudgeplot_arrangements/TN_Right/kmerpairs.smu",
    },
    {
        "title": "T2 merge",
        "subtitle": "TN Left T2 + Right T2",
        "path": ROOT / "analysis/smudgeplot_arrangements/TN_T2_merge/kmerpairs.smu",
    },
    {
        "title": "T3 merge",
        "subtitle": "TN Left T3 + Right T3",
        "path": ROOT / "analysis/smudgeplot_arrangements/TN_T3_merge/kmerpairs.smu",
    },
]


# ============================================================
# Helpers
# ============================================================

def weighted_quantile(values, weights, q):
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]

    cumulative = np.cumsum(weights)

    return float(
        values[
            np.searchsorted(
                cumulative,
                q * cumulative[-1],
            )
        ]
    )


def read_smu(path):
    if not path.is_file():
        raise FileNotFoundError(path)

    d = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "coverage_A",
            "coverage_B",
            "pairs",
        ],
    )

    for c in d.columns:
        d[c] = pd.to_numeric(
            d[c],
            errors="raise",
        )

    d = d.loc[
        (d["coverage_A"] >= 0)
        & (d["coverage_B"] >= 0)
        & (d["pairs"] > 0)
    ].copy()

    A = d["coverage_A"].to_numpy(dtype=float)
    B = d["coverage_B"].to_numpy(dtype=float)
    weights = d["pairs"].to_numpy(dtype=float)

    total = A + B

    keep = (
        np.isfinite(total)
        & np.isfinite(weights)
        & (total > 0)
        & (weights > 0)
    )

    A = A[keep]
    B = B[keep]
    weights = weights[keep]
    total = total[keep]

    minor = np.minimum(A, B)

    x = minor / total
    y = total

    return {
        "x": x,
        "y": y,
        "weights": weights,
        "rows": len(d),
        "total_pairs": int(weights.sum()),
    }


# ============================================================
# Read all inputs first
# ============================================================

for panel in PANELS:
    panel["data"] = read_smu(
        panel["path"]
    )


# ============================================================
# Common y axis
#
# Use global weighted 99.95th percentile, rounded to 10.
# This gives every panel exactly the same coordinate system.
# ============================================================

all_y = np.concatenate(
    [
        p["data"]["y"]
        for p in PANELS
    ]
)

all_w = np.concatenate(
    [
        p["data"]["weights"]
        for p in PANELS
    ]
)

ymax = weighted_quantile(
    all_y,
    all_w,
    0.9995,
)

ymax = max(
    50.0,
    np.ceil(ymax / 10.0) * 10.0,
)

print(
    "COMMON_YMAX",
    ymax,
    sep="\t",
)


# ============================================================
# Density grids
# ============================================================

X_BINS = 180
Y_BINS = 220

all_positive = []

for panel in PANELS:

    x = panel["data"]["x"]
    y = panel["data"]["y"]
    w = panel["data"]["weights"]

    density, x_edges, y_edges = np.histogram2d(
        x,
        y,
        bins=[
            X_BINS,
            Y_BINS,
        ],
        range=[
            [0.0, 0.5],
            [0.0, ymax],
        ],
        weights=w,
    )

    density = density.T

    density = gaussian_filter(
        density,
        sigma=1.15,
    )

    panel["density"] = density
    panel["x_edges"] = x_edges
    panel["y_edges"] = y_edges

    positive = density[
        density > 0
    ]

    if len(positive):
        all_positive.append(
            positive
        )


if not all_positive:
    raise RuntimeError(
        "No positive density in any panel."
    )


# ============================================================
# Shared density scale
# ============================================================

positive = np.concatenate(
    all_positive
)

vmax = float(
    positive.max()
)

vmin = max(
    1.0,
    float(
        np.quantile(
            positive,
            0.04,
        )
    ),
)

if vmin >= vmax:
    vmin = max(
        1.0,
        vmax / 100.0,
    )

levels = np.geomspace(
    vmin,
    vmax,
    32,
)

print(
    "GLOBAL_VMIN",
    vmin,
    sep="\t",
)

print(
    "GLOBAL_VMAX",
    vmax,
    sep="\t",
)


# ============================================================
# Style
# ============================================================

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)

cmap = plt.get_cmap(
    "viridis"
).copy()

cmap.set_under(
    "white"
)


# ============================================================
# Figure
# ============================================================

fig, axes = plt.subplots(
    3,
    2,
    figsize=(7.2, 9.0),
    sharex=True,
    sharey=True,
)

axes = axes.ravel()

last_filled = None

panel_letters = [
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
]


for ax, panel, letter in zip(
    axes,
    PANELS,
    panel_letters,
):

    density = panel["density"]

    x_edges = panel["x_edges"]
    y_edges = panel["y_edges"]

    x_centres = (
        x_edges[:-1]
        + x_edges[1:]
    ) / 2

    y_centres = (
        y_edges[:-1]
        + y_edges[1:]
    ) / 2

    last_filled = ax.contourf(
        x_centres,
        y_centres,
        density,
        levels=levels,
        norm=LogNorm(
            vmin=vmin,
            vmax=vmax,
        ),
        cmap=cmap,
        extend="max",
        antialiased=True,
    )

    ax.set_xlim(
        0,
        0.5,
    )

    ax.set_ylim(
        0,
        ymax,
    )

    ax.spines[
        "top"
    ].set_visible(False)

    ax.spines[
        "right"
    ].set_visible(False)

    ax.tick_params(
        labelsize=8,
        length=3,
    )

    ax.text(
        -0.11,
        1.04,
        letter,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        ha="left",
        va="bottom",
    )

    ax.text(
        0.5,
        1.055,
        panel["title"],
        transform=ax.transAxes,
        fontsize=10.5,
        fontweight="bold",
        ha="center",
        va="bottom",
    )

    ax.text(
        0.5,
        1.005,
        panel["subtitle"],
        transform=ax.transAxes,
        fontsize=7.4,
        ha="center",
        va="bottom",
    )

    print(
        "PANEL",
        letter,
        panel["title"],
        panel["path"],
        panel["data"]["rows"],
        panel["data"]["total_pairs"],
        sep="\t",
    )


# ============================================================
# Shared labels
# ============================================================

fig.supxlabel(
    "Normalised minor k-mer coverage: B / (A + B)",
    fontsize=10,
    y=0.045,
)

fig.supylabel(
    "Total coverage of the k-mer pair: A + B",
    fontsize=10,
    x=0.035,
)


# ============================================================
# Shared colourbar
# ============================================================

cbar = fig.colorbar(
    last_filled,
    ax=axes.tolist(),
    orientation="vertical",
    fraction=0.025,
    pad=0.025,
)

cbar.set_label(
    "k-mer pairs",
    fontsize=8.5,
)

cbar.ax.tick_params(
    labelsize=7.5,
    length=2,
)


# ============================================================
# Finish
# ============================================================

fig.subplots_adjust(
    left=0.11,
    right=0.88,
    bottom=0.09,
    top=0.95,
    hspace=0.27,
    wspace=0.13,
)

prefix = (
    OUTDIR
    / "gynandromorph_controls_sixpanel_viridis"
)

pdf_file = prefix.with_suffix(
    ".pdf"
)

svg_file = prefix.with_suffix(
    ".svg"
)

png_file = prefix.with_suffix(
    ".png"
)

fig.savefig(
    pdf_file,
    bbox_inches="tight",
)

fig.savefig(
    svg_file,
    bbox_inches="tight",
)

fig.savefig(
    png_file,
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)

print(
    "WROTE",
    pdf_file,
    sep="\t",
)

print(
    "WROTE",
    svg_file,
    sep="\t",
)

print(
    "WROTE",
    png_file,
    sep="\t",
)
