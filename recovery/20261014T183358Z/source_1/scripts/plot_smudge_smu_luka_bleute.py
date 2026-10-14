#!/usr/bin/env python3

from pathlib import Path
import argparse
import json
import re

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from scipy.ndimage import gaussian_filter


# ============================================================
# Palette
# ============================================================

PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

PALETTE_FILE = (
    PROJECT
    / "palettes"
    / "barfbow"
    / "luka_bleute.json"
)


def is_hex_colour(value):
    return (
        isinstance(value, str)
        and re.fullmatch(
            r"#[0-9A-Fa-f]{6}",
            value,
        )
        is not None
    )


def find_colour_list(obj):
    """
    Recover the first plausible array of hex colours from the
    existing luka_bleute JSON, irrespective of its metadata
    wrapper.
    """

    if isinstance(obj, list):
        colours = [
            value
            for value in obj
            if is_hex_colour(value)
        ]

        if len(colours) >= 3:
            return colours

        for value in obj:
            found = find_colour_list(value)
            if found is not None:
                return found

    if isinstance(obj, dict):

        # Prefer obvious palette keys first.
        for key in (
            "colors",
            "colours",
            "palette",
            "hex",
            "values",
        ):
            if key in obj:
                found = find_colour_list(
                    obj[key]
                )
                if found is not None:
                    return found

        for value in obj.values():
            found = find_colour_list(value)
            if found is not None:
                return found

    return None


if not PALETTE_FILE.is_file():
    raise FileNotFoundError(
        f"Existing luka_bleute palette not found: "
        f"{PALETTE_FILE}"
    )

palette_json = json.loads(
    PALETTE_FILE.read_text()
)

PALETTE = find_colour_list(
    palette_json
)

if PALETTE is None:
    raise RuntimeError(
        "Could not recover hex colours from "
        f"{PALETTE_FILE}"
    )

print(
    "PALETTE",
    PALETTE_FILE,
    f"n={len(PALETTE)}",
    sep="\t",
)


# ============================================================
# Helpers
# ============================================================

def weighted_quantile(values, weights, q):
    order = np.argsort(values)

    values = values[order]
    weights = weights[order]

    cumulative = np.cumsum(weights)

    target = (
        q * cumulative[-1]
    )

    return float(
        values[
            np.searchsorted(
                cumulative,
                target,
            )
        ]
    )


def estimate_haploid_coverage(
    x,
    y,
    weight,
):
    """
    Approximate the underlying 1-copy coverage from the
    dominant balanced AB smudge.

    The AB class has:
       x ~ 0.5
       total coverage ~ 2C
    """

    mask = (
        (x >= 0.43)
        & (x <= 0.50)
        & (y >= 20)
    )

    if mask.sum() < 10:
        return None

    yy = np.rint(
        y[mask]
    ).astype(int)

    ww = weight[mask]

    totals = np.bincount(
        yy,
        weights=ww,
    )

    if not totals.any():
        return None

    # Ignore very lowest residual region.
    lower = 20

    mode_total = (
        np.argmax(
            totals[lower:]
        )
        + lower
    )

    return mode_total / 2.0


# ============================================================
# Arguments
# ============================================================

parser = argparse.ArgumentParser()

parser.add_argument(
    "--input",
    required=True,
)

parser.add_argument(
    "--output-prefix",
    required=True,
)

parser.add_argument(
    "--title",
    default="4 legs combined",
)

parser.add_argument(
    "--ymax",
    type=float,
    default=None,
)

parser.add_argument(
    "--no-labels",
    action="store_true",
)

args = parser.parse_args()


# ============================================================
# Read .smu
#
# Columns are:
#   A k-mer coverage
#   B k-mer coverage
#   number of k-mer pairs
# ============================================================

input_path = Path(
    args.input
)

d = pd.read_csv(
    input_path,
    sep="\t",
    header=None,
    names=[
        "coverage_A",
        "coverage_B",
        "pairs",
    ],
)

for column in d.columns:
    d[column] = pd.to_numeric(
        d[column],
        errors="raise",
    )

d = d.loc[
    (d["coverage_A"] >= 0)
    & (d["coverage_B"] >= 0)
    & (d["pairs"] > 0)
].copy()

A = d["coverage_A"].to_numpy(
    dtype=float
)

B = d["coverage_B"].to_numpy(
    dtype=float
)

weights = d["pairs"].to_numpy(
    dtype=float
)

total = A + B

keep = (
    np.isfinite(total)
    & (total > 0)
)

A = A[keep]
B = B[keep]
weights = weights[keep]
total = total[keep]

minor = np.minimum(
    A,
    B,
)

x = minor / total
y = total


# ============================================================
# Plot limits
# ============================================================

if args.ymax is None:
    ymax = weighted_quantile(
        y,
        weights,
        0.9995,
    )

    ymax = max(
        50.0,
        np.ceil(
            ymax / 10.0
        ) * 10.0,
    )

else:
    ymax = args.ymax


# ============================================================
# Weighted density grid
# ============================================================

X_BINS = 180
Y_BINS = 220

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
    weights=weights,
)

# Display orientation: y by x.
density = density.T

# Mild smoothing reproduces the visual structure of the
# original smudgeplot while retaining the actual weighted data.
density = gaussian_filter(
    density,
    sigma=1.15,
)

positive = density[
    density > 0
]

if not len(positive):
    raise RuntimeError(
        "No positive density after binning."
    )

vmax = float(
    positive.max()
)

# Suppress only the very lowest-density visual noise.
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

x_centres = (
    x_edges[:-1]
    + x_edges[1:]
) / 2

y_centres = (
    y_edges[:-1]
    + y_edges[1:]
) / 2


# ============================================================
# Existing luka_bleute colormap
# ============================================================

cmap = LinearSegmentedColormap.from_list(
    "luka_bleute",
    PALETTE,
    N=256,
)

cmap.set_under(
    "white"
)


# ============================================================
# Figure
# ============================================================

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)

fig = plt.figure(
    figsize=(4.6, 4.35)
)

ax = fig.add_axes(
    [
        0.16,
        0.16,
        0.76,
        0.74,
    ]
)

filled = ax.contourf(
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

ax.set_xlabel(
    "Normalised minor k-mer coverage: B / (A + B)",
    fontsize=9.5,
)

ax.set_ylabel(
    "Total coverage of the k-mer pair: A + B",
    fontsize=9.5,
)

ax.set_title(
    args.title,
    fontsize=11,
    fontweight="bold",
    pad=7,
)

ax.spines[
    "top"
].set_visible(False)

ax.spines[
    "right"
].set_visible(False)


# ============================================================
# Theoretical smudge labels
#
# Infer C from the dominant AB cloud, then place the common
# diploid/polyploid reference compositions at their theoretical
# positions. These labels are reference positions, not fits.
# ============================================================

haploid = estimate_haploid_coverage(
    x,
    y,
    weights,
)

if (
    haploid is not None
    and not args.no_labels
):
    label_specs = [
        ("AB",   1, 1),
        ("AAB",  2, 1),
        ("AAAB", 3, 1),
        ("AABB", 2, 2),
    ]

    for label, major_copies, minor_copies in label_specs:

        copies = (
            major_copies
            + minor_copies
        )

        lx = (
            minor_copies
            / copies
        )

        ly = (
            copies
            * haploid
        )

        if (
            0 < ly < ymax
        ):
            ax.text(
                lx,
                ly,
                label,
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold",
                color="#17212B",
                alpha=0.80,
            )

    print(
        "INFERRED_HAPLOID_COVERAGE",
        f"{haploid:.3f}",
        sep="\t",
    )


# ============================================================
# Compact colourbar
# ============================================================

cax = fig.add_axes(
    [
        0.205,
        0.27,
        0.038,
        0.20,
    ]
)

cbar = fig.colorbar(
    filled,
    cax=cax,
)

cbar.set_label(
    "k-mer pairs",
    fontsize=7.5,
)

cbar.ax.tick_params(
    labelsize=7,
    length=2,
)


# ============================================================
# Vector outputs
# ============================================================

prefix = Path(
    args.output_prefix
)

prefix.parent.mkdir(
    parents=True,
    exist_ok=True,
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
)

fig.savefig(
    svg_file,
)

fig.savefig(
    png_file,
    dpi=300,
)

plt.close(fig)

print(
    "INPUT",
    input_path.resolve(),
    sep="\t",
)

print(
    "ROWS",
    len(d),
    sep="\t",
)

print(
    "TOTAL_KMER_PAIRS",
    int(weights.sum()),
    sep="\t",
)

print(
    "YMAX",
    ymax,
    sep="\t",
)

print(
    "WROTE",
    pdf_file.resolve(),
    sep="\t",
)

print(
    "WROTE",
    svg_file.resolve(),
    sep="\t",
)

print(
    "WROTE",
    png_file.resolve(),
    sep="\t",
)
