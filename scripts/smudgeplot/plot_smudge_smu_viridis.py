#!/usr/bin/env python3

from pathlib import Path
import argparse

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter


def weighted_quantile(values, weights, q):
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    target = q * cumulative[-1]
    return float(values[np.searchsorted(cumulative, target)])


def estimate_haploid_coverage(x, y, weight):
    """
    Approximate the underlying one-copy coverage from the dominant
    balanced AB-like cloud.

    This is a plotting aid only, not a ploidy fit.
    """
    mask = (x >= 0.43) & (x <= 0.50) & (y >= 20)

    if mask.sum() < 10:
        return None

    yy = np.rint(y[mask]).astype(int)
    ww = weight[mask]

    totals = np.bincount(yy, weights=ww)

    if not totals.any():
        return None

    lower = 20
    mode_total = np.argmax(totals[lower:]) + lower

    return mode_total / 2.0


parser = argparse.ArgumentParser()

parser.add_argument("--input", required=True)
parser.add_argument("--output-prefix", required=True)
parser.add_argument("--title", default="Smudgeplot density")
parser.add_argument("--ymax", type=float, default=None)
parser.add_argument("--no-labels", action="store_true")

args = parser.parse_args()


input_path = Path(args.input)

d = pd.read_csv(
    input_path,
    sep="\t",
    header=None,
    names=["coverage_A", "coverage_B", "pairs"],
)

for column in d.columns:
    d[column] = pd.to_numeric(d[column], errors="raise")

d = d.loc[
    (d["coverage_A"] >= 0)
    & (d["coverage_B"] >= 0)
    & (d["pairs"] > 0)
].copy()

A = d["coverage_A"].to_numpy(dtype=float)
B = d["coverage_B"].to_numpy(dtype=float)
weights = d["pairs"].to_numpy(dtype=float)

total = A + B

keep = np.isfinite(total) & (total > 0)

A = A[keep]
B = B[keep]
weights = weights[keep]
total = total[keep]

minor = np.minimum(A, B)

x = minor / total
y = total


if args.ymax is None:
    ymax = weighted_quantile(y, weights, 0.9995)
    ymax = max(50.0, np.ceil(ymax / 10.0) * 10.0)
else:
    ymax = args.ymax


X_BINS = 180
Y_BINS = 220

density, x_edges, y_edges = np.histogram2d(
    x,
    y,
    bins=[X_BINS, Y_BINS],
    range=[[0.0, 0.5], [0.0, ymax]],
    weights=weights,
)

density = density.T
density = gaussian_filter(density, sigma=1.15)

positive = density[density > 0]

if not len(positive):
    raise RuntimeError("No positive density after binning.")

vmax = float(positive.max())
vmin = max(1.0, float(np.quantile(positive, 0.04)))

if vmin >= vmax:
    vmin = max(1.0, vmax / 100.0)

levels = np.geomspace(vmin, vmax, 32)

x_centres = (x_edges[:-1] + x_edges[1:]) / 2
y_centres = (y_edges[:-1] + y_edges[1:]) / 2

cmap = plt.get_cmap("viridis").copy()
cmap.set_under("white")


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

fig = plt.figure(figsize=(4.6, 4.35))
ax = fig.add_axes([0.16, 0.16, 0.76, 0.74])

filled = ax.contourf(
    x_centres,
    y_centres,
    density,
    levels=levels,
    norm=LogNorm(vmin=vmin, vmax=vmax),
    cmap=cmap,
    extend="max",
    antialiased=True,
)

ax.set_xlim(0, 0.5)
ax.set_ylim(0, ymax)

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

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)


haploid = estimate_haploid_coverage(x, y, weights)

if haploid is not None and not args.no_labels:
    label_specs = [
        ("AB", 1, 1),
        ("AAB", 2, 1),
        ("AAAB", 3, 1),
        ("AABB", 2, 2),
    ]

    for label, major_copies, minor_copies in label_specs:
        copies = major_copies + minor_copies
        lx = minor_copies / copies
        ly = copies * haploid

        if 0 < ly < ymax:
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


cax = fig.add_axes([0.205, 0.27, 0.038, 0.20])

cbar = fig.colorbar(filled, cax=cax)

cbar.set_label(
    "k-mer pairs",
    fontsize=7.5,
)

cbar.ax.tick_params(
    labelsize=7,
    length=2,
)


prefix = Path(args.output_prefix)
prefix.parent.mkdir(parents=True, exist_ok=True)

pdf_file = prefix.with_suffix(".pdf")
svg_file = prefix.with_suffix(".svg")
png_file = prefix.with_suffix(".png")

fig.savefig(pdf_file)
fig.savefig(svg_file)
fig.savefig(png_file, dpi=300)

plt.close(fig)

print("INPUT", input_path.resolve(), sep="\t")
print("ROWS", len(d), sep="\t")
print("TOTAL_KMER_PAIRS", int(weights.sum()), sep="\t")
print("YMAX", ymax, sep="\t")
print("WROTE", pdf_file.resolve(), sep="\t")
print("WROTE", svg_file.resolve(), sep="\t")
print("WROTE", png_file.resolve(), sep="\t")
