#!/usr/bin/env python3
"""Figure 2H — observed paired AF density using the established main-panel aesthetic."""

from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from common import add_historical_predictive_region, informative_density_vmax, panel_output, project_root

BASE = Path("analysis/figures/pairwise_AF_cli")
DEFAULT_BINS = BASE / "TN_merged_control_sites_hex40_batlowW.AF_bins.tsv.gz"
DEFAULT_HIST_ROOT = BASE / "AF_narrative_frames/original_GT_histograms"
DEFAULT_COUNTS = DEFAULT_HIST_ROOT / "original_GT_counts.npz"
DEFAULT_AUDIT = DEFAULT_HIST_ROOT / "original_GT_counts.audit.json"
EXPECTED_N = 31_126_758
GT_COLOURS = ("#2C9BCB", "#F28E2B", "#B61F72")
MODEL_COLOUR = "#E600A9"
CAP = 0.015


def loess(x, y, span=0.20):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    neighbours = min(len(x), max(5, int(np.ceil(span * len(x)))))
    fitted = np.empty_like(y, dtype=float)
    for index, centre in enumerate(x):
        distance = np.abs(x - centre)
        bandwidth = np.partition(distance, neighbours - 1)[neighbours - 1]
        bandwidth = max(float(bandwidth) * 1.000001, 1e-12)
        weights = np.maximum(0.0, 1.0 - (distance / bandwidth) ** 3) ** 3
        design = np.column_stack([np.ones(len(x)), x - centre])
        root_weights = np.sqrt(weights)
        coefficients = np.linalg.lstsq(
            design * root_weights[:, None], y * root_weights, rcond=None
        )[0]
        fitted[index] = coefficients[0]
    return np.maximum(fitted, 0.0)


def marginal(axis, counts, x, limits, vertical=False):
    fraction = counts.astype(float) / counts.sum()
    smooth = loess(x[1:-1], fraction[1:-1])
    colour = "#555555"

    if vertical:
        axis.fill_betweenx(x[1:-1], 0, smooth, color=colour, alpha=0.17)
        axis.plot(smooth, x[1:-1], color=colour, lw=0.9)
        for index in (0, 49):
            value = min(float(fraction[index]), CAP)
            axis.plot([0, value], [x[index], x[index]], color=colour, lw=0.9)
            if fraction[index] > CAP:
                axis.plot(value, x[index], marker="<", ms=3, color=colour, clip_on=False)
                axis.annotate(
                    f"{fraction[index]:.1%}",
                    (value, x[index]),
                    xytext=(-3, 0),
                    textcoords="offset points",
                    ha="right",
                    va="center",
                    fontsize=6,
                    clip_on=False,
                )
        axis.set(ylim=limits, xlim=(CAP * 1.05, 0))
        axis.set_xticks([0, CAP])
        axis.xaxis.set_major_formatter(PercentFormatter(1, decimals=1))
        axis.set_yticks([])
        axis.set_xlabel("Sites / bin", fontsize=6, labelpad=2)
    else:
        axis.fill_between(x[1:-1], 0, smooth, color=colour, alpha=0.17)
        axis.plot(x[1:-1], smooth, color=colour, lw=0.9)
        for index in (0, 49):
            value = min(float(fraction[index]), CAP)
            axis.plot([x[index], x[index]], [0, value], color=colour, lw=0.9)
            if fraction[index] > CAP:
                axis.plot(x[index], value, marker="v", ms=3, color=colour, clip_on=False)
                axis.annotate(
                    f"{fraction[index]:.1%}",
                    (x[index], value),
                    xytext=(0, -4),
                    textcoords="offset points",
                    ha="center",
                    va="top",
                    fontsize=6,
                    clip_on=False,
                )
        axis.set(xlim=limits, ylim=(CAP * 1.05, 0))
        axis.set_yticks([0, CAP])
        axis.yaxis.set_major_formatter(PercentFormatter(1, decimals=1))
        axis.set_xticks([])
        axis.set_ylabel("Sites / bin", fontsize=6, labelpad=2)

    axis.tick_params(labelsize=6, length=2, pad=2)
    for spine in axis.spines.values():
        spine.set_visible(False)


def coloured_marginal(axis, counts_by_genotype, x, limits, vertical=False):
    counts_by_genotype = np.asarray(counts_by_genotype, dtype=float)
    if counts_by_genotype.ndim != 2 or counts_by_genotype.shape[0] < 3:
        raise RuntimeError(
            "Expected genotype-by-AF-bin histogram; "
            f"observed shape {counts_by_genotype.shape}"
        )

    callable_counts = counts_by_genotype[:3]
    total = float(callable_counts.sum())
    if total <= 0:
        raise RuntimeError("Marginal genotype histogram is empty")

    pre_collection_ids = {id(artist) for artist in axis.collections}
    pre_line_ids = {id(artist) for artist in axis.lines}
    marginal(axis, callable_counts.sum(axis=0), x, limits, vertical=vertical)

    for artist in list(axis.collections):
        if id(artist) not in pre_collection_ids:
            artist.remove()
    for artist in list(axis.lines):
        if id(artist) not in pre_line_ids:
            artist.remove()

    upper_limit = axis.get_xlim()[1] if vertical else axis.get_ylim()[1]
    endpoint_marker_level = 0.985 * upper_limit
    endpoint_trigger = 0.92 * upper_limit

    for genotype_index, colour in enumerate(GT_COLOURS):
        fraction = callable_counts[genotype_index] / total
        smooth = loess(x[1:-1], fraction[1:-1])

        if vertical:
            axis.fill_betweenx(
                x[1:-1], 0, smooth,
                color=colour, alpha=0.15, linewidth=0, zorder=2,
            )
            axis.plot(
                smooth, x[1:-1],
                color=colour, linewidth=0.85, alpha=0.95, zorder=3,
            )
        else:
            axis.fill_between(
                x[1:-1], 0, smooth,
                color=colour, alpha=0.15, linewidth=0, zorder=2,
            )
            axis.plot(
                x[1:-1], smooth,
                color=colour, linewidth=0.85, alpha=0.95, zorder=3,
            )

        endpoint0 = float(fraction[0])
        endpoint1 = float(fraction[-1])
        if vertical:
            if endpoint0 >= endpoint_trigger:
                axis.scatter(endpoint_marker_level, 0.0, marker=">", s=22,
                             color=colour, linewidths=0, clip_on=False, zorder=5)
            if endpoint1 >= endpoint_trigger:
                axis.scatter(endpoint_marker_level, 1.0, marker=">", s=22,
                             color=colour, linewidths=0, clip_on=False, zorder=5)
        else:
            if endpoint0 >= endpoint_trigger:
                axis.scatter(0.0, endpoint_marker_level, marker="^", s=22,
                             color=colour, linewidths=0, clip_on=False, zorder=5)
            if endpoint1 >= endpoint_trigger:
                axis.scatter(1.0, endpoint_marker_level, marker="^", s=22,
                             color=colour, linewidths=0, clip_on=False, zorder=5)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--bins", default=str(DEFAULT_BINS))
    parser.add_argument("--counts", default=str(DEFAULT_COUNTS))
    parser.add_argument("--audit", default=str(DEFAULT_AUDIT))
    parser.add_argument("--output", default=None)
    parser.add_argument("--density-quantile", type=float, default=0.995)
    parser.add_argument("--origin-exclude", type=float, default=0.05)
    parser.add_argument("--density-blend", type=float, default=0.50, help="0 = fully compressed informative scale; 1 = original raw-max scale")
    args = parser.parse_args()

    root = project_root(args.project_root)
    paths = [Path(args.bins), Path(args.counts), Path(args.audit)]
    paths = [path if path.is_absolute() else root / path for path in paths]
    bins_path, counts_path, audit_path = paths
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)

    output = panel_output(root, "H", args.output)
    palette_path = root / "palettes/barfbow/luka_bleute.json"
    if not palette_path.is_file():
        raise FileNotFoundError(palette_path)
    colours = json.loads(palette_path.read_text())
    if not isinstance(colours, list) or len(colours) < 2:
        raise RuntimeError(f"Invalid Luka-Bleute palette: {palette_path}")
    cmap = ListedColormap(colours)

    bins = pd.read_csv(bins_path, sep="\t")
    required = {"left_centre", "right_centre", "width", "height", "count"}
    missing = required - set(bins.columns)
    if missing:
        raise RuntimeError(f"Hex table missing columns: {sorted(missing)}")

    with np.load(counts_path, allow_pickle=False) as cache:
        hist = cache["hist"].copy()
    metadata = json.loads(audit_path.read_text())
    samples = metadata["signature"]["samples"]
    if hist.shape != (3, 4, 50):
        raise RuntimeError(f"Unexpected histogram shape: {hist.shape}")

    left_index = samples.index("TN_Left_merged")
    right_index = samples.index("TN_Right_merged")
    left_counts = hist[left_index]
    right_counts = hist[right_index]
    n = int(bins["count"].sum())
    if n != int(left_counts.sum()) or n != int(right_counts.sum()):
        raise RuntimeError(
            "Hex plot and marginal histograms have different site totals: "
            f"hex={n}, left={left_counts.sum()}, right={right_counts.sum()}"
        )
    if n != EXPECTED_N:
        raise RuntimeError(f"Expected {EXPECTED_N:,} retained sites; found {n:,}")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 7,
        "axes.labelsize": 8,
        "axes.linewidth": 0.6,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "text.usetex": False,
        "savefig.facecolor": "white",
    })

    fig = plt.figure(figsize=(5.3, 5.3))
    fig.text(0.015, 0.992, "H", fontsize=15, fontweight="bold", ha="left", va="top")
    ax = fig.add_axes([0.275, 0.275, 0.60, 0.60])
    top = fig.add_axes([0.275, 0.885, 0.60, 0.065])
    side = fig.add_axes([0.120, 0.275, 0.065, 0.60])
    cax = fig.add_axes([0.905, 0.405, 0.022, 0.34])

    occupied = bins.loc[bins["count"] > 0].copy()
    unit_hex = np.array([
        [0.0, -0.5],
        [0.5, -0.25],
        [0.5, 0.25],
        [0.0, 0.5],
        [-0.5, 0.25],
        [-0.5, -0.25],
    ])
    centres = occupied[["left_centre", "right_centre"]].to_numpy()
    sizes = occupied[["width", "height"]].to_numpy()
    vertices = centres[:, None, :] + unit_hex[None, :, :] * sizes[:, None, :]

    counts = occupied["count"].to_numpy(dtype=float)
    compressed_vmax = informative_density_vmax(
        counts,
        occupied["left_centre"].to_numpy(dtype=float),
        occupied["right_centre"].to_numpy(dtype=float),
        origin_exclude=args.origin_exclude,
        quantile=args.density_quantile,
    )
    raw_max = float(counts.max())
    if not (0.0 <= args.density_blend <= 1.0):
        raise ValueError("--density-blend must be between 0 and 1")
    # Interpolate in log space because the colour normalization itself is log.
    # 0.0 reproduces the aggressively compressed scale; 1.0 reproduces the
    # original raw-maximum scale. 0.5 is the exact visual midpoint in log space.
    display_vmax = float(np.exp(
        (1.0 - args.density_blend) * np.log(compressed_vmax)
        + args.density_blend * np.log(raw_max)
    ))

    collection = PolyCollection(
        vertices,
        array=counts,
        cmap=cmap,
        norm=LogNorm(vmin=1, vmax=display_vmax, clip=True),
        edgecolors="none",
        linewidths=0,
        rasterized=False,
    )
    ax.add_collection(collection)

    padding = max(float(bins["width"].max()), float(bins["height"].max())) / 2 + 0.003
    limits = (-padding, 1 + padding)
    ax.set(xlim=limits, ylim=limits, aspect="equal")
    ax.plot([0, 1], [0, 1], color="white", lw=0.5, alpha=0.45)

    mixing = 0.20
    error = 0.001
    for p_left in (0.0, 0.5, 1.0):
        for p_right in (0.0, 0.5, 1.0):
            q_left = error + (1 - 2 * error) * (
                (1 - mixing) * p_left + mixing * p_right
            )
            q_right = error + (1 - 2 * error) * (
                mixing * p_left + (1 - mixing) * p_right
            )

            add_historical_predictive_region(
                ax,
                q_left,
                q_right,
                depth=40,
                concentration=22.0,
                mass=0.50,
                colour=MODEL_COLOUR,
                linewidth=0.70,
                alpha=0.80,
            )

            ax.scatter(
                [q_left], [q_right],
                s=14,
                color=MODEL_COLOUR,
                alpha=0.68,
                edgecolors="white",
                linewidths=0.45,
                zorder=10,
            )

    ticks = [0, 0.25, 0.5, 0.75, 1]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xlabel("Left allele fraction", labelpad=3)
    ax.set_ylabel("Right allele fraction", labelpad=3)
    ax.tick_params(length=2.5, pad=2)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    cb = fig.colorbar(collection, cax=cax, extend="max")
    cb.set_label("Sites per hexagon", fontsize=7, labelpad=5)
    cb.ax.tick_params(labelsize=6, length=2)
    cb.outline.set_linewidth(0.5)

    x = (np.arange(50) + 0.5) / 50
    coloured_marginal(top, left_counts, x, limits)
    coloured_marginal(side, right_counts, x, limits, vertical=True)

    fig.text(0.275, 0.905, f"{n:,} sites", fontsize=7)
    fig.text(
        0.275, 0.035,
        "20% mixing model  •  Centroids and 50% predictive regions",
        fontsize=6.5,
    )

    fig.savefig(output, dpi=600)
    plt.close(fig)

    print(f"BINS\t{bins_path}")
    print(f"COUNTS\t{counts_path}")
    print("HISTORICAL_OVERLAY\tbeta-binomial depth=40 concentration=22 mass=0.50")
    print(f"SITES\t{n}")
    print(f"RAW_MAX_SITES_PER_HEX\t{raw_max}")
    print(f"COMPRESSED_VMAX_SITES_PER_HEX\t{compressed_vmax}")
    print(f"DISPLAY_VMAX_SITES_PER_HEX\t{display_vmax}")
    print(f"DENSITY_SCALE\tLogNorm log-space blend={args.density_blend}; compressed_q={args.density_quantile}; origin_exclude={args.origin_exclude}")
    print(f"WROTE\t{output}")


if __name__ == "__main__":
    main()
