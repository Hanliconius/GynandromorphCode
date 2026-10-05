#!/usr/bin/env python3
"""Figure 2E — genotype-stratified AF profiles in the established LOESS style."""

from __future__ import annotations
import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from common import panel_output, project_root

DEFAULT_CACHE = Path(
    "analysis/figures/pairwise_AF_cli/AF_narrative_frames/"
    "original_GT_histograms/ERR2298208_vs_TNLeft_original_GT_counts.npz"
)

COLOURS = ["#2C9BCB", "#F28E2B", "#B61F72", "#AAAAAA"]
LABELS = ["0/0", "0/1", "1/1", "Missing/non-diploid GT"]
CAP = 0.015


def loess(x, y, span=0.20):
    """Original local-linear tricube smoother used in the earlier AF figures."""
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


def draw_profiles(
    ax,
    counts,
    centres,
    *,
    linestyle,
    linewidth,
    alpha,
    zorder,
    show_clipped_labels,
):
    counts = np.asarray(counts, dtype=float)
    n = float(counts.sum())
    if n <= 0:
        raise RuntimeError("Empty histogram")
    fraction = counts / n

    for index, colour in enumerate(COLOURS):
        fitted = loess(centres[1:-1], fraction[index, 1:-1])
        ax.plot(
            centres[1:-1], fitted,
            color=colour,
            alpha=alpha,
            linewidth=linewidth,
            linestyle=linestyle,
            zorder=zorder,
        )

        for endpoint in (0, len(centres) - 1):
            value = float(fraction[index, endpoint])
            if value <= 0:
                continue
            ax.plot(
                [centres[endpoint], centres[endpoint]],
                [0, min(value, CAP)],
                color=colour,
                alpha=alpha,
                linewidth=linewidth,
                linestyle=linestyle,
                zorder=zorder,
            )

        for bin_id in np.flatnonzero(fraction[index] > CAP):
            value = float(fraction[index, bin_id])
            ax.plot(
                centres[bin_id], CAP * 0.98, "^",
                color=colour,
                markersize=5,
                zorder=zorder + 0.2,
                clip_on=False,
            )
            if show_clipped_labels:
                ax.annotate(
                    f"{100 * value:.2f}%",
                    (centres[bin_id], CAP),
                    xytext=((index - 1.5) * 9, 5),
                    textcoords="offset points",
                    rotation=90,
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    color=colour,
                    annotation_clip=False,
                )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--cache", default=str(DEFAULT_CACHE))
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    root = project_root(args.project_root)
    cache_path = Path(args.cache)
    if not cache_path.is_absolute():
        cache_path = root / cache_path
    if not cache_path.is_file():
        raise FileNotFoundError(cache_path)

    output = panel_output(root, "E", args.output)

    with np.load(cache_path, allow_pickle=False) as cache:
        hist = (cache["hist"] if "hist" in cache else cache["counts"]).copy()
        edges = cache["edges"].copy() if "edges" in cache else np.linspace(0.0, 1.0, hist.shape[-1] + 1)

    if hist.shape != (2, 4, 50):
        raise RuntimeError(f"Unexpected AF histogram shape: {hist.shape}")

    centres = (edges[:-1] + edges[1:]) / 2
    control_counts = hist[0]
    left_counts = hist[1]

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })

    fig = plt.figure(figsize=(8.2, 6.2))
    ax = fig.add_axes([0.14, 0.19, 0.78, 0.66])

    draw_profiles(
        ax,
        control_counts,
        centres,
        linestyle=(0, (4, 2)),
        linewidth=1.45,
        alpha=0.43,
        zorder=2,
        show_clipped_labels=False,
    )
    draw_profiles(
        ax,
        left_counts,
        centres,
        linestyle="-",
        linewidth=2.15,
        alpha=0.92,
        zorder=4,
        show_clipped_labels=True,
    )

    ax.axvline(0.5, color="#555555", ls="--", lw=0.8, alpha=0.7, zorder=1)
    ax.set_xlim(-0.015, 1.015)
    ax.set_ylim(0, CAP)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticks([0, 0.005, 0.010, 0.015])
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=1))
    ax.set_xlabel("Observed ALT-read fraction", fontsize=12)
    ax.set_ylabel("Sites per 0.02 AF bin", fontsize=12)
    ax.tick_params(labelsize=10.5)
    ax.grid(axis="y", color="#E5E5E5", lw=0.65)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    genotype_handles = [
        Patch(facecolor=colour, edgecolor="none", label=label)
        for colour, label in zip(COLOURS, LABELS)
    ]
    genotype_legend = ax.legend(
        handles=genotype_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.20),
        frameon=False,
        ncol=4,
        fontsize=9.5,
        title="Original VCF genotype",
        title_fontsize=10,
    )
    ax.add_artist(genotype_legend)

    style_handles = [
        Line2D([0], [0], color="#555555", lw=1.45, alpha=0.50,
               linestyle=(0, (4, 2)), label="Control: ERR2298208"),
        Line2D([0], [0], color="#555555", lw=2.15, alpha=0.95,
               linestyle="-", label="Gynandromorph Left"),
    ]
    ax.legend(handles=style_handles, loc="upper right", frameon=False, fontsize=9.5)

    ax.text(
        -0.12, 1.08, "E",
        transform=ax.transAxes,
        fontsize=17,
        fontweight="bold",
        ha="left",
        va="bottom",
    )
    ax.text(
        0.02, 0.94,
        f"Control n = {int(control_counts.sum()):,}   |   Left n = {int(left_counts.sum()):,}",
        transform=ax.transAxes,
        fontsize=9.5,
        color="#444444",
        ha="left",
        va="top",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.80),
    )

    fig.savefig(output)
    plt.close(fig)

    print(f"INPUT\t{cache_path}")
    print(f"FOCAL_SITES\t{int(left_counts.sum())}")
    print(f"CONTROL_SITES\t{int(control_counts.sum())}")
    print(f"WROTE\t{output}")


if __name__ == "__main__":
    main()
