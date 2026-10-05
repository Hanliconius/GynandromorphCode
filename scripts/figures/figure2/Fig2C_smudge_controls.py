#!/usr/bin/env python3
"""Figure 2C — Smudgeplot controls using the established six-panel aesthetic."""

from __future__ import annotations
import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from common import panel_output, project_root, read_smu, shared_smudge_scale, smudge_density

YMAX = 120.0
ALL_INPUTS = {
    "diploid": Path("analysis/smudgeplot_standard_rerun_20261002/ERR1051012_diploid/L6/kmerpairs.smu"),
    "pseudo": Path("analysis/smudgeplot_lowerL_20261002/pseudoTetraploid_ERR1307111_ERR1307112/L6/kmerpairs.smu"),
    "left": Path("analysis/smudgeplot_arrangements/TN_Left/kmerpairs.smu"),
    "right": Path("analysis/smudgeplot_arrangements/TN_Right/kmerpairs.smu"),
    "t2": Path("analysis/smudgeplot_lowerL_20261002/TN_T2_merge/L6/kmerpairs.smu"),
    "t3": Path("analysis/smudgeplot_lowerL_20261002/TN_T3_merge/L6/kmerpairs.smu"),
}


def prepare(root: Path):
    prepared = {}
    densities = []
    for key, rel in ALL_INPUTS.items():
        path = root / rel
        data = read_smu(path)
        density, x_edges, y_edges = smudge_density(data, ymax=YMAX, x_bins=180, y_bins=220, sigma=1.15)
        prepared[key] = (path, data, density, x_edges, y_edges)
        densities.append(density)
    vmin, vmax, levels = shared_smudge_scale(densities)
    return prepared, vmin, vmax, levels


def draw(ax, density, x_edges, y_edges, *, vmin, vmax, levels, title, subtitle):
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_under("white")
    x_centres = (x_edges[:-1] + x_edges[1:]) / 2
    y_centres = (y_edges[:-1] + y_edges[1:]) / 2
    masked = np.ma.masked_less_equal(density, 0)
    artist = ax.contourf(
        x_centres,
        y_centres,
        masked,
        levels=levels,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        cmap=cmap,
        extend="max",
        antialiased=True,
    )
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, YMAX)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8, length=3)
    ax.text(0.5, 1.055, title, transform=ax.transAxes, fontsize=10.5,
            fontweight="bold", ha="center", va="bottom")
    ax.text(0.5, 1.005, subtitle, transform=ax.transAxes, fontsize=7.4,
            ha="center", va="bottom")
    return artist


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    root = project_root(args.project_root)
    output = panel_output(root, "C", args.output)
    prepared, vmin, vmax, levels = prepare(root)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })

    fig, axes = plt.subplots(2, 1, figsize=(3.7, 6.2), sharex=True, sharey=True)
    fig.text(0.015, 0.99, "C", fontsize=15, fontweight="bold", ha="left", va="top")

    specs = [
        ("diploid", "Diploid control", "ERR1051012  •  L=6  •  AB"),
        ("pseudo", "Pseudo-tetraploid", "ERR1307111 + ERR1307112  •  L=6"),
    ]

    artist = None
    for ax, (key, title, subtitle) in zip(axes, specs):
        path, data, density, x_edges, y_edges = prepared[key]
        artist = draw(
            ax, density, x_edges, y_edges,
            vmin=vmin, vmax=vmax, levels=levels,
            title=title, subtitle=subtitle,
        )
        print("PANEL", key, path, data["rows"], data["total_pairs"], sep="\t")

    axes[-1].set_xlabel("Normalised minor k-mer coverage: B / (A + B)", fontsize=9)
    for ax in axes:
        ax.set_ylabel("Total coverage of the k-mer pair: A + B", fontsize=9)

    cbar = fig.colorbar(artist, ax=axes.tolist(), orientation="vertical", fraction=0.035, pad=0.035)
    cbar.set_label("k-mer pairs", fontsize=8.5)
    cbar.ax.tick_params(labelsize=7.5, length=2)

    fig.subplots_adjust(left=0.22, right=0.84, bottom=0.10, top=0.94, hspace=0.30)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)

    print(f"SHARED_DENSITY_SCALE\t{vmin}\t{vmax}")
    print(f"WROTE\t{output}")


if __name__ == "__main__":
    main()
