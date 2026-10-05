#!/usr/bin/env python3
"""Figure 2G — simulated paired AF distributions in the established founder-model style."""

from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm
import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from common import MODEL_COLOUR, add_historical_predictive_region, expected_mixed_af, informative_density_vmax, load_simulation_pair, panel_output, project_root, validate_fraction_arrays

NARRATIVE = Path("analysis/figures/pairwise_AF_cli/AF_narrative_frames")
SCENARIOS = [
    {
        "name": "AB / A",
        "subtitle": "Diploid / haploid",
        "source": NARRATIVE / "AB_A_digametic/AB_A_simulated_AF.npz",
        "prefix": "mix20",
        "states": [(0., 0.), (.5, 0.), (.5, 1.), (1., 1.)],
    },
    {
        "name": "AA / BC",
        "subtitle": "Duplicated haplotype",
        "source": NARRATIVE / "gyno_three_founder/gyno_simulated_AF.npz",
        "prefix": "mix20",
        "states": [(l, r) for l in (0., 1.) for r in (0., .5, 1.)],
    },
    {
        "name": "AB / AC",
        "subtitle": "Shared haplotype",
        "source": NARRATIVE / "simulated_founder_empirical_depth/simulated_AF.npz",
        "prefix": "mix20_model3",
        "states": [
            (l, r)
            for l in (0., .5, 1.)
            for r in (0., .5, 1.)
            if (l, r) not in {(0., 1.), (1., 0.)}
        ],
    },
    {
        "name": "AB / CD",
        "subtitle": "Four haplotypes",
        "source": NARRATIVE / "simulated_founder_empirical_depth/simulated_AF.npz",
        "prefix": "mix20_model4",
        "states": [(l, r) for l in (0., .5, 1.) for r in (0., .5, 1.)],
    },
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--output", default=None)
    parser.add_argument("--gridsize", type=int, default=40)
    parser.add_argument("--density-quantile", type=float, default=0.995)
    parser.add_argument("--origin-exclude", type=float, default=0.05)
    parser.add_argument("--density-blend", type=float, default=0.50, help="0 = fully compressed informative scale; 1 = original raw-max scale")
    args = parser.parse_args()

    root = project_root(args.project_root)
    output = panel_output(root, "G", args.output)

    palette_path = root / "palettes/barfbow/luka_bleute.json"
    if not palette_path.is_file():
        raise FileNotFoundError(palette_path)
    palette = json.loads(palette_path.read_text())
    if not isinstance(palette, list) or len(palette) < 2:
        raise RuntimeError(f"Invalid Luka-Bleute palette: {palette_path}")
    cmap = ListedColormap(palette, name="luka_bleute")

    loaded = []
    sizes = set()
    for scenario in SCENARIOS:
        path = root / scenario["source"]
        left, right = load_simulation_pair(path, scenario["prefix"])
        validate_fraction_arrays(left, right)
        loaded.append((scenario, path, left, right))
        sizes.add(len(left))
    if len(sizes) != 1:
        raise RuntimeError(f"Simulation panel sizes differ: {sorted(sizes)}")
    n = sizes.pop()

    temporary, temporary_ax = plt.subplots()
    maximum = 1.0
    scale_counts = []
    scale_x = []
    scale_y = []
    for _, _, left, right in loaded:
        artist = temporary_ax.hexbin(
            left, right,
            gridsize=args.gridsize,
            extent=(0, 1, 0, 1),
            mincnt=1,
        )
        counts = np.asarray(artist.get_array(), dtype=float)
        offsets = np.asarray(artist.get_offsets(), dtype=float)
        maximum = max(maximum, float(counts.max()))
        scale_counts.append(counts)
        scale_x.append(offsets[:, 0])
        scale_y.append(offsets[:, 1])
        artist.remove()
    plt.close(temporary)

    scale_counts = np.concatenate(scale_counts)
    scale_x = np.concatenate(scale_x)
    scale_y = np.concatenate(scale_y)
    compressed_vmax = informative_density_vmax(
        scale_counts,
        scale_x,
        scale_y,
        origin_exclude=args.origin_exclude,
        quantile=args.density_quantile,
    )
    if not (0.0 <= args.density_blend <= 1.0):
        raise ValueError("--density-blend must be between 0 and 1")
    # Interpolate in log space because the colour normalization itself is log.
    # 0.0 reproduces the aggressively compressed scale; 1.0 reproduces the
    # original raw-maximum scale. 0.5 is the exact visual midpoint in log space.
    display_vmax = float(np.exp(
        (1.0 - args.density_blend) * np.log(compressed_vmax)
        + args.density_blend * np.log(maximum)
    ))
    norm = LogNorm(vmin=1, vmax=display_vmax, clip=True)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "text.usetex": False,
        "savefig.facecolor": "white",
    })

    fig, axes = plt.subplots(2, 2, figsize=(7.15, 7.0), sharex=True, sharey=True)
    fig.text(0.015, 0.992, "G", fontsize=15, fontweight="bold", ha="left", va="top")

    last_artist = None
    for ax, (scenario, path, left, right) in zip(axes.ravel(), loaded):
        last_artist = ax.hexbin(
            left,
            right,
            gridsize=args.gridsize,
            extent=(0, 1, 0, 1),
            mincnt=1,
            cmap=cmap,
            norm=norm,
            linewidths=0,
            rasterized=False,
        )
        ax.set_xlim(-0.018, 1.018)
        ax.set_ylim(-0.022, 1.022)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([0, .5, 1])
        ax.set_yticks([0, .5, 1])
        ax.tick_params(labelsize=9, length=3)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.plot([0, 1], [0, 1], color="white", alpha=.35, lw=.7)

        ax.text(
            0.5, 1.075, scenario["name"],
            transform=ax.transAxes,
            fontsize=11,
            fontweight="bold",
            ha="center",
            va="bottom",
        )
        ax.text(
            0.5, 1.025, scenario["subtitle"],
            transform=ax.transAxes,
            fontsize=8,
            ha="center",
            va="bottom",
        )

        for p_left, p_right in scenario["states"]:
            q_left, q_right = expected_mixed_af(
                p_left, p_right, c_left=0.20, c_right=0.20, error=0.001
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
            ax.plot(
                q_left,
                q_right,
                "o",
                color=MODEL_COLOUR,
                markersize=4,
                markeredgecolor="white",
                markeredgewidth=.6,
                zorder=8,
                clip_on=False,
            )

        print(f"SOURCE\t{scenario['name']}\t{path}\t{len(left)}")

    for ax in axes[-1, :]:
        ax.set_xlabel("Left ALT-read fraction", fontsize=10)
    for ax in axes[:, 0]:
        ax.set_ylabel("Right ALT-read fraction", fontsize=10)

    cbar = fig.colorbar(
        last_artist,
        ax=axes.ravel().tolist(),
        orientation="vertical",
        fraction=0.025,
        pad=0.025,
        extend="max",
    )
    cbar.ax.tick_params(labelsize=8, length=2)
    cbar.set_label("Sites / hexagon", fontsize=9)

    fig.text(
        0.11,
        0.025,
        "20% read mixing  •  dots: predicted centroids  •  outlines: analytical 50% predictive regions",
        fontsize=8,
        color="#555555",
        ha="left",
    )

    fig.subplots_adjust(left=0.11, right=0.88, bottom=0.09, top=0.93, hspace=0.27, wspace=0.13)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)

    print("HISTORICAL_OVERLAY\tbeta-binomial depth=40 concentration=22 mass=0.50")
    print(f"SITES_PER_MODEL\t{n}")
    print(f"RAW_MAX_SITES_PER_HEX\t{maximum}")
    print(f"COMPRESSED_VMAX_SITES_PER_HEX\t{compressed_vmax}")
    print(f"DISPLAY_VMAX_SITES_PER_HEX\t{display_vmax}")
    print(f"DENSITY_SCALE\tLogNorm log-space blend={args.density_blend}; compressed_q={args.density_quantile}; origin_exclude={args.origin_exclude}")
    print(f"WROTE\t{output}")


if __name__ == "__main__":
    main()
