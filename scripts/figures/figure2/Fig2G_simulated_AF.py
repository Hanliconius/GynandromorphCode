#!/usr/bin/env python3
"""Figure 2G — simulated paired AF distributions in the established founder-model style."""

from __future__ import annotations
import argparse
import ast
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

from common import load_simulation_pair, panel_output, project_root, validate_fraction_arrays

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


def load_historical_model_namespace(root: Path):
    """Load the exact established expected()/predictive_region() definitions."""
    original_script = root / "scripts/plot_simulated_founder_slides.py"
    corrected_script = root / "scripts/plot_founder_slides_empirical_depth.py"
    for path in (original_script, corrected_script):
        if not path.is_file():
            raise FileNotFoundError(path)

    original_source = original_script.read_text()
    original_tree = ast.parse(original_source)
    nodes = []
    for node in original_tree.body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == "rng"
            for target in node.targets
        ):
            break
        if isinstance(node, (ast.Import, ast.ImportFrom, ast.Assign)):
            nodes.append(node)
    nodes.extend(node for node in original_tree.body if isinstance(node, ast.FunctionDef))

    ns = {"__name__": "figure2_model_overlay", "__file__": str(original_script)}
    exec(
        compile(ast.Module(body=nodes, type_ignores=[]), str(original_script), "exec"),
        ns,
    )

    corrected_source = corrected_script.read_text()
    corrected_tree = ast.parse(corrected_source)
    functions = [node for node in corrected_tree.body if isinstance(node, ast.FunctionDef)]
    predictive = [node for node in functions if node.name == "predictive_region"]
    if len(predictive) != 1:
        raise RuntimeError("Expected exactly one corrected predictive_region function")
    predictive_text = ast.get_source_segment(corrected_source, predictive[0]) or ""
    if "histogram2d" in predictive_text or "region_cache" in predictive_text:
        raise RuntimeError("Historical predictive_region is the obsolete binned implementation")

    safe_nodes = [
        node for node in corrected_tree.body
        if isinstance(node, (ast.Import, ast.ImportFrom))
    ]
    for node in corrected_tree.body:
        if isinstance(node, ast.Assign):
            try:
                ast.literal_eval(node.value)
            except (ValueError, TypeError, SyntaxError):
                continue
            safe_nodes.append(node)
    safe_nodes.extend(functions)
    exec(
        compile(ast.Module(body=safe_nodes, type_ignores=[]), str(corrected_script), "exec"),
        ns,
    )

    if "expected" not in ns or "predictive_region" not in ns:
        raise RuntimeError("Historical founder-model drawing functions were not recovered")

    ns["MODEL_COLOUR"] = "#F000B8"
    ns["ns"] = ns
    return ns, corrected_script


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--output", default=None)
    parser.add_argument("--gridsize", type=int, default=40)
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

    model_ns, corrected_model = load_historical_model_namespace(root)

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
    for _, _, left, right in loaded:
        artist = temporary_ax.hexbin(
            left, right,
            gridsize=args.gridsize,
            extent=(0, 1, 0, 1),
            mincnt=1,
        )
        maximum = max(maximum, float(artist.get_array().max()))
        artist.remove()
    plt.close(temporary)
    norm = LogNorm(vmin=1, vmax=maximum)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
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
            q_left, q_right = model_ns["expected"](p_left, p_right, 0.20)
            model_ns["predictive_region"](ax, q_left, q_right)
            ax.plot(
                q_left,
                q_right,
                "o",
                color="#F000B8",
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

    print(f"HISTORICAL_OVERLAY_SOURCE\t{corrected_model}")
    print(f"SITES_PER_MODEL\t{n}")
    print(f"SHARED_MAX_SITES_PER_HEX\t{maximum}")
    print(f"WROTE\t{output}")


if __name__ == "__main__":
    main()
