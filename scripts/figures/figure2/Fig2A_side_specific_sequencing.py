#!/usr/bin/env python3
"""Figure 2A — side-specific sequencing schematic."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import Circle, PathPatch

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from common import (
    LEFT_T2,
    LEFT_T3,
    RIGHT_T2,
    RIGHT_T3,
    add_panel_letter,
    panel_output,
    project_root,
    set_publication_style,
)


def wing_patch(points, facecolor, edgecolor="black", linewidth=0.8, alpha=1.0, zorder=1):
    verts = [points[0]]
    codes = [MplPath.MOVETO]
    for control1, control2, endpoint in points[1:]:
        verts.extend([control1, control2, endpoint])
        codes.extend([MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
    verts.append(points[0])
    codes.append(MplPath.CLOSEPOLY)
    return PathPatch(
        MplPath(verts, codes),
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        alpha=alpha,
        zorder=zorder,
    )


def add_butterfly(ax):
    # Stylised vector butterfly: intended as a clean schematic, not a taxonomic illustration.
    left_fore = [
        ((0.49, 0.64)),
        ((0.40, 0.82), (0.22, 0.92), (0.10, 0.82)),
        ((0.03, 0.74), (0.10, 0.57), (0.28, 0.55)),
        ((0.36, 0.54), (0.44, 0.58), (0.49, 0.64)),
    ]
    left_hind = [
        ((0.48, 0.60)),
        ((0.34, 0.57), (0.14, 0.53), (0.12, 0.36)),
        ((0.11, 0.23), (0.28, 0.20), (0.40, 0.34)),
        ((0.47, 0.43), (0.49, 0.52), (0.48, 0.60)),
    ]
    right_fore = [
        ((0.51, 0.64)),
        ((0.60, 0.82), (0.78, 0.92), (0.90, 0.82)),
        ((0.97, 0.74), (0.90, 0.57), (0.72, 0.55)),
        ((0.64, 0.54), (0.56, 0.58), (0.51, 0.64)),
    ]
    right_hind = [
        ((0.52, 0.60)),
        ((0.66, 0.57), (0.86, 0.53), (0.88, 0.36)),
        ((0.89, 0.23), (0.72, 0.20), (0.60, 0.34)),
        ((0.53, 0.43), (0.51, 0.52), (0.52, 0.60)),
    ]

    for spec in (left_fore, left_hind, right_fore, right_hind):
        ax.add_patch(wing_patch(spec, "#151515", linewidth=0.7, zorder=1))

    # Red/orange postman bands.
    left_band = [
        ((0.14, 0.46)),
        ((0.20, 0.42), (0.33, 0.40), (0.42, 0.47)),
        ((0.39, 0.39), (0.29, 0.30), (0.18, 0.31)),
        ((0.13, 0.33), (0.12, 0.39), (0.14, 0.46)),
    ]
    right_band = [
        ((0.86, 0.46)),
        ((0.80, 0.42), (0.67, 0.40), (0.58, 0.47)),
        ((0.61, 0.39), (0.71, 0.30), (0.82, 0.31)),
        ((0.87, 0.33), (0.88, 0.39), (0.86, 0.46)),
    ]
    for spec in (left_band, right_band):
        ax.add_patch(wing_patch(spec, "#D94E2C", edgecolor="none", zorder=2))

    # Pale forewing patches.
    for x in (0.28, 0.72):
        ax.add_patch(Circle((x, 0.72), 0.055, facecolor="#F1EEE6", edgecolor="none", zorder=3))

    # Body/head/antennae.
    ax.plot([0.5, 0.5], [0.35, 0.73], color="#111111", linewidth=5.5, solid_capstyle="round", zorder=4)
    ax.add_patch(Circle((0.5, 0.755), 0.025, facecolor="#111111", edgecolor="none", zorder=4))
    ax.plot([0.494, 0.45], [0.77, 0.84], color="#111111", linewidth=0.8, zorder=4)
    ax.plot([0.506, 0.55], [0.77, 0.84], color="#111111", linewidth=0.8, zorder=4)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    root = project_root(args.project_root)
    output = panel_output(root, "A", args.output)
    set_publication_style(8.0)

    fig = plt.figure(figsize=(3.35, 4.35))
    add_panel_letter(fig, "A", x=0.015, y=0.985, size=14)
    fig.text(0.17, 0.965, "side-specific sequencing", ha="left", va="top", fontsize=8.5)

    ax = fig.add_axes([0.12, 0.49, 0.76, 0.42])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    add_butterfly(ax)
    ax.plot([0.5, 0.5], [0.12, 0.93], color="#666666", linewidth=0.8, linestyle=(0, (2, 2)), zorder=5)
    ax.text(0.16, 0.12, "♀", fontsize=17, ha="center", va="center")
    ax.text(0.84, 0.12, "♂", fontsize=17, ha="center", va="center")

    # Four sequencing libraries beneath the specimen.
    ax2 = fig.add_axes([0.10, 0.08, 0.80, 0.37])
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.axis("off")

    rows = [
        (0.70, "T2 leg", LEFT_T2, RIGHT_T2),
        (0.32, "T3 leg", LEFT_T3, RIGHT_T3),
    ]
    for y, label, left_colour, right_colour in rows:
        ax2.text(0.02, y, label, ha="left", va="center", fontsize=8, color="#555555")
        ax2.add_patch(Circle((0.38, y), 0.072, facecolor=left_colour, edgecolor="white", linewidth=0.8))
        ax2.add_patch(Circle((0.70, y), 0.072, facecolor=right_colour, edgecolor="white", linewidth=0.8))
        ax2.plot([0.30, 0.22], [y, y + 0.08], color=left_colour, linewidth=1.1)
        ax2.plot([0.78, 0.86], [y, y + 0.08], color=right_colour, linewidth=1.1)

    ax2.text(0.38, 0.04, "left", ha="center", va="bottom", fontsize=8)
    ax2.text(0.70, 0.04, "right", ha="center", va="bottom", fontsize=8)
    ax2.text(0.27, 0.90, "left", ha="center", va="center", fontsize=7.2, color=LEFT_T2)
    ax2.text(0.81, 0.90, "right", ha="center", va="center", fontsize=7.2, color=RIGHT_T2)

    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE\t{output}")


if __name__ == "__main__":
    main()
