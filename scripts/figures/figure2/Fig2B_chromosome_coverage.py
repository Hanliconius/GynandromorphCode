#!/usr/bin/env python3
"""Figure 2B — mirrored chromosome-by-chromosome relative coverage.

Uses the established chromosome-centred layout, with the sample colours,
bar outlines, guide lines and typography restored from the earlier coverage
figures.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from common import add_panel_letter, panel_output, project_root

DEFAULT_INPUT = Path(
    "analysis/chromosome_coverage_bars/"
    "run_20260923T162015902039Z/chromosome_coverage.tsv"
)

LEFT_T2 = "#4C78A8"
LEFT_T3 = "#72B7B2"
RIGHT_T2 = "#F58518"
RIGHT_T3 = "#E45756"


def chromosome_rank(label: str):
    text = str(label).strip().upper()
    if text.isdigit():
        return (0, int(text), text)
    if text == "W":
        return (1, 0, text)
    if text == "Z":
        return (1, 1, text)
    return (2, 0, text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    root = project_root(args.project_root)
    input_path = Path(args.input)
    if not input_path.is_absolute():
        input_path = root / input_path
    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    output = panel_output(root, "B", args.output)

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })

    df = pd.read_csv(input_path, sep="\t")
    required = {
        "label",
        "left_t2_relative_depth",
        "left_t3_relative_depth",
        "right_t2_relative_depth",
        "right_t3_relative_depth",
    }
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"Coverage table missing columns: {sorted(missing)}")

    df = df.copy()
    df["_rank"] = df["label"].map(chromosome_rank)
    df = df.sort_values("_rank").reset_index(drop=True)

    labels = df["label"].astype(str).tolist()
    y = np.arange(len(df), dtype=float)
    columns = [
        "left_t2_relative_depth",
        "left_t3_relative_depth",
        "right_t2_relative_depth",
        "right_t3_relative_depth",
    ]
    values = np.concatenate([df[c].to_numpy(float) for c in columns])
    finite = values[np.isfinite(values)]
    if not len(finite):
        raise RuntimeError("No finite relative-depth values")

    xmax = max(1.10, float(np.nanmax(finite)) * 1.08)

    fig = plt.figure(figsize=(4.25, 6.15))
    add_panel_letter(fig, "B", x=0.018, y=0.988, size=15)
    ax = fig.add_axes([0.09, 0.12, 0.86, 0.81])

    bar_h = 0.15
    offsets = (-0.10, +0.10)
    bar_kw = dict(height=bar_h, edgecolor="black", linewidth=0.35, zorder=3)

    ax.barh(
        y + offsets[0],
        -df["left_t2_relative_depth"].to_numpy(float),
        color=LEFT_T2,
        label="T2L",
        **bar_kw,
    )
    ax.barh(
        y + offsets[1],
        -df["left_t3_relative_depth"].to_numpy(float),
        color=LEFT_T3,
        label="T3L",
        **bar_kw,
    )
    ax.barh(
        y + offsets[0],
        df["right_t2_relative_depth"].to_numpy(float),
        color=RIGHT_T2,
        label="T2R",
        **bar_kw,
    )
    ax.barh(
        y + offsets[1],
        df["right_t3_relative_depth"].to_numpy(float),
        color=RIGHT_T3,
        label="T3R",
        **bar_kw,
    )

    ax.axvline(-1.0, color="#999999", linestyle="--", linewidth=0.8, zorder=1)
    ax.axvline(+1.0, color="#999999", linestyle="--", linewidth=0.8, zorder=1)

    centre_half_width = 0.105 * xmax
    ax.axvspan(-centre_half_width, centre_half_width, color="white", zorder=4)
    ax.axvline(0.0, color="#444444", linewidth=0.8, zorder=5)

    for yi, label in zip(y, labels):
        special = label.upper() in {"W", "Z"}
        ax.text(
            0.0,
            yi,
            label,
            ha="center",
            va="center",
            fontsize=7.4,
            color="#B2182B" if special else "#333333",
            fontweight="bold" if special else "normal",
            zorder=6,
        )

    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-0.6, len(df) - 0.4)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{abs(value):.1f}"))
    ax.set_xlabel("Relative coverage", fontsize=10, labelpad=5)
    ax.tick_params(axis="x", labelsize=8, length=3)
    ax.grid(axis="x", color="#E5E5E5", linewidth=0.65, zorder=0)
    ax.set_axisbelow(True)

    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    fig.text(0.24, 0.075, "Left", ha="center", va="center", fontsize=9)
    fig.text(0.79, 0.075, "Right", ha="center", va="center", fontsize=9)
    fig.text(0.50, 0.952, "Chromosome", ha="center", va="bottom", fontsize=9.5)

    handles = [
        Patch(facecolor=LEFT_T2, edgecolor="black", linewidth=0.35, label="T2L"),
        Patch(facecolor=LEFT_T3, edgecolor="black", linewidth=0.35, label="T3L"),
        Patch(facecolor=RIGHT_T2, edgecolor="black", linewidth=0.35, label="T2R"),
        Patch(facecolor=RIGHT_T3, edgecolor="black", linewidth=0.35, label="T3R"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.52, 0.002),
        ncol=4,
        frameon=False,
        fontsize=7.4,
        handlelength=1.3,
        columnspacing=1.1,
    )

    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(f"INPUT\t{input_path}")
    print(f"WROTE\t{output}")


if __name__ == "__main__":
    main()
