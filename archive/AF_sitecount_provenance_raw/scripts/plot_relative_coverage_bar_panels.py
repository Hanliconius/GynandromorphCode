#!/usr/bin/env python3

from pathlib import Path
import argparse

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

SAMPLE_ORDER = ["t2l", "t3l", "t2r", "t3r"]
SAMPLE_LABELS = {
    "t2l": "T2L",
    "t3l": "T3L",
    "t2r": "T2R",
    "t3r": "T3R",
}

SAMPLE_ALIASES = {
    "t2l": "t2l",
    "left_t2": "t2l",
    "tn_left_t2": "t2l",
    "tnleftt2": "t2l",

    "t3l": "t3l",
    "left_t3": "t3l",
    "tn_left_t3": "t3l",
    "tnleftt3": "t3l",

    "t2r": "t2r",
    "right_t2": "t2r",
    "tn_right_t2": "t2r",
    "tnrightt2": "t2r",

    "t3r": "t3r",
    "right_t3": "t3r",
    "tn_right_t3": "t3r",
    "tnrightt3": "t3r",
}

CLASS_ORDER = ["Autosomes", "W", "Z"]
CLASS_ALIASES = {
    "auto": "Autosomes",
    "autos": "Autosomes",
    "autosome": "Autosomes",
    "autosomes": "Autosomes",
    "a": "Autosomes",

    "w": "W",
    "chrw": "W",

    "z": "Z",
    "chrz": "Z",
}

SAMPLE_COLOURS = {
    "t2l": "#4C78A8",
    "t3l": "#72B7B2",
    "t2r": "#F58518",
    "t3r": "#E45756",
}


def canonicalise_sample(value):
    key = (
        str(value)
        .strip()
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
    )
    key = key.replace("__", "_")
    if key not in SAMPLE_ALIASES:
        raise ValueError(f"Unrecognised sample label: {value}")
    return SAMPLE_ALIASES[key]


def canonicalise_class(value):
    key = (
        str(value)
        .strip()
        .lower()
        .replace("-", "")
        .replace("_", "")
        .replace(" ", "")
    )
    if key not in CLASS_ALIASES:
        raise ValueError(f"Unrecognised chromosome class: {value}")
    return CLASS_ALIASES[key]


def first_present(columns, candidates, label):
    for candidate in candidates:
        if candidate in columns:
            return candidate
    raise ValueError(
        f"Could not find a {label} column. "
        f"Tried: {candidates}"
    )


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot observed versus expected relative coverage "
            "for Autosomes, W, and Z across T2L/T3L/T2R/T3R."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input TSV in long format.",
    )

    parser.add_argument(
        "--output-prefix",
        required=True,
        help="Output prefix for PDF/PNG.",
    )

    parser.add_argument(
        "--title",
        default="Observed versus expected relative coverage",
        help="Figure title.",
    )

    parser.add_argument(
        "--subtitle",
        default=(
            "Bars show observed relative coverage; open circles and black line show expected coverage; "
            "vertical segments show deviation."
        ),
        help="Figure subtitle.",
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    output_prefix = Path(args.output_prefix)

    df = pd.read_csv(input_path, sep="\t")

    sample_col = first_present(
        df.columns,
        [
            "sample",
            "sample_id",
            "sample_key",
            "track",
            "library",
        ],
        "sample",
    )

    class_col = first_present(
        df.columns,
        [
            "chrom_class",
            "class",
            "chromosome_class",
            "compartment",
            "target",
        ],
        "chromosome-class",
    )

    observed_col = first_present(
        df.columns,
        [
            "observed",
            "observed_coverage",
            "observed_relative_coverage",
            "relative_coverage",
        ],
        "observed",
    )

    expected_col = first_present(
        df.columns,
        [
            "expected",
            "expected_coverage",
            "expected_relative_coverage",
        ],
        "expected",
    )

    df = df.copy()
    df["sample_key"] = df[sample_col].map(canonicalise_sample)
    df["class_key"] = df[class_col].map(canonicalise_class)
    df["observed_value"] = df[observed_col].astype(float)
    df["expected_value"] = df[expected_col].astype(float)

    # Collapse safely if duplicate rows exist.
    grouped = (
        df.groupby(
            ["class_key", "sample_key"],
            as_index=False,
        )[["observed_value", "expected_value"]]
        .mean()
    )

    expected_combinations = {
        (chrom_class, sample_key)
        for chrom_class in CLASS_ORDER
        for sample_key in SAMPLE_ORDER
    }

    observed_combinations = set(
        zip(grouped["class_key"], grouped["sample_key"])
    )

    missing = sorted(expected_combinations - observed_combinations)

    if missing:
        raise RuntimeError(
            "Missing class/sample combinations: "
            + ", ".join(
                f"{chrom_class}:{sample_key}"
                for chrom_class, sample_key in missing
            )
        )

    ymax = max(
        grouped["observed_value"].max(),
        grouped["expected_value"].max(),
    )
    ymax = max(1.05, ymax * 1.15)

    fig = plt.figure(figsize=(12.5, 4.8))

    fig.text(
        0.06, 0.95,
        args.title,
        fontsize=18,
        fontweight="bold",
        ha="left",
    )

    fig.text(
        0.06, 0.90,
        args.subtitle,
        fontsize=10.5,
        color="#555555",
        ha="left",
    )

    axes = [
        fig.add_axes([0.08, 0.18, 0.24, 0.62]),
        fig.add_axes([0.39, 0.18, 0.24, 0.62]),
        fig.add_axes([0.70, 0.18, 0.24, 0.62]),
    ]

    for panel_index, chrom_class in enumerate(CLASS_ORDER):
        ax = axes[panel_index]

        sub = (
            grouped.loc[grouped["class_key"] == chrom_class]
            .set_index("sample_key")
            .loc[SAMPLE_ORDER]
            .reset_index()
        )

        x = np.arange(len(SAMPLE_ORDER))
        observed = sub["observed_value"].to_numpy()
        expected = sub["expected_value"].to_numpy()

        ax.bar(
            x,
            observed,
            width=0.72,
            color=[SAMPLE_COLOURS[s] for s in SAMPLE_ORDER],
            edgecolor="black",
            linewidth=0.4,
            zorder=2,
        )

        # Deviation segments and expected markers/line.
        for xi, obs, exp in zip(x, observed, expected):
            ax.plot(
                [xi, xi],
                [min(obs, exp), max(obs, exp)],
                color="#444444",
                linewidth=1.2,
                zorder=4,
            )

        ax.plot(
            x,
            expected,
            color="black",
            linewidth=1.0,
            zorder=5,
        )

        ax.scatter(
            x,
            expected,
            s=28,
            facecolors="white",
            edgecolors="black",
            linewidths=1.0,
            zorder=6,
        )

        ax.axhline(
            1.0,
            color="#999999",
            linestyle="--",
            linewidth=0.8,
            zorder=1,
        )

        ax.set_title(
            chrom_class,
            fontsize=13,
            pad=8,
            fontweight="bold",
        )

        ax.set_xticks(x)
        ax.set_xticklabels(
            [SAMPLE_LABELS[s] for s in SAMPLE_ORDER],
            fontsize=10.5,
        )

        ax.set_ylim(0, ymax)
        ax.grid(axis="y", color="#E5E5E5", linewidth=0.7)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        if panel_index == 0:
            ax.set_ylabel(
                "Relative coverage",
                fontsize=11.5,
            )
        else:
            ax.tick_params(labelleft=False)

    sample_handles = [
        Patch(
            facecolor=SAMPLE_COLOURS[s],
            edgecolor="black",
            linewidth=0.4,
            label=SAMPLE_LABELS[s],
        )
        for s in SAMPLE_ORDER
    ]

    style_handles = [
        Patch(
            facecolor="#BBBBBB",
            edgecolor="black",
            linewidth=0.4,
            label="Observed",
        ),
        Line2D(
            [0], [0],
            color="black",
            linewidth=1.0,
            marker="o",
            markerfacecolor="white",
            markeredgecolor="black",
            markersize=5,
            label="Expected",
        ),
        Line2D(
            [0], [0],
            color="#444444",
            linewidth=1.2,
            label="Deviation",
        ),
    ]

    fig.legend(
        handles=sample_handles,
        loc="upper left",
        bbox_to_anchor=(0.06, 0.02),
        frameon=False,
        ncol=4,
        fontsize=10,
        title="Samples",
        title_fontsize=10.5,
    )

    fig.legend(
        handles=style_handles,
        loc="upper right",
        bbox_to_anchor=(0.94, 0.02),
        frameon=False,
        ncol=3,
        fontsize=10,
        title="Display",
        title_fontsize=10.5,
    )

    pdf_file = output_prefix.with_suffix(".pdf")
    png_file = output_prefix.with_suffix(".png")

    fig.savefig(pdf_file)
    fig.savefig(png_file, dpi=220)

    print("INPUT", input_path.resolve(), sep="\t")
    print("OUTPUT_PDF", pdf_file.resolve(), sep="\t")
    print("OUTPUT_PNG", png_file.resolve(), sep="\t")
    print("SAMPLE_COLUMN", sample_col, sep="\t")
    print("CLASS_COLUMN", class_col, sep="\t")
    print("OBSERVED_COLUMN", observed_col, sep="\t")
    print("EXPECTED_COLUMN", expected_col, sep="\t")


if __name__ == "__main__":
    main()
