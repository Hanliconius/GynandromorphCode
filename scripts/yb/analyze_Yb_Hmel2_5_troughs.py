#!/usr/bin/env python3
"""Interrogate the normalized Hmel2.5 Yb tracks at 50-bp resolution."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FEATURES = {
    "Yb_d": (1_403_500, 1_405_500),
    "Krista_2024_color_region": (1_408_879, 1_429_719),
    "cortex_proximal_TSS": (1_424_680, 1_424_680),
    "Yb_v": (1_429_500, 1_430_500),
}


def contiguous_segments(frame: pd.DataFrame, threshold: float) -> pd.DataFrame:
    rows = []
    for sample, data in frame.groupby("sample", sort=False):
        data = data.sort_values("start").copy()
        low = data["normalized_coverage"].le(threshold)
        group = (low.ne(low.shift(fill_value=False))).cumsum()

        for _, block in data.loc[low].groupby(group[low], sort=False):
            rows.append(
                {
                    "sample": sample,
                    "side": block["side"].iloc[0],
                    "threshold": threshold,
                    "start": int(block["start"].min()),
                    "end": int(block["end"].max()),
                    "width_bp": int(block["end"].max() - block["start"].min()),
                    "n_bins": int(len(block)),
                    "median_normalized_depth": float(
                        block["normalized_coverage"].median()
                    ),
                    "min_normalized_depth": float(
                        block["normalized_coverage"].min()
                    ),
                }
            )

    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--long", required=True, type=Path)
    ap.add_argument("--side-means", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    long = pd.read_csv(args.long, sep="\t")
    side = pd.read_csv(args.side_means, sep="\t")

    focus_start = 1_398_000
    focus_end = 1_435_000
    focus = long.loc[
        long["midpoint"].ge(focus_start)
        & long["midpoint"].le(focus_end)
    ].copy()

    segments = pd.concat(
        [
            contiguous_segments(focus, 0.50),
            contiguous_segments(focus, 0.25),
            contiguous_segments(focus, 0.10),
        ],
        ignore_index=True,
    )

    seg_path = args.output_dir / "Yb_Hmel2.5_low_coverage_segments.tsv"
    segments.to_csv(seg_path, sep="\t", index=False)

    feature_rows = []
    for name, (start, end) in FEATURES.items():
        if start == end:
            select = focus.loc[
                (focus["start"] <= start)
                & (focus["end"] > start)
            ]
        else:
            select = focus.loc[
                (focus["end"] > start)
                & (focus["start"] < end)
            ]

        for sample, data in select.groupby("sample", sort=False):
            feature_rows.append(
                {
                    "feature": name,
                    "start": start,
                    "end": end,
                    "sample": sample,
                    "side": data["side"].iloc[0],
                    "n_bins": int(len(data)),
                    "median_normalized_depth": float(
                        data["normalized_coverage"].median()
                    ),
                    "mean_normalized_depth": float(
                        data["normalized_coverage"].mean()
                    ),
                    "min_normalized_depth": float(
                        data["normalized_coverage"].min()
                    ),
                    "fraction_below_0.5": float(
                        data["normalized_coverage"].lt(0.5).mean()
                    ),
                }
            )

    feature_table = pd.DataFrame(feature_rows)
    feature_path = args.output_dir / "Yb_Hmel2.5_feature_depth_summary.tsv"
    feature_table.to_csv(feature_path, sep="\t", index=False)

    fig, axes = plt.subplots(
        2, 1,
        figsize=(11.5, 6.4),
        sharex=True,
        sharey=True,
    )
    fig.subplots_adjust(
        left=0.10,
        right=0.985,
        bottom=0.12,
        top=0.88,
        hspace=0.16,
    )

    for ax, side_name in zip(axes, ("Left", "Right")):
        d = focus.loc[focus["side"].eq(side_name)]
        for sample, sample_data in d.groupby("sample", sort=False):
            ax.plot(
                sample_data["midpoint"] / 1e6,
                sample_data["normalized_coverage"],
                linewidth=0.75,
                alpha=0.55,
                label=sample,
            )

        sm = side.loc[
            side["side"].eq(side_name)
            & side["midpoint"].ge(focus_start)
            & side["midpoint"].le(focus_end)
        ]
        ax.plot(
            sm["midpoint"] / 1e6,
            sm["side_mean_normalized_coverage"],
            linewidth=1.8,
            label=f"{side_name} mean",
        )

        ax.axhline(1.0, linestyle="--", linewidth=0.8)
        ax.axhline(0.5, linestyle=":", linewidth=0.8)

        for name, (start, end) in FEATURES.items():
            if start == end:
                ax.axvline(start / 1e6, linewidth=0.7, alpha=0.5)
            else:
                ax.axvspan(start / 1e6, end / 1e6, alpha=0.06)

        ax.set_ylabel("Normalized\ncoverage")
        ax.set_title(side_name, loc="left", fontweight="bold")
        ax.set_ylim(0, 1.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", linewidth=0.4, alpha=0.25)
        ax.legend(frameon=False, fontsize=7, loc="upper right")

    axes[-1].set_xlim(focus_start / 1e6, focus_end / 1e6)
    axes[-1].set_xlabel("Hmel215003o position (Mb)")

    fig.text(
        0.10,
        0.965,
        "Figure-8 Yb interval: normalized Hmel2.5 coverage",
        ha="left",
        va="top",
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.10,
        0.925,
        (
            "50-bp coverage from Krista historical bigWigs, divided by each "
            "library's mean coverage across Hmel215003o."
        ),
        ha="left",
        va="top",
        fontsize=8.5,
    )

    pdf = args.output_dir / "Yb_Hmel2.5_Figure8_targeted_zoom.pdf"
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    print("YB_HMEL2_5_TROUGH_ANALYSIS_COMPLETE")
    print(f"SEGMENTS\t{seg_path.resolve()}")
    print(f"FEATURES\t{feature_path.resolve()}")
    print(f"PDF\t{pdf.resolve()}")
    print()
    print("DEEPEST_SEGMENTS")
    if segments.empty:
        print("none")
    else:
        print(
            segments.sort_values(
                ["threshold", "min_normalized_depth", "width_bp"],
                ascending=[True, True, False],
            )
            .head(30)
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()
