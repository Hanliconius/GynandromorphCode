#!/usr/bin/env python3
"""Normalized coverage replication across the H. melpomene Yb/cortex region.

Inputs are bedGraph slices converted from the four historical Krista-Pipho
gynandromorph_2024 bigWig tracks.  Coverage is normalized within each library
to the median coverage of two flanking windows outside the annotated Yb/ivory
features.  This removes library-depth differences while preserving local copy-
number signal.

Coordinates are Hmel2.5 scaffold Hmel215003o, matching the historical bigWigs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


CHROM = "Hmel215003o"
REGION_START = 1_380_000
REGION_END = 1_600_000
BIN_SIZE = 250

# Flanks deliberately avoid the known mapped Yb intervals and the ivory
# deletion.  They provide a local diploid-depth baseline for each library.
NORMALIZATION_FLANKS = (
    (1_380_000, 1_395_000),
    (1_585_000, 1_600_000),
)

FEATURES = {
    "Yb_d": (1_403_500, 1_405_500),
    "Krista_2024_color_region": (1_408_879, 1_429_719),
    "Yb_v": (1_429_500, 1_430_500),
    "ivory_delta78k": (1_494_902, 1_573_177),
}

SAMPLES = (
    ("TN2L_S1", "T2L", "Left"),
    ("TN3L_S3", "T3L", "Left"),
    ("TN2R_S2", "T2R", "Right"),
    ("TN3R_S4", "T3R", "Right"),
)


def load_bedgraph(path: Path) -> pd.DataFrame:
    table = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "coverage"],
    )

    table = table.loc[
        table["chrom"].astype(str).eq(CHROM)
    ].copy()

    for column in ("start", "end", "coverage"):
        table[column] = pd.to_numeric(
            table[column],
            errors="raise",
        )

    table = table.loc[
        table["end"].gt(REGION_START)
        & table["start"].lt(REGION_END)
    ].copy()

    if table.empty:
        raise RuntimeError(
            f"No {CHROM}:{REGION_START}-{REGION_END} coverage in {path}"
        )

    return table


def interval_mean(
    table: pd.DataFrame,
    start: int,
    end: int,
) -> float:
    left = np.maximum(
        table["start"].to_numpy(dtype=float),
        float(start),
    )
    right = np.minimum(
        table["end"].to_numpy(dtype=float),
        float(end),
    )
    overlap = np.maximum(
        0.0,
        right - left,
    )
    values = table["coverage"].to_numpy(dtype=float)

    denom = overlap.sum()

    if denom <= 0:
        return np.nan

    return float(
        np.sum(overlap * values) / denom
    )


def binned_profile(
    table: pd.DataFrame,
) -> pd.DataFrame:
    starts = np.arange(
        REGION_START,
        REGION_END,
        BIN_SIZE,
        dtype=int,
    )
    ends = np.minimum(
        starts + BIN_SIZE,
        REGION_END,
    )

    coverage = np.array(
        [
            interval_mean(
                table,
                int(start),
                int(end),
            )
            for start, end in zip(
                starts,
                ends,
            )
        ],
        dtype=float,
    )

    return pd.DataFrame(
        {
            "start": starts,
            "end": ends,
            "midpoint": (
                starts + ends
            ) / 2.0,
            "raw_coverage": coverage,
        }
    )


def flank_baseline(
    profile: pd.DataFrame,
) -> float:
    midpoint = profile["midpoint"].to_numpy()
    mask = np.zeros(
        len(profile),
        dtype=bool,
    )

    for start, end in NORMALIZATION_FLANKS:
        mask |= (
            (midpoint >= start)
            & (midpoint < end)
        )

    values = profile.loc[
        mask,
        "raw_coverage",
    ].to_numpy(dtype=float)

    values = values[
        np.isfinite(values)
        & (values > 0)
    ]

    if values.size == 0:
        raise RuntimeError(
            "No finite positive flank coverage available for normalization"
        )

    baseline = float(
        np.median(values)
    )

    if not np.isfinite(baseline) or baseline <= 0:
        raise RuntimeError(
            f"Invalid normalization baseline: {baseline}"
        )

    return baseline


def feature_mean(
    profile: pd.DataFrame,
    start: int,
    end: int,
) -> float:
    midpoint = profile["midpoint"].to_numpy()
    mask = (
        (midpoint >= start)
        & (midpoint < end)
    )

    values = profile.loc[
        mask,
        "normalized_coverage",
    ].to_numpy(dtype=float)

    values = values[
        np.isfinite(values)
    ]

    return (
        float(np.mean(values))
        if values.size
        else np.nan
    )


def draw(
    profiles: dict[str, pd.DataFrame],
    output_pdf: Path,
):
    fig, axes = plt.subplots(
        nrows=4,
        ncols=1,
        figsize=(11.5, 8.2),
        sharex=True,
        sharey=True,
    )

    fig.subplots_adjust(
        left=0.11,
        right=0.985,
        bottom=0.10,
        top=0.91,
        hspace=0.12,
    )

    ymax = max(
        1.5,
        max(
            float(
                np.nanpercentile(
                    profile["normalized_coverage"],
                    99.5,
                )
            )
            for profile in profiles.values()
        )
        * 1.05,
    )

    for axis, (
        sample,
        short_label,
        side,
    ) in zip(
        axes,
        SAMPLES,
    ):
        profile = profiles[sample]

        x = (
            profile["midpoint"].to_numpy(dtype=float)
            / 1_000_000
        )
        y = profile[
            "normalized_coverage"
        ].to_numpy(dtype=float)

        axis.plot(
            x,
            y,
            linewidth=1.0,
        )

        axis.axhline(
            1.0,
            linewidth=0.8,
            linestyle="--",
        )
        axis.axhline(
            0.5,
            linewidth=0.6,
            linestyle=":",
        )

        # The ivory deletion is a copy-number diagnostic and therefore gets
        # a light span.  Yb(d)/Yb(v) are marked as narrow vertical intervals.
        d0, d1 = FEATURES[
            "ivory_delta78k"
        ]
        axis.axvspan(
            d0 / 1_000_000,
            d1 / 1_000_000,
            alpha=0.08,
        )

        for feature in (
            "Yb_d",
            "Yb_v",
        ):
            f0, f1 = FEATURES[
                feature
            ]
            axis.axvspan(
                f0 / 1_000_000,
                f1 / 1_000_000,
                alpha=0.10,
            )

        axis.set_ylim(
            0.0,
            ymax,
        )
        axis.set_ylabel(
            short_label,
            rotation=0,
            ha="right",
            va="center",
            labelpad=18,
            fontweight="bold",
        )
        axis.grid(
            axis="y",
            linewidth=0.4,
            alpha=0.25,
        )
        axis.spines["top"].set_visible(
            False
        )
        axis.spines["right"].set_visible(
            False
        )

    axes[-1].set_xlabel(
        "Hmel215003o position (Mb)"
    )

    fig.text(
        0.11,
        0.965,
        "Yb/cortex normalized coverage",
        ha="left",
        va="top",
        fontsize=16,
        fontweight="bold",
    )

    fig.text(
        0.11,
        0.935,
        (
            "Historical Krista coverage tracks; each library normalized to "
            "its own local flanking median (baseline = 1). "
            "Dashed line = 1× local depth; dotted line = 0.5×."
        ),
        ha="left",
        va="top",
        fontsize=9,
    )

    fig.text(
        0.985,
        0.015,
        (
            "Shaded broad interval: ivory Δ78k deletion "
            "(Hmel215003o:1,494,902–1,573,177). "
            "Narrow spans: mapped Yb(d) and Yb(v) intervals."
        ),
        ha="right",
        va="bottom",
        fontsize=7.5,
    )

    fig.savefig(
        output_pdf,
        bbox_inches="tight",
    )
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bedgraph-dir",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
    )
    args = parser.parse_args()

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    profiles = {}
    summary_rows = []
    long_rows = []

    for sample, short_label, side in SAMPLES:
        path = (
            args.bedgraph_dir
            / f"{sample}.dedup.Yb_region.bedGraph"
        )

        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing/empty bedGraph: {path}"
            )

        table = load_bedgraph(
            path
        )
        profile = binned_profile(
            table
        )
        baseline = flank_baseline(
            profile
        )

        profile[
            "normalized_coverage"
        ] = (
            profile[
                "raw_coverage"
            ]
            / baseline
        )

        profile["sample"] = sample
        profile["label"] = short_label
        profile["side"] = side
        profile[
            "normalization_baseline"
        ] = baseline

        profiles[sample] = profile
        long_rows.append(
            profile
        )

        record = {
            "sample": sample,
            "label": short_label,
            "side": side,
            "flank_median_raw_coverage":
                baseline,
        }

        for feature, (
            start,
            end,
        ) in FEATURES.items():
            record[
                f"{feature}_mean_normalized_coverage"
            ] = feature_mean(
                profile,
                start,
                end,
            )

        summary_rows.append(
            record
        )

    long_table = pd.concat(
        long_rows,
        ignore_index=True,
    )
    summary = pd.DataFrame(
        summary_rows
    )

    long_path = (
        args.output_dir
        / "Yb_normalized_coverage_250bp.long.tsv.gz"
    )
    summary_path = (
        args.output_dir
        / "Yb_normalized_coverage.summary.tsv"
    )
    pdf_path = (
        args.output_dir
        / "Yb_normalized_coverage.pdf"
    )
    audit_path = (
        args.output_dir
        / "Yb_normalized_coverage.audit.json"
    )

    long_table.to_csv(
        long_path,
        sep="\t",
        index=False,
        compression="gzip",
    )
    summary.to_csv(
        summary_path,
        sep="\t",
        index=False,
    )

    draw(
        profiles,
        pdf_path,
    )

    audit = {
        "coordinate_system":
            "Hmel2.5 / Hmel215003o",
        "region": {
            "chrom": CHROM,
            "start": REGION_START,
            "end": REGION_END,
            "bin_size": BIN_SIZE,
        },
        "normalization": {
            "method":
                "per-library median raw coverage in two local flanking windows",
            "flanks": [
                {
                    "start": start,
                    "end": end,
                }
                for start, end
                in NORMALIZATION_FLANKS
            ],
            "baseline_after_normalization":
                1.0,
        },
        "features": {
            key: {
                "start": value[0],
                "end": value[1],
            }
            for key, value
            in FEATURES.items()
        },
        "outputs": {
            "pdf": str(
                pdf_path.resolve()
            ),
            "summary": str(
                summary_path.resolve()
            ),
            "long": str(
                long_path.resolve()
            ),
        },
    }

    audit_path.write_text(
        json.dumps(
            audit,
            indent=2,
        )
        + "\n"
    )

    print(
        "YB_NORMALIZED_COVERAGE_COMPLETE"
    )
    print(
        summary.to_string(
            index=False
        )
    )
    print(
        f"PDF\t{pdf_path.resolve()}"
    )
    print(
        f"SUMMARY\t{summary_path.resolve()}"
    )
    print(
        f"AUDIT\t{audit_path.resolve()}"
    )


if __name__ == "__main__":
    main()
