#!/usr/bin/env python3
"""Normalize Krista's historical Hmel2.5 coverage tracks like Figure 8.

The available historical inputs are bigWigs, not the original BAMs used by
Livraghi et al. (2021). Therefore this reproduces the key Figure-8
normalization exactly (divide each individual's local coverage by that
individual's mean depth across the full cortex scaffold), while using the
coverage signal already stored in Krista's bigWigs.

The paper's original BAM-level analysis used mosdepth, 50-bp windows and
MAPQ>=30. Those BAM-level filters cannot be retroactively imposed on an
already-created bigWig, so the audit records that distinction explicitly.
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
PLOT_START = 1_395_000
PLOT_END = 1_440_000
BIN_SIZE = 50

FEATURES = {
    "Yb_d": (1_403_500, 1_405_500),
    "Krista_2024_color_region": (1_408_879, 1_429_719),
    "Yb_v": (1_429_500, 1_430_500),
}

SAMPLES = (
    ("TN2L_S1", "Left T2", "Left"),
    ("TN3L_S3", "Left T3", "Left"),
    ("TN2R_S2", "Right T2", "Right"),
    ("TN3R_S4", "Right T3", "Right"),
)


def load_bg(path: Path) -> pd.DataFrame:
    table = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "coverage"],
    )
    table = table.loc[table["chrom"].astype(str).eq(CHROM)].copy()
    for col in ("start", "end", "coverage"):
        table[col] = pd.to_numeric(table[col], errors="raise")
    if table.empty:
        raise RuntimeError(f"No {CHROM} records in {path}")
    return table


def make_bins(table: pd.DataFrame) -> pd.DataFrame:
    starts = np.arange(PLOT_START, PLOT_END, BIN_SIZE, dtype=int)
    ends = np.minimum(starts + BIN_SIZE, PLOT_END)
    numer = np.zeros(len(starts), dtype=float)

    target = table.loc[
        table["end"].gt(PLOT_START)
        & table["start"].lt(PLOT_END)
    ]

    for row in target.itertuples(index=False):
        s = max(int(row.start), PLOT_START)
        e = min(int(row.end), PLOT_END)
        if e <= s:
            continue

        first = (s - PLOT_START) // BIN_SIZE
        last = (e - 1 - PLOT_START) // BIN_SIZE

        for idx in range(first, last + 1):
            bs = int(starts[idx])
            be = int(ends[idx])
            overlap = max(0, min(e, be) - max(s, bs))
            if overlap:
                numer[idx] += overlap * float(row.coverage)

    widths = (ends - starts).astype(float)

    return pd.DataFrame(
        {
            "chrom": CHROM,
            "start": starts,
            "end": ends,
            "midpoint": (starts + ends) / 2.0,
            # Missing bigWig spans are explicitly treated as zero coverage.
            "raw_coverage": numer / widths,
        }
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bedgraph-dir", required=True, type=Path)
    parser.add_argument("--scaffold-stats", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    stats = pd.read_csv(args.scaffold_stats, sep="\t")
    required = {"sample", "scaffold_length", "scaffold_mean_depth"}
    missing = required - set(stats.columns)
    if missing:
        raise RuntimeError(f"Missing scaffold-stat columns: {sorted(missing)}")

    profiles = {}
    long_rows = []
    summary_rows = []

    for sample, label, side in SAMPLES:
        bg = args.bedgraph_dir / f"{sample}.Hmel215003o.bedGraph"
        if not bg.is_file() or bg.stat().st_size == 0:
            raise RuntimeError(f"Missing/empty bedGraph: {bg}")

        row = stats.loc[stats["sample"].eq(sample)]
        if len(row) != 1:
            raise RuntimeError(f"Expected one scaffold-stat row for {sample}")

        scaffold_mean = float(row.iloc[0]["scaffold_mean_depth"])
        scaffold_length = int(row.iloc[0]["scaffold_length"])

        if not np.isfinite(scaffold_mean) or scaffold_mean <= 0:
            raise RuntimeError(
                f"Invalid Hmel215003o mean depth for {sample}: {scaffold_mean}"
            )

        profile = make_bins(load_bg(bg))
        profile["normalized_coverage"] = (
            profile["raw_coverage"] / scaffold_mean
        )
        profile["sample"] = sample
        profile["label"] = label
        profile["side"] = side
        profile["scaffold_mean_depth"] = scaffold_mean
        profile["scaffold_length"] = scaffold_length

        profiles[sample] = profile
        long_rows.append(profile)

        k0, k1 = FEATURES["Krista_2024_color_region"]
        yb = profile.loc[
            profile["midpoint"].ge(k0)
            & profile["midpoint"].le(k1)
        ]

        summary_rows.append(
            {
                "sample": sample,
                "label": label,
                "side": side,
                "scaffold_length": scaffold_length,
                "scaffold_mean_depth": scaffold_mean,
                "plot_region_median_normalized_depth":
                    float(profile["normalized_coverage"].median()),
                "Krista_color_region_median_normalized_depth":
                    float(yb["normalized_coverage"].median()),
                "Krista_color_region_min_normalized_depth":
                    float(yb["normalized_coverage"].min()),
            }
        )

    long = pd.concat(long_rows, ignore_index=True)
    summary = pd.DataFrame(summary_rows)

    side = (
        long.groupby(
            ["side", "chrom", "start", "end", "midpoint"],
            as_index=False,
            sort=False,
        )["normalized_coverage"]
        .mean()
        .rename(
            columns={"normalized_coverage": "side_mean_normalized_coverage"}
        )
    )

    long_path = args.output_dir / "Yb_Hmel2.5_Figure8_normalized_50bp.long.tsv.gz"
    side_path = args.output_dir / "Yb_Hmel2.5_Figure8_normalized_50bp.side_means.tsv.gz"
    summary_path = args.output_dir / "Yb_Hmel2.5_Figure8_normalized.summary.tsv"
    pdf_path = args.output_dir / "Yb_Hmel2.5_Figure8_style_normalized.pdf"
    audit_path = args.output_dir / "Yb_Hmel2.5_Figure8_style_normalized.audit.json"

    long.to_csv(long_path, sep="\t", index=False, compression="gzip")
    side.to_csv(side_path, sep="\t", index=False, compression="gzip")
    summary.to_csv(summary_path, sep="\t", index=False)

    fig, axes = plt.subplots(
        2, 1,
        figsize=(11.5, 6.5),
        sharex=True,
        sharey=True,
    )
    fig.subplots_adjust(
        left=0.10,
        right=0.985,
        bottom=0.12,
        top=0.87,
        hspace=0.18,
    )

    for ax, side_name in zip(axes, ("Left", "Right")):
        for sample, label, side_name_sample in SAMPLES:
            if side_name_sample != side_name:
                continue
            d = profiles[sample]
            ax.plot(
                d["midpoint"] / 1e6,
                d["normalized_coverage"],
                linewidth=0.8,
                alpha=0.55,
                label=label,
            )

        sm = side.loc[side["side"].eq(side_name)]
        ax.plot(
            sm["midpoint"] / 1e6,
            sm["side_mean_normalized_coverage"],
            linewidth=1.7,
            label=f"{side_name} mean",
        )

        ax.axhline(1.0, linewidth=0.8, linestyle="--")
        ax.axhline(0.5, linewidth=0.7, linestyle=":")

        k0, k1 = FEATURES["Krista_2024_color_region"]
        ax.axvspan(k0 / 1e6, k1 / 1e6, alpha=0.055)

        for key in ("Yb_d", "Yb_v"):
            f0, f1 = FEATURES[key]
            ax.axvspan(f0 / 1e6, f1 / 1e6, alpha=0.12)

        ax.set_ylabel("Normalized\ncoverage")
        ax.set_title(side_name, loc="left", fontweight="bold")
        ax.set_ylim(0, 2.0)
        ax.grid(axis="y", linewidth=0.4, alpha=0.25)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.legend(loc="upper right", frameon=False, fontsize=8)

    axes[-1].set_xlim(PLOT_START / 1e6, PLOT_END / 1e6)
    axes[-1].set_xlabel("Hmel215003o position (Mb)")

    fig.text(
        0.10, 0.965,
        "Yb coverage on Hmel2.5",
        ha="left", va="top",
        fontsize=16, fontweight="bold",
    )
    fig.text(
        0.10, 0.925,
        (
            "Krista historical bigWigs; 50-bp bins. Each library is divided "
            "by its own mean coverage across the full Hmel215003o scaffold, "
            "matching the normalization denominator used in Figure 8."
        ),
        ha="left", va="top",
        fontsize=8.7,
    )
    fig.text(
        0.985, 0.02,
        (
            "Dashed = scaffold mean (1.0); dotted = half depth (0.5). "
            "Broad shade = Krista 2024 colour-region interval; "
            "narrow shades = historical Yb(d)/Yb(v) intervals."
        ),
        ha="right", va="bottom",
        fontsize=7.4,
    )

    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    audit = {
        "reference": "Hmel2.5",
        "reference_subspecies": "Heliconius melpomene melpomene",
        "chromosome_scaffold": CHROM,
        "plot_region": [PLOT_START, PLOT_END],
        "bin_size_bp": BIN_SIZE,
        "normalization": (
            "each sample 50-bp local bigWig coverage divided by its own "
            "length-weighted mean coverage across the full Hmel215003o scaffold"
        ),
        "Figure8_method_comparison": {
            "paper":
                "BWA-mapped BAMs; mosdepth 50-bp windows; MAPQ>=30; "
                "window median depth divided by individual mean full-scaffold depth",
            "this_historical_replication":
                "Krista's already-created bigWigs; 50-bp rebinning; "
                "same full-scaffold mean-depth normalization. "
                "The original BAM-level MAPQ>=30 filter and per-base medians "
                "cannot be imposed retrospectively from the bigWigs.",
        },
        "upstream_bigwig_commit":
            "99376deebc25f4e0bb7aa9927c7969c693591ef4",
        "features": {
            key: {"start": val[0], "end": val[1]}
            for key, val in FEATURES.items()
        },
        "outputs": {
            "pdf": str(pdf_path.resolve()),
            "summary": str(summary_path.resolve()),
            "long": str(long_path.resolve()),
            "side_means": str(side_path.resolve()),
        },
    }
    audit_path.write_text(json.dumps(audit, indent=2) + "\n")

    print("YB_HMEL2_5_FIGURE8_NORMALIZATION_COMPLETE")
    print(summary.to_string(index=False))
    print(f"PDF\t{pdf_path.resolve()}")
    print(f"SUMMARY\t{summary_path.resolve()}")
    print(f"AUDIT\t{audit_path.resolve()}")


if __name__ == "__main__":
    main()
