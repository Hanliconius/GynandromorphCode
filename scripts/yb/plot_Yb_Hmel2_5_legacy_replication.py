#!/usr/bin/env python3
"""Render Hmel2.5 Yb tracks from the frozen 50-bp normalized table.

The track order, sample colors, stacked layout, and ascending genomic axis follow
Krista Pipho's 2024 Optix coverage R Markdown in the upstream repository. The
input signal is the historical Hmel2.5 bigWig data normalized by each sample's
mean depth across the full Hmel215003o scaffold. It is not re-aligned BAM-level
coverage and therefore does not claim to reproduce a MAPQ>=30 BAM filter.
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
START = 1_395_000
END = 1_440_000
BIN_SIZE = 50

# This is the order used by arrangeGrob(L2, R2, L3, R3, ncol=1) upstream.
SAMPLES = (
    ("TN2L_S1", "TN2L_S1", "#0000FF"),
    ("TN2R_S2", "TN2R_S2", "#FF0000"),
    ("TN3L_S3", "TN3L_S3", "#228B22"),
    ("TN3R_S4", "TN3R_S4", "#FFA500"),
)

FEATURES = (
    ("Yb(d)", 1_403_500, 1_405_500, 0.14),
    ("Krista color interval", 1_408_879, 1_429_719, 0.055),
    ("Yb(v)", 1_429_500, 1_430_500, 0.14),
)


def load_profiles(path: Path) -> dict[str, pd.DataFrame]:
    table = pd.read_csv(path, sep="\t", compression="infer")
    required = {
        "chrom", "start", "end", "midpoint", "normalized_coverage", "sample"
    }
    missing = required - set(table.columns)
    if missing:
        raise RuntimeError(f"Input table lacks columns: {sorted(missing)}")

    table = table.loc[
        table["chrom"].astype(str).eq(CHROM)
        & pd.to_numeric(table["midpoint"], errors="coerce").between(START, END)
    ].copy()
    for column in ("start", "end", "midpoint", "normalized_coverage"):
        table[column] = pd.to_numeric(table[column], errors="raise")

    profiles = {}
    for sample, _, _ in SAMPLES:
        profile = table.loc[table["sample"].astype(str).eq(sample)].copy()
        if profile.empty:
            raise RuntimeError(f"No rows for required sample {sample}")
        profile = profile.sort_values("midpoint", kind="stable")
        if profile["midpoint"].duplicated().any():
            raise RuntimeError(f"Duplicate bin midpoints for {sample}")
        steps = np.diff(profile["midpoint"].to_numpy(dtype=float))
        if steps.size and not np.allclose(steps, BIN_SIZE):
            raise RuntimeError(f"Expected {BIN_SIZE}-bp bins for {sample}")
        if not np.isfinite(profile["normalized_coverage"]).all():
            raise RuntimeError(f"Non-finite normalized coverage for {sample}")
        profiles[sample] = profile

    return profiles


def add_features(axis: plt.Axes) -> None:
    for label, left, right, alpha in FEATURES:
        axis.axvspan(
            left / 1e6,
            right / 1e6,
            color="#777777",
            alpha=alpha,
            linewidth=0,
            zorder=0,
        )


def render(
    profiles: dict[str, pd.DataFrame],
    output: Path,
    scale: str,
) -> None:
    fig, axes = plt.subplots(
        nrows=4,
        ncols=1,
        figsize=(7.2, 5.4),
        sharex=True,
        sharey=(scale == "linear"),
        gridspec_kw={"hspace": 0.08},
    )
    fig.subplots_adjust(left=0.16, right=0.985, top=0.91, bottom=0.12)

    for axis, (sample, label, color) in zip(axes, SAMPLES):
        data = profiles[sample]
        x = data["midpoint"].to_numpy(dtype=float) / 1e6
        y = data["normalized_coverage"].to_numpy(dtype=float)

        add_features(axis)
        if scale == "log":
            y = np.where(y > 0, np.log(y), np.nan)
        axis.plot(x, y, color=color, linewidth=0.75, solid_capstyle="butt", zorder=2)
        axis.axhline(
            0.0 if scale == "log" else 1.0,
            color="#333333",
            linewidth=0.55,
            linestyle="--",
            zorder=1,
        )

        axis.text(
            0.012, 0.82, label,
            transform=axis.transAxes,
            ha="left", va="top",
            fontsize=8,
            color=color,
            fontweight="bold",
        )
        axis.set_xlim(START / 1e6, END / 1e6)
        axis.tick_params(axis="both", labelsize=7, width=0.55, length=2.5)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_linewidth(0.55)
        axis.spines["bottom"].set_linewidth(0.55)
        axis.grid(False)

    axes[-1].set_xlabel("Genomic position on Hmel215003o (Mb)", fontsize=8)
    if scale == "log":
        fig.supylabel("ln(normalized coverage)", x=0.035, fontsize=8)
        fig.suptitle("Yb locus coverage on Hmel2.5", x=0.16, ha="left", fontsize=10)
    else:
        fig.supylabel("Coverage / full-scaffold mean", x=0.035, fontsize=8)
        fig.suptitle("Yb locus coverage on Hmel2.5", x=0.16, ha="left", fontsize=10)

    fig.savefig(output, format="pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    if not args.input.is_file() or args.input.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty normalized input: {args.input}")

    profiles = load_profiles(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    linear_pdf = args.output_dir / "Yb_Hmel2.5_Figure8_legacy_tracks_linear.pdf"
    log_pdf = args.output_dir / "Yb_Hmel2.5_Figure8_legacy_tracks_log.pdf"
    render(profiles, linear_pdf, "linear")
    render(profiles, log_pdf, "log")

    summary_rows = []
    for sample, label, _ in SAMPLES:
        profile = profiles[sample]
        row = {
            "sample": sample,
            "display_label": label,
            "bins": int(len(profile)),
            "window_median_normalized_coverage": float(
                profile["normalized_coverage"].median()
            ),
        }
        for feature, left, right, _ in FEATURES:
            inside = profile.loc[
                profile["midpoint"].between(left, right, inclusive="left"),
                "normalized_coverage",
            ]
            row[f"{feature.replace(' ', '_')}_median_normalized_coverage"] = (
                float(inside.median()) if not inside.empty else np.nan
            )
        summary_rows.append(row)

    summary_path = args.output_dir / "Yb_Hmel2.5_Figure8_legacy_tracks.summary.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)

    audit = {
        "reference": "Hmel2.5",
        "scaffold": CHROM,
        "region_1based_inclusive_for_display": [START, END],
        "bin_size_bp": BIN_SIZE,
        "normalization": "50-bp coverage divided by each sample's full-scaffold mean",
        "source_signal": "historical Krista Pipho 2024 bigWig tracks",
        "bam_mapq_filter_applied": False,
        "layout_source": (
            "Krista Pipho 2024 Optix coverage plot R Markdown: "
            "four stacked tracks in TN2L, TN2R, TN3L, TN3R order; "
            "ascending genomic coordinate; blue, red, forest green, orange"
        ),
        "log_view": "natural log of normalized coverage; zero-depth bins omitted",
        "feature_intervals_hmel2_5": [
            {"name": name, "start": left, "end": right}
            for name, left, right, _ in FEATURES
        ],
        "outputs": {
            "linear_pdf": str(linear_pdf.resolve()),
            "log_pdf": str(log_pdf.resolve()),
            "summary_tsv": str(summary_path.resolve()),
        },
    }
    audit_path = args.output_dir / "Yb_Hmel2.5_Figure8_legacy_tracks.audit.json"
    audit_path.write_text(json.dumps(audit, indent=2) + "\n")

    print("YB_HMEL2_5_LEGACY_TRACK_RENDER_COMPLETE")
    print(f"LINEAR_PDF\t{linear_pdf.resolve()}")
    print(f"LOG_PDF\t{log_pdf.resolve()}")
    print(f"SUMMARY\t{summary_path.resolve()}")
    print(f"AUDIT\t{audit_path.resolve()}")


if __name__ == "__main__":
    main()
