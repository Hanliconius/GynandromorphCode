#!/usr/bin/env python3
"""Render the four raw-read HmRos Yb tracks in Krista's stacked log style."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SEQID = "NC_140403.1"
START = 1_620_000
END = 1_670_000
BIN_SIZE = 50
SAMPLES = (
    ("TN2L_S1", "TN2L_S1", "#0000FF"),
    ("TN2R_S2", "TN2R_S2", "#FF0000"),
    ("TN3L_S3", "TN3L_S3", "#228B22"),
    ("TN3R_S4", "TN3R_S4", "#FFA500"),
)
FEATURES = (
    ("Yb(d)", 1_633_041, 1_635_646),
    ("Krista interval start envelope", 1_638_549, 1_640_834),
    ("Yb(v)", 1_658_079, 1_659_072),
)


def read_bigwig(path: Path) -> tuple[np.ndarray, np.ndarray]:
    try:
        import pyBigWig
    except ImportError as exc:
        raise RuntimeError(
            "pyBigWig is required in the gynandrocalling Python environment"
        ) from exc

    handle = pyBigWig.open(str(path))
    try:
        chroms = handle.chroms()
        if SEQID not in chroms:
            raise RuntimeError(f"{SEQID} is absent from {path}; found {len(chroms)} contigs")
        if chroms[SEQID] < END:
            raise RuntimeError(f"{SEQID} is shorter than plot end {END}: {chroms[SEQID]}")
        values = np.asarray(handle.values(SEQID, START - 1, END, numpy=True), dtype=float)
    finally:
        handle.close()

    if values.size != END - START + 1:
        raise RuntimeError(f"Unexpected BigWig span length in {path}: {values.size}")
    values[~np.isfinite(values)] = 0.0
    # Align aggregation to the same 0-based genomic bins bamCoverage uses.
    # The plot window clips the first and last bins, so those edge means use
    # only the bases present inside the requested window.
    coordinates = np.arange(START - 1, END, dtype=np.int64)
    bin_ids, inverse = np.unique(coordinates // BIN_SIZE, return_inverse=True)
    totals = np.bincount(inverse, weights=values)
    counts = np.bincount(inverse)
    binned = totals / counts
    positions = bin_ids * BIN_SIZE + BIN_SIZE / 2.0 + 0.5
    return positions, binned


def render(profiles: dict[str, tuple[np.ndarray, np.ndarray]], output: Path) -> None:
    fig, axes = plt.subplots(
        nrows=4,
        ncols=1,
        figsize=(7.2, 7.0),
        sharex=True,
        gridspec_kw={"hspace": 0.08},
    )
    fig.subplots_adjust(left=0.15, right=0.985, top=0.92, bottom=0.12)

    for axis, (sample, label, color) in zip(axes, SAMPLES):
        x, raw = profiles[sample]
        y = np.full(raw.shape, np.nan, dtype=float)
        positive = raw > 0
        y[positive] = np.log(raw[positive])

        axis.axvspan(FEATURES[0][1] / 1e6, FEATURES[0][2] / 1e6,
                     color="#777777", alpha=0.14, linewidth=0, zorder=0)
        axis.axvspan(FEATURES[1][1] / 1e6, FEATURES[1][2] / 1e6,
                     color="#777777", alpha=0.05, linewidth=0, zorder=0)
        axis.axvspan(FEATURES[2][1] / 1e6, FEATURES[2][2] / 1e6,
                     color="#777777", alpha=0.14, linewidth=0, zorder=0)
        axis.axvline(1_658_299 / 1e6, color="#777777", linewidth=0.55,
                     linestyle=":", zorder=1)
        for boundary in FEATURES[1][1:]:
            axis.axvline(boundary / 1e6, color="#777777", linewidth=0.5,
                         linestyle=":", zorder=1)

        axis.plot(x / 1e6, y, color=color, linewidth=0.8, solid_capstyle="butt")
        axis.text(0.012, 0.84, label, transform=axis.transAxes,
                  ha="left", va="top", fontsize=8, color=color, fontweight="bold")
        axis.tick_params(axis="both", labelsize=7, width=0.55, length=2.5)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_linewidth(0.55)
        axis.spines["bottom"].set_linewidth(0.55)
        axis.grid(False)

    axes[-1].set_xlim(START / 1e6, END / 1e6)
    axes[-1].set_xlabel(f"Genomic position on {SEQID} (Mb)", fontsize=8)
    fig.supylabel("ln(coverage score)", x=0.035, fontsize=8)
    fig.suptitle("Yb locus coverage on H. melpomene rosina",
                 x=0.15, ha="left", fontsize=10)
    fig.text(0.985, 0.025,
             "Raw 50-bp bigWig scores; natural-log y transform. Shading marks Yb(d), "
             "the transferred Krista interval envelope, and Yb(v).",
             ha="right", va="bottom", fontsize=6.6)
    fig.savefig(output, format="pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    profiles: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    rows = []
    for sample, label, _ in SAMPLES:
        path = args.input_dir / f"{sample}.dedup.50bp.bw"
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing/empty BigWig: {path}")
        x, y = read_bigwig(path)
        profiles[sample] = (x, y)
        rows.append({
            "sample": sample,
            "display_label": label,
            "contig": SEQID,
            "window_start_1based": START,
            "window_end_1based": END,
            "bin_size_bp": BIN_SIZE,
            "bins": int(y.size),
            "positive_bins": int((y > 0).sum()),
            "window_median_bigwig_coverage": float(np.median(y)),
        })

    pdf = args.output_dir / "Yb_HmRos_Krista_style_raw_log_tracks.pdf"
    summary = args.output_dir / "Yb_HmRos_Krista_style_raw.summary.tsv"
    audit = args.output_dir / "Yb_HmRos_Krista_style_raw.audit.json"
    render(profiles, pdf)
    pd.DataFrame(rows).to_csv(summary, sep="\t", index=False)
    audit.write_text(json.dumps({
        "reference": "HmRos RefSeq GCF_056956725.1",
        "contig": SEQID,
        "window_1based_inclusive": [START, END],
        "coverage_source": "fresh BWA-MEM alignments of raw PRJNA1473931 PE reads; Picard MarkDuplicates",
        "bigwig": "deepTools bamCoverage, 50-bp bins; no normalization, MAPQ filter, or duplicate exclusion",
        "plot_transform": "natural logarithm of positive BigWig coverage scores",
        "track_order_and_colors": [
            {"sample": sample, "display_label": label, "color": color}
            for sample, label, color in SAMPLES
        ],
        "transferred_features_hmros": {
            name: [left, right] for name, left, right in FEATURES
        },
        "outputs": {"pdf": str(pdf), "summary": str(summary)},
    }, indent=2) + "\n")
    print("YB_RAW_KRISTA_STYLE_PLOT_COMPLETE")
    print(f"PDF\t{pdf}")
    print(f"SUMMARY\t{summary}")
    print(f"AUDIT\t{audit}")


if __name__ == "__main__":
    main()
