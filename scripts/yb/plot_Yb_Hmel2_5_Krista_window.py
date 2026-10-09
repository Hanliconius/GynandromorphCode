#!/usr/bin/env python3
"""Plot the full historical Hmel2.5 Yb window directly from Krista bigWigs.

The existing bedGraphs are extracted from Krista's historical bigWigs. This
renderer rebins their stored values to 50 bp and applies no additional
per-library normalization. It produces linear-depth and log-y views on the
broader Hmel215003o:1.380-1.600 Mb window for comparison with the historical
plots.
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
START = 1_380_000
END = 1_600_000
BIN_SIZE = 50
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


def read_bedgraph(path: Path) -> pd.DataFrame:
    data = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "coverage"],
    )
    data = data.loc[data["chrom"].astype(str).eq(CHROM)].copy()
    for column in ("start", "end", "coverage"):
        data[column] = pd.to_numeric(data[column], errors="raise")
    if data.empty:
        raise RuntimeError(f"No {CHROM} records in {path}")
    return data


def rebin(data: pd.DataFrame) -> pd.DataFrame:
    starts = np.arange(START, END, BIN_SIZE, dtype=int)
    ends = np.minimum(starts + BIN_SIZE, END)
    sums = np.zeros(starts.size, dtype=float)
    target = data.loc[data["end"].gt(START) & data["start"].lt(END)]

    for row in target.itertuples(index=False):
        left = max(int(row.start), START)
        right = min(int(row.end), END)
        if right <= left:
            continue
        first = (left - START) // BIN_SIZE
        last = (right - 1 - START) // BIN_SIZE
        for index in range(first, last + 1):
            overlap = max(
                0,
                min(right, int(ends[index])) - max(left, int(starts[index])),
            )
            sums[index] += overlap * float(row.coverage)

    widths = (ends - starts).astype(float)
    return pd.DataFrame(
        {
            "chrom": CHROM,
            "start": starts,
            "end": ends,
            "midpoint": (starts + ends) / 2,
            "coverage": sums / widths,
        }
    )


def render(profiles: dict[str, pd.DataFrame], output: Path, log_y: bool) -> None:
    fig, axes = plt.subplots(
        nrows=4,
        ncols=1,
        figsize=(7.2, 7.0),
        sharex=True,
        gridspec_kw={"hspace": 0.08},
    )
    fig.subplots_adjust(left=0.16, right=0.985, top=0.93, bottom=0.12)

    for axis, (sample, label, color) in zip(axes, SAMPLES):
        profile = profiles[sample]
        y = profile["coverage"].to_numpy(dtype=float)
        if log_y:
            y = np.where(y > 0, y, np.nan)

        for _, left, right, alpha in FEATURES:
            axis.axvspan(
                left / 1e6,
                right / 1e6,
                color="#777777",
                alpha=alpha,
                linewidth=0,
                zorder=0,
            )
        axis.plot(
            profile["midpoint"].to_numpy(dtype=float) / 1e6,
            y,
            color=color,
            linewidth=0.65,
            solid_capstyle="butt",
            zorder=2,
        )
        if log_y:
            axis.set_yscale("log")
        axis.text(
            0.012,
            0.82,
            label,
            transform=axis.transAxes,
            ha="left",
            va="top",
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

    axes[-1].set_xlabel(f"Genomic position on {CHROM} (Mb)", fontsize=8)
    fig.supylabel("Coverage (bigWig values)", x=0.035, fontsize=8)
    fig.suptitle(
        "Yb coverage on Hmel2.5 - historical window",
        x=0.16,
        ha="left",
        fontsize=10,
    )
    fig.savefig(output, format="pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bedgraph-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    profiles: dict[str, pd.DataFrame] = {}
    summaries = []
    for sample, label, _ in SAMPLES:
        path = args.bedgraph_dir / f"{sample}.Hmel215003o.bedGraph"
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing/empty historical bedGraph: {path}")
        profile = rebin(read_bedgraph(path))
        profiles[sample] = profile
        summaries.append(
            {
                "sample": sample,
                "display_label": label,
                "window_start": START,
                "window_end": END,
                "bin_size_bp": BIN_SIZE,
                "bins": len(profile),
                "window_median_bigwig_coverage": float(profile["coverage"].median()),
            }
        )

    linear_pdf = args.output_dir / "Yb_Hmel2.5_Krista_window_linear.pdf"
    log_pdf = args.output_dir / "Yb_Hmel2.5_Krista_window_log_y.pdf"
    render(profiles, linear_pdf, log_y=False)
    render(profiles, log_pdf, log_y=True)

    summary = args.output_dir / "Yb_Hmel2.5_Krista_window.summary.tsv"
    pd.DataFrame(summaries).to_csv(summary, sep="\t", index=False)
    audit = {
        "reference": "Hmel2.5",
        "scaffold": CHROM,
        "window_1based_inclusive_for_display": [START, END],
        "bin_size_bp": BIN_SIZE,
        "source": "Krista historical bigWigs, via full-scaffold bedGraphs",
        "additional_normalization": "none; plot stored bigWig values directly",
        "views": {
            "linear": "linear y-axis of stored coverage values",
            "log_y": "logarithmic y-axis of positive stored coverage values",
        },
        "features": [
            {"name": name, "start": left, "end": right}
            for name, left, right, _ in FEATURES
        ],
        "outputs": {
            "linear_pdf": str(linear_pdf.resolve()),
            "log_y_pdf": str(log_pdf.resolve()),
            "summary_tsv": str(summary.resolve()),
        },
    }
    audit_path = args.output_dir / "Yb_Hmel2.5_Krista_window.audit.json"
    audit_path.write_text(json.dumps(audit, indent=2) + "\n")

    print("YB_HMEL2_5_KRISTA_WINDOW_COMPLETE")
    print(f"LINEAR_PDF\t{linear_pdf.resolve()}")
    print(f"LOG_Y_PDF\t{log_pdf.resolve()}")
    print(f"SUMMARY\t{summary.resolve()}")
    print(f"AUDIT\t{audit_path.resolve()}")


if __name__ == "__main__":
    main()
