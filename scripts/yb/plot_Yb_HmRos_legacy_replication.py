#!/usr/bin/env python3
"""Render a legacy-style HmRos view of the Yb coverage tracks.

This is a visual companion to the historical Hmel2.5 figure, not a strict
normalization-level replication. The HmRos input is the existing 250-bp,
MAPQ>=20/baseQ>=20 panel-normalized coverage table. Each focal track is shown
relative to the median of the 18-control panel at the same bin. The plot uses
native HmRos coordinates and marks the uncertain start boundary of Krista's
historical Hmel2.5 color interval as a bracket.
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


SEQID = "NC_140403.1"
START = 1_620_000
END = 1_670_000
BIN_SIZE = 250

# Preserve Krista's historical panel order and palette.
SAMPLES = (
    ("TN_Left_T2", "TN2L_S1", "#0000FF"),
    ("TN_Right_T2", "TN2R_S2", "#FF0000"),
    ("TN_Left_T3", "TN3L_S3", "#228B22"),
    ("TN_Right_T3", "TN3R_S4", "#FFA500"),
)

# Exact HmRos marker intervals from the Hmel2.5-to-HmRos transfer audit.
YB_D = (1_633_041, 1_635_646)
YB_V = (1_658_079, 1_659_072)
KRISTA_START_BRACKET = (1_638_549, 1_640_834)
KRISTA_END_EXACT = 1_658_299
CORTEX_TSS = 1_653_288


def load_profiles(path: Path) -> dict[str, pd.DataFrame]:
    table = pd.read_csv(path, sep="\t", compression="infer")
    required = {"sample", "seqid", "midpoint", "panel_normalized_depth"}
    missing = required - set(table.columns)
    if missing:
        raise RuntimeError(f"Input table lacks columns: {sorted(missing)}")

    table["midpoint"] = pd.to_numeric(table["midpoint"], errors="raise")
    table["panel_normalized_depth"] = pd.to_numeric(
        table["panel_normalized_depth"], errors="coerce"
    )
    table = table.loc[
        table["seqid"].astype(str).eq(SEQID)
        & table["midpoint"].between(START, END)
    ].copy()

    profiles: dict[str, pd.DataFrame] = {}
    for sample, _, _ in SAMPLES:
        profile = table.loc[table["sample"].astype(str).eq(sample)].copy()
        if profile.empty:
            raise RuntimeError(f"No rows for required HmRos sample {sample}")
        profile = profile.sort_values("midpoint", kind="stable")
        if profile["midpoint"].duplicated().any():
            raise RuntimeError(f"Duplicate bin midpoints for {sample}")
        steps = np.diff(profile["midpoint"].to_numpy(dtype=float))
        if steps.size and not np.allclose(steps, BIN_SIZE):
            raise RuntimeError(f"Expected {BIN_SIZE}-bp bins for {sample}")
        if profile["panel_normalized_depth"].notna().sum() < 20:
            raise RuntimeError(f"Too few stable control-panel bins for {sample}")
        profiles[sample] = profile

    return profiles


def add_features(axis: plt.Axes) -> None:
    # The historical color interval start is not exact in HmRos. Shade the
    # transferred envelope lightly and draw both start-bracket limits.
    axis.axvspan(
        KRISTA_START_BRACKET[0] / 1e6,
        KRISTA_END_EXACT / 1e6,
        color="#777777",
        alpha=0.045,
        linewidth=0,
        zorder=0,
    )
    for boundary in KRISTA_START_BRACKET:
        axis.axvline(
            boundary / 1e6,
            color="#777777",
            linewidth=0.5,
            linestyle=":",
            zorder=1,
        )
    axis.axvline(
        KRISTA_END_EXACT / 1e6,
        color="#777777",
        linewidth=0.5,
        linestyle=":",
        zorder=1,
    )
    for left, right in (YB_D, YB_V):
        axis.axvspan(
            left / 1e6,
            right / 1e6,
            color="#777777",
            alpha=0.14,
            linewidth=0,
            zorder=0,
        )
    axis.axvline(
        CORTEX_TSS / 1e6,
        color="#555555",
        linewidth=0.55,
        linestyle="--",
        zorder=1,
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
        y = data["panel_normalized_depth"].to_numpy(dtype=float)

        add_features(axis)
        if scale == "log":
            y = np.where(y > 0, np.log(y), np.nan)
        axis.plot(
            x,
            y,
            color=color,
            linewidth=0.75,
            solid_capstyle="butt",
            zorder=2,
        )
        axis.axhline(
            0.0 if scale == "log" else 1.0,
            color="#333333",
            linewidth=0.55,
            linestyle="--",
            zorder=1,
        )
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

    axes[-1].set_xlabel(f"Genomic position on {SEQID} (Mb)", fontsize=8)
    if scale == "log":
        fig.supylabel("ln(panel-normalized coverage)", x=0.035, fontsize=8)
        fig.suptitle("Yb locus coverage on HmRos", x=0.16, ha="left", fontsize=10)
    else:
        fig.supylabel("Coverage / HmRos control-panel median", x=0.035, fontsize=8)
        fig.suptitle("Yb locus coverage on HmRos", x=0.16, ha="left", fontsize=10)

    fig.text(
        0.985,
        0.025,
        "Shaded: Yb(d), historical color interval envelope, and Yb(v); "
        "dotted lines bracket the uncertain interval start; dashed line: cortex TSS.",
        ha="right",
        va="bottom",
        fontsize=6.6,
    )
    fig.savefig(output, format="pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    if not args.input.is_file() or args.input.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty HmRos normalized input: {args.input}")

    profiles = load_profiles(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    linear_pdf = args.output_dir / "Yb_HmRos_Figure8_legacy_tracks_linear.pdf"
    log_pdf = args.output_dir / "Yb_HmRos_Figure8_legacy_tracks_log.pdf"
    render(profiles, linear_pdf, "linear")
    render(profiles, log_pdf, "log")

    summary_rows = []
    for sample, label, _ in SAMPLES:
        values = profiles[sample]["panel_normalized_depth"].dropna()
        summary_rows.append(
            {
                "sample": sample,
                "display_label": label,
                "bins_total": int(len(profiles[sample])),
                "stable_bins": int(values.size),
                "window_median_panel_normalized_coverage": (
                    float(values.median()) if not values.empty else np.nan
                ),
            }
        )

    summary_path = args.output_dir / "Yb_HmRos_Figure8_legacy_tracks.summary.tsv"
    pd.DataFrame(summary_rows).to_csv(summary_path, sep="\t", index=False)

    audit = {
        "reference": "HmRos GCF_056956725.1 / ASM5695672v1",
        "seqid": SEQID,
        "region_1based_inclusive_for_display": [START, END],
        "bin_size_bp": BIN_SIZE,
        "source_signal": (
            "Existing HmRos focal tracks from the MAPQ>=20/baseQ>=20 250-bp "
            "panel-normalized coverage analysis; no realignment performed."
        ),
        "normalization": (
            "Each focal library is scaled by mapped bases, then divided by "
            "the median library-normalized depth of 18 controls at each bin. "
            "This differs from the historical Hmel2.5 full-scaffold-mean normalization."
        ),
        "track_order_and_colors": [
            {"sample": sample, "legacy_label": label, "color": color}
            for sample, label, color in SAMPLES
        ],
        "transferred_features_hmros": {
            "Yb_d_exact": list(YB_D),
            "historical_Krista_color_region_start_bracket": list(KRISTA_START_BRACKET),
            "historical_Krista_color_region_end_exact": KRISTA_END_EXACT,
            "cortex_proximal_TSS_exact": CORTEX_TSS,
            "Yb_v_exact": list(YB_V),
        },
        "log_view": "natural log of panel-normalized coverage; zero/negative bins omitted",
        "outputs": {
            "linear_pdf": str(linear_pdf.resolve()),
            "log_pdf": str(log_pdf.resolve()),
            "summary_tsv": str(summary_path.resolve()),
        },
    }
    audit_path = args.output_dir / "Yb_HmRos_Figure8_legacy_tracks.audit.json"
    audit_path.write_text(json.dumps(audit, indent=2) + "\n")

    print("YB_HMROS_LEGACY_TRACK_RENDER_COMPLETE")
    print(f"LINEAR_PDF\t{linear_pdf.resolve()}")
    print(f"LOG_PDF\t{log_pdf.resolve()}")
    print(f"SUMMARY\t{summary_path.resolve()}")
    print(f"AUDIT\t{audit_path.resolve()}")


if __name__ == "__main__":
    main()
