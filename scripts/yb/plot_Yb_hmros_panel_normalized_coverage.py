#!/usr/bin/env python3
"""Panel-of-normals normalized coverage across the HmRos Yb/cortex locus."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


FOCAL_ORDER = [
    "TN_Left_T2",
    "TN_Left_T3",
    "TN_Right_T2",
    "TN_Right_T3",
]

DISPLAY = {
    "TN_Left_T2": "Left T2",
    "TN_Left_T3": "Left T3",
    "TN_Right_T2": "Right T2",
    "TN_Right_T3": "Right T3",
}


def robust_mad(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.nan
    med = np.median(values)
    return float(np.median(np.abs(values - med)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bins", required=True, type=Path)
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--coordinates", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--control-rmad-max", type=float, default=0.35)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    bins = pd.read_csv(args.bins, sep="\t")
    samples = pd.read_csv(args.samples, sep="\t")
    coords = pd.read_csv(args.coordinates, sep="\t").iloc[0]

    required_bins = {
        "sample", "seqid", "start", "end", "midpoint",
        "mean_depth", "mapped_bases"
    }
    missing = required_bins - set(bins.columns)
    if missing:
        raise RuntimeError(f"Missing bin columns: {sorted(missing)}")

    required_samples = {"sample", "role", "side", "leg"}
    missing = required_samples - set(samples.columns)
    if missing:
        raise RuntimeError(f"Missing sample columns: {sorted(missing)}")

    bins = bins.merge(
        samples[["sample", "role", "side", "leg"]],
        on="sample",
        how="left",
        validate="many_to_one",
    )

    if bins["role"].isna().any():
        raise RuntimeError("Some coverage rows lack sample metadata.")

    mapped_bases = bins[["sample", "mapped_bases"]].drop_duplicates()
    target = float(mapped_bases["mapped_bases"].median())

    if not np.isfinite(target) or target <= 0:
        raise RuntimeError("Invalid median mapped-bases normalization target.")

    bins["library_scale"] = target / bins["mapped_bases"]
    bins["depth_library_normalized"] = (
        bins["mean_depth"] * bins["library_scale"]
    )

    controls = bins.loc[bins["role"].eq("control")].copy()
    focal = bins.loc[bins["role"].eq("focal")].copy()

    if controls["sample"].nunique() < 10:
        raise RuntimeError(
            f"Too few controls for panel normalization: "
            f"{controls['sample'].nunique()}"
        )

    keys = ["seqid", "start", "end", "midpoint"]

    control_stats = (
        controls.groupby(keys, as_index=False)
        ["depth_library_normalized"]
        .agg(
            control_median="median",
            control_n="count",
        )
    )

    mad_rows = []
    for key, group in controls.groupby(keys, sort=False):
        values = group["depth_library_normalized"].to_numpy(dtype=float)
        med = float(np.nanmedian(values))
        mad = robust_mad(values)
        rmad = mad / med if np.isfinite(med) and med > 0 else np.nan
        mad_rows.append(
            dict(
                zip(keys, key),
                control_mad=mad,
                control_rmad=rmad,
            )
        )

    control_mad = pd.DataFrame(mad_rows)
    control_stats = control_stats.merge(
        control_mad,
        on=keys,
        how="left",
        validate="one_to_one",
    )

    typical_control_depth = float(
        control_stats["control_median"].median()
    )

    control_stats["stable_control_bin"] = (
        control_stats["control_n"].ge(10)
        & control_stats["control_median"].ge(
            0.20 * typical_control_depth
        )
        & control_stats["control_rmad"].le(
            args.control_rmad_max
        )
    )

    focal = focal.merge(
        control_stats,
        on=keys,
        how="left",
        validate="many_to_one",
    )

    focal["panel_normalized_depth"] = (
        focal["depth_library_normalized"]
        / focal["control_median"]
    )

    focal.loc[
        ~focal["stable_control_bin"],
        "panel_normalized_depth",
    ] = np.nan

    focal = focal.sort_values(
        ["sample", "midpoint"],
        kind="stable",
    )

    focal["rolling_1kb"] = (
        focal.groupby("sample", sort=False)["panel_normalized_depth"]
        .transform(
            lambda s: s.rolling(
                window=4,
                center=True,
                min_periods=2,
            ).median()
        )
    )

    side_mean = (
        focal.groupby(
            ["side", *keys],
            as_index=False,
            sort=False,
        )["panel_normalized_depth"]
        .mean()
        .rename(
            columns={
                "panel_normalized_depth": "side_mean_normalized_depth"
            }
        )
    )

    out_long = (
        args.output_dir
        / "Yb_hmros_panel_normalized_250bp.long.tsv.gz"
    )
    out_side = (
        args.output_dir
        / "Yb_hmros_side_means_250bp.tsv.gz"
    )
    out_control = (
        args.output_dir
        / "Yb_hmros_control_panel_qc.tsv.gz"
    )
    out_pdf = (
        args.output_dir
        / "Yb_hmros_panel_normalized_coverage.pdf"
    )
    out_summary = (
        args.output_dir
        / "Yb_hmros_panel_normalized_coverage.summary.tsv"
    )
    out_audit = (
        args.output_dir
        / "Yb_hmros_panel_normalized_coverage.audit.json"
    )

    focal.to_csv(
        out_long,
        sep="\t",
        index=False,
        compression="gzip",
    )
    side_mean.to_csv(
        out_side,
        sep="\t",
        index=False,
        compression="gzip",
    )
    control_stats.to_csv(
        out_control,
        sep="\t",
        index=False,
        compression="gzip",
    )

    cortex_start = int(coords["cortex_start"])
    cortex_end = int(coords["cortex_end"])

    summary_rows = []
    for sample in FOCAL_ORDER:
        d = focal.loc[focal["sample"].eq(sample)].copy()
        c = d.loc[
            d["midpoint"].ge(cortex_start)
            & d["midpoint"].le(cortex_end)
        ]
        summary_rows.append(
            {
                "sample": sample,
                "display": DISPLAY[sample],
                "mapped_bases": int(d["mapped_bases"].iloc[0]),
                "stable_bins_total": int(
                    d["panel_normalized_depth"].notna().sum()
                ),
                "cortex_bins": int(
                    c["panel_normalized_depth"].notna().sum()
                ),
                "cortex_median_panel_normalized_depth": (
                    float(c["panel_normalized_depth"].median())
                    if c["panel_normalized_depth"].notna().any()
                    else np.nan
                ),
            }
        )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out_summary, sep="\t", index=False)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(11.5, 6.7),
        sharex=True,
        sharey=True,
    )
    fig.subplots_adjust(
        left=0.10,
        right=0.985,
        bottom=0.12,
        top=0.88,
        hspace=0.18,
    )

    x0 = int(coords["plot_start"]) / 1e6
    x1 = int(coords["plot_end"]) / 1e6
    c0 = cortex_start / 1e6
    c1 = cortex_end / 1e6

    for axis, side in zip(axes, ("Left", "Right")):
        sd = focal.loc[focal["side"].eq(side)]

        for sample in FOCAL_ORDER:
            if not sample.startswith("TN_" + side):
                continue
            d = sd.loc[sd["sample"].eq(sample)]
            axis.plot(
                d["midpoint"] / 1e6,
                d["panel_normalized_depth"],
                linewidth=0.55,
                alpha=0.35,
                label=DISPLAY[sample] + " bins",
            )
            axis.plot(
                d["midpoint"] / 1e6,
                d["rolling_1kb"],
                linewidth=1.0,
                alpha=0.85,
                label=DISPLAY[sample] + " 1-kb median",
            )

        sm = side_mean.loc[side_mean["side"].eq(side)]
        axis.plot(
            sm["midpoint"] / 1e6,
            sm["side_mean_normalized_depth"],
            linewidth=1.8,
            label=side + " mean",
        )

        axis.axhline(1.0, linewidth=0.8, linestyle="--")
        axis.axhline(0.5, linewidth=0.7, linestyle=":")
        axis.axvspan(c0, c1, alpha=0.08)

        unstable = control_stats.loc[
            ~control_stats["stable_control_bin"]
        ]
        for row in unstable.itertuples(index=False):
            axis.axvspan(
                row.start / 1e6,
                row.end / 1e6,
                alpha=0.035,
            )

        axis.set_xlim(x0, x1)
        axis.set_ylim(0, 1.8)
        axis.set_ylabel(
            "Panel-normalized\ncoverage"
        )
        axis.set_title(side, loc="left", fontweight="bold")
        axis.legend(
            loc="upper right",
            frameon=False,
            fontsize=7,
            ncol=2,
        )
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(axis="y", linewidth=0.4, alpha=0.25)

    axes[-1].set_xlabel(
        f"{coords['seqid']} position (Mb)"
    )

    fig.text(
        0.10,
        0.965,
        "Yb/cortex coverage on HmRos",
        ha="left",
        va="top",
        fontsize=16,
        fontweight="bold",
    )
    fig.text(
        0.10,
        0.925,
        (
            "GCF_056956725.1; 250-bp depth bins, MAPQ≥20 and baseQ≥20. "
            "Each library is scaled by mapped bases, then divided by the "
            "median of 18 independently normalized controls at each bin."
        ),
        ha="left",
        va="top",
        fontsize=8.6,
    )
    fig.text(
        0.985,
        0.02,
        (
            "Dashed = diploid expectation (1.0); dotted = half depth (0.5). "
            "Cortex gene shaded. Unstable control-panel bins are faintly shaded."
        ),
        ha="right",
        va="bottom",
        fontsize=7.4,
    )

    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    audit = {
        "assembly": "GCF_056956725.1",
        "assembly_name": "ASM5695672v1",
        "normalization": {
            "depth_filters": {
                "minimum_mapping_quality": 20,
                "minimum_base_quality": 20,
            },
            "bin_size_bp": 250,
            "library_size_factor":
                "median mapped bases across samples / sample mapped bases",
            "local_expected_depth":
                "median library-normalized depth across 18 controls per bin",
            "reported_ratio":
                "focal library-normalized depth / control median at same bin",
            "control_bin_filter":
                f"n>=10; median>=20% typical local control depth; relative MAD<={args.control_rmad_max}",
            "side_mean":
                "mean of T2 and T3 after each library is normalized independently",
        },
        "historical_coordinate_note":
            "Hmel2.5 Yb intervals are provenance only and are not directly reused as HmRos coordinates.",
        "coordinates_file": str(args.coordinates.resolve()),
        "outputs": {
            "pdf": str(out_pdf.resolve()),
            "summary": str(out_summary.resolve()),
            "long": str(out_long.resolve()),
            "side_means": str(out_side.resolve()),
            "control_qc": str(out_control.resolve()),
        },
    }

    out_audit.write_text(json.dumps(audit, indent=2) + "\n")

    print("YB_HMROS_NORMALIZED_COVERAGE_COMPLETE")
    print(summary.to_string(index=False))
    print(f"PDF\t{out_pdf.resolve()}")
    print(f"SUMMARY\t{out_summary.resolve()}")
    print(f"AUDIT\t{out_audit.resolve()}")


if __name__ == "__main__":
    main()
