#!/usr/bin/env python3
"""Fit descriptive 1n/3n AF centroids and plot genome-wide normalized depth."""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap, LogNorm
import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar
from scipy.special import betaln, gammaln, logsumexp


VERSION = "1.2.0"
BLUE = "#0072B2"
ORANGE = "#D55E00"
PURPLE = "#7B3294"
GREY = "#666666"


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", required=True, type=Path,
                   help="Completed NDJ_LR_difference_scan_v1 run directory")
    p.add_argument("--contig", default="NC_140403.1")
    p.add_argument("--chromosome-label", default="Chromosome 15")
    p.add_argument("--start", type=int, default=200_000,
                   help="0-based, inclusive regional start")
    p.add_argument("--end", type=int, default=2_100_000,
                   help="0-based, exclusive regional end")
    p.add_argument("--error-rate", type=float, default=0.001)
    p.add_argument("--bins", type=int, default=20)
    p.add_argument("--validate-only", action="store_true")
    return p


def stamp(path):
    p = Path(path).resolve()
    s = p.stat()
    if s.st_size == 0:
        raise ValueError(f"Empty input: {p}")
    return {"path": str(p), "bytes": s.st_size, "mtime_ns": s.st_mtime_ns}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    partial = path.with_name(path.name + ".partial")
    partial.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    partial.replace(path)


def json_clean(value):
    """Convert NumPy/Pandas scalar values into strict JSON-compatible types."""
    if isinstance(value, dict):
        return {str(key): json_clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_clean(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_clean(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    return value


def load_inputs(args):
    run = args.run.resolve()
    config_path = run / "config.json"
    complete_path = run / "COMPLETE.json"
    windows_path = run / "windows.tsv"
    for path in [config_path, complete_path, windows_path]:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty required input: {path}")
    config = json.loads(config_path.read_text())
    complete = json.loads(complete_path.read_text())
    if complete.get("status") != "exploratory_scan_complete_not_biological_validation":
        raise RuntimeError("The supplied discovery run is not marked complete")
    matches = [
        (i, row) for i, row in enumerate(config["chromosomes"])
        if row["contig"] == args.contig
    ]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one config entry for {args.contig}; found {len(matches)}")
    index, chromosome = matches[0]
    counts_path = run / "chromosomes" / f"{index:02d}" / "counts.npz"
    if not counts_path.is_file() or counts_path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing chromosome counts: {counts_path}")
    windows = pd.read_csv(windows_path, sep="\t")
    required_windows = {
        "chrom", "start", "end", "left_depth_relative", "right_depth_relative",
        "usable_depth", "contig_background",
    }
    missing = required_windows - set(windows.columns)
    if missing:
        raise RuntimeError(f"windows.tsv lacks columns: {sorted(missing)}")
    with np.load(counts_path, allow_pickle=False) as handle:
        required_counts = {
            "pos", "left_ref", "left_alt", "right_ref", "right_alt",
        }
        missing = required_counts - set(handle.files)
        if missing:
            raise RuntimeError(f"counts.npz lacks arrays: {sorted(missing)}")
        counts = {key: handle[key] for key in required_counts}
    lengths = {key: len(value) for key, value in counts.items()}
    if len(set(lengths.values())) != 1:
        raise RuntimeError(f"Count array length mismatch: {lengths}")
    keep = (counts["pos"] > args.start) & (counts["pos"] <= args.end)
    counts = {key: value[keep] for key, value in counts.items()}
    if len(counts["pos"]) < 100:
        raise RuntimeError(f"Only {len(counts['pos'])} eligible regional sites; need at least 100")
    nl = counts["left_ref"] + counts["left_alt"]
    nr = counts["right_ref"] + counts["right_alt"]
    if np.any(nl <= 0) or np.any(nr <= 0):
        raise RuntimeError("Zero AD depth passed the discovery scan")
    counts["left_total"] = nl
    counts["right_total"] = nr
    counts["left_AF"] = counts["left_alt"] / nl
    counts["right_AF"] = counts["right_alt"] / nr
    regional_windows = windows.loc[
        windows["chrom"].eq(args.contig)
        & windows["start"].ge(args.start)
        & windows["end"].le(args.end)
        & windows["usable_depth"]
    ].copy()
    if regional_windows.empty:
        raise RuntimeError("No usable depth windows in selected region")
    inputs = [config_path, complete_path, windows_path, counts_path]
    palette = Path(__file__).with_name("224_batlowW_r.json")
    if not palette.is_file():
        raise FileNotFoundError(f"Missing bundled palette: {palette}")
    inputs.append(palette)
    return run, config, chromosome, windows, regional_windows, counts, palette, inputs


def beta_binomial_logpmf(k, n, q, concentration, error_rate):
    q = error_rate + (1 - 2 * error_rate) * q
    q = np.clip(q, 1e-8, 1 - 1e-8)
    alpha = q * concentration
    beta = (1 - q) * concentration
    return (
        gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
        + betaln(k + alpha, n - k + beta)
        - betaln(alpha, beta)
    )


def centroids(ploidy_left, ploidy_right, c_left, c_right):
    rows = []
    for left_dosage in range(ploidy_left + 1):
        p_left = left_dosage / ploidy_left
        for right_dosage in range(ploidy_right + 1):
            p_right = right_dosage / ploidy_right
            q_left = (1 - c_left) * p_left + c_left * p_right
            q_right = c_right * p_left + (1 - c_right) * p_right
            rows.append((left_dosage, right_dosage, q_left, q_right))
    return np.asarray(rows, dtype=float)


def model_log_likelihood(counts, ploidy_left, ploidy_right,
                         c_left, c_right, concentration, error_rate):
    centres = centroids(ploidy_left, ploidy_right, c_left, c_right)
    ql = centres[:, 2][None, :]
    qr = centres[:, 3][None, :]
    left = beta_binomial_logpmf(
        counts["left_alt"][:, None], counts["left_total"][:, None],
        ql, concentration, error_rate,
    )
    right = beta_binomial_logpmf(
        counts["right_alt"][:, None], counts["right_total"][:, None],
        qr, concentration, error_rate,
    )
    # Equal latent dosage-state weights make this a transparent centroid model,
    # not a pedigree/genotype-frequency model.
    return float(np.sum(logsumexp(left + right, axis=1) - math.log(len(centres))))


def fit_model(counts, ploidy_left, ploidy_right, error_rate):
    starts = [
        (0.01, 0.01, math.log(25)),
        (0.05, 0.05, math.log(25)),
        (0.10, 0.10, math.log(25)),
        (0.15, 0.15, math.log(25)),
        (0.20, 0.20, math.log(25)),
        (0.10, 0.20, math.log(25)),
        (0.20, 0.10, math.log(25)),
        (0.30, 0.15, math.log(50)),
        (0.15, 0.30, math.log(50)),
    ]

    def objective(theta):
        return -model_log_likelihood(
            counts, ploidy_left, ploidy_right,
            float(theta[0]), float(theta[1]), float(np.exp(theta[2])), error_rate,
        )

    results = []
    for start in starts:
        fit = minimize(
            objective, start, method="L-BFGS-B",
            bounds=[(0, 0.45), (0, 0.45), (math.log(2), math.log(2000))],
            options={"maxiter": 500, "ftol": 1e-10},
        )
        if fit.success and np.isfinite(fit.fun):
            results.append(fit)
    if not results:
        raise RuntimeError(f"No converged fit for {ploidy_left}n/{ploidy_right}n")
    best = min(results, key=lambda item: item.fun)
    c_left, c_right, log_concentration = best.x
    log_likelihood = -float(best.fun)
    n = len(counts["pos"])
    return {
        "model": f"{ploidy_left}n_left__{ploidy_right}n_right",
        "ploidy_left": ploidy_left,
        "ploidy_right": ploidy_right,
        "c_left": float(c_left),
        "c_right": float(c_right),
        "concentration": float(np.exp(log_concentration)),
        "log_likelihood": log_likelihood,
        "parameters": 3,
        "BIC": float(-2 * log_likelihood + 3 * math.log(n)),
    }


def fit_no_mixing(counts, error_rate):
    def objective(log_concentration):
        return -model_log_likelihood(
            counts, 1, 3, 0, 0, float(np.exp(log_concentration)), error_rate,
        )
    result = minimize_scalar(
        objective, bounds=(math.log(2), math.log(2000)), method="bounded",
        options={"xatol": 1e-8},
    )
    if not result.success or not np.isfinite(result.fun):
        raise RuntimeError("No converged no-mixing 1n/3n fit")
    return {
        "c_left": 0.0,
        "c_right": 0.0,
        "concentration": float(np.exp(result.x)),
        "log_likelihood": -float(result.fun),
    }


def nearest_distances(counts, centres):
    observed = np.column_stack([counts["left_AF"], counts["right_AF"]])
    expected = centres[:, 2:4]
    squared = np.sum((observed[:, None, :] - expected[None, :, :]) ** 2, axis=2)
    return np.sqrt(np.min(squared, axis=1))


def fit_all_models(counts, error_rate):
    rows = []
    for ploidy_left in [1, 2, 3]:
        for ploidy_right in [1, 2, 3]:
            row = fit_model(counts, ploidy_left, ploidy_right, error_rate)
            centres = centroids(
                ploidy_left, ploidy_right, row["c_left"], row["c_right"],
            )
            distance = nearest_distances(counts, centres)
            row["median_nearest_AF_distance"] = float(np.median(distance))
            row["AF_distance_90th_percentile"] = float(np.quantile(distance, 0.9))
            rows.append(row)
            print(
                f"FIT\t{row['model']}\tcL={row['c_left']:.4f}"
                f"\tcR={row['c_right']:.4f}\tBIC={row['BIC']:.2f}",
                flush=True,
            )
    table = pd.DataFrame(rows).sort_values("BIC").reset_index(drop=True)
    table["delta_BIC"] = table["BIC"] - table["BIC"].min()
    table["BIC_rank"] = np.arange(1, len(table) + 1)
    return table


def colourmap(palette):
    colors = np.asarray(json.loads(Path(palette).read_text()), dtype=float) / 255
    cmap = ListedColormap(colors[10:])
    cmap.set_bad("white")
    return cmap


def density_limits(left, right, bins):
    fig, ax = plt.subplots()
    collection = ax.hexbin(
        left, right, gridsize=bins, extent=(0, 1, 0, 1), mincnt=1,
    )
    maximum = max(2, int(np.max(collection.get_array())))
    plt.close(fig)
    return maximum


def draw_density(ax, counts, args, palette, vmax, title):
    collection = ax.hexbin(
        counts["left_AF"], counts["right_AF"],
        gridsize=args.bins, extent=(0, 1, 0, 1), mincnt=1,
        cmap=colourmap(palette), norm=LogNorm(vmin=1, vmax=vmax),
        linewidths=0, rasterized=True,
    )
    ax.set(
        xlim=(-0.025, 1.025), ylim=(-0.025, 1.025),
        xlabel="Left ALT-read fraction", ylabel="Right ALT-read fraction",
        title=title,
    )
    ax.set_aspect("equal")
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_yticks(np.linspace(0, 1, 6))
    return collection


def overlay_centroids(ax, centres, color, label_prefix, filled=False):
    for left_dosage, right_dosage, q_left, q_right in centres:
        ax.scatter(
            q_left, q_right, s=55,
            facecolor=color if filled else "white",
            edgecolor="white" if filled else color,
            linewidth=1.4, zorder=8,
        )
        ax.annotate(
            f"{int(left_dosage)}/1,{int(right_dosage)}/3",
            (q_left, q_right), xytext=(4, 4), textcoords="offset points",
            fontsize=6.5, color=color, zorder=9,
        )
    ax.plot([], [], "o", markerfacecolor=color if filled else "white",
            markeredgecolor=color, label=label_prefix)


def make_af_figure(counts, model_table, no_mixing, args, palette):
    target = model_table.loc[
        model_table["model"].eq("1n_left__3n_right")
    ].iloc[0]
    unfitted = centroids(1, 3, 0, 0)
    fitted = centroids(1, 3, target.c_left, target.c_right)
    vmax = density_limits(counts["left_AF"], counts["right_AF"], args.bins)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), constrained_layout=True)
    for ax, title in zip(
        axes,
        ["All eligible regional SNPs", "Unmixed 1n Left / 3n Right",
         "Fitted reciprocal mixing"],
    ):
        image = draw_density(ax, counts, args, palette, vmax, title)
    overlay_centroids(axes[1], unfitted, GREY, "Unmixed dosage centroids")
    overlay_centroids(axes[2], fitted, PURPLE, "Fitted dosage centroids", filled=True)
    axes[1].legend(frameon=False, fontsize=8, loc="lower right")
    axes[2].legend(frameon=False, fontsize=8, loc="lower right")
    colorbar = fig.colorbar(image, ax=axes, shrink=0.73, pad=0.012)
    colorbar.set_label("Sites per hexagon (log scale)")
    rank = int(target.BIC_rank)
    fig.suptitle(
        f"{args.chromosome_label}: {args.start/1e6:.1f}-{args.end/1e6:.1f} Mb | "
        "1n Left / 3n Right centroid test",
        fontsize=16, fontweight="bold",
    )
    fig.text(
        0.5, 0.005,
        f"n={len(counts['pos']):,} unselected eligible SNPs | "
        f"fitted cL={target.c_left:.3f}, cR={target.c_right:.3f} | "
        f"1n/3n BIC rank={rank}/9 | "
        f"delta logL vs unmixed 1n/3n={target.log_likelihood-no_mixing['log_likelihood']:.1f}. "
        "Equal latent dosage-state weights; descriptive centroid model, not a pedigree likelihood.",
        ha="center", fontsize=8.5, color="#444444",
    )
    return fig, target, unfitted, fitted


def make_depth_figure(windows, regional_windows, config, args):
    data = windows.loc[
        windows["usable_depth"]
        & windows["left_depth_relative"].gt(0)
        & windows["right_depth_relative"].gt(0)
    ].copy()
    background = data["contig_background"].astype(bool)
    # Work in linear relative-coverage space. First center each side on its
    # genome-wide autosomal mean, then subtract Left from Right window by
    # window. Shared coverage fluctuations consequently cancel. An ideal
    # 1n-Left/3n-Right region has (1.5 - 1.0) - (0.5 - 1.0) = +1.0; the
    # reverse orientation has delta -1.0.
    left_center = float(data.loc[background, "left_depth_relative"].mean())
    right_center = float(data.loc[background, "right_depth_relative"].mean())
    data["left_delta"] = data["left_depth_relative"] - left_center
    data["right_delta"] = data["right_depth_relative"] - right_center
    data["right_minus_left_delta"] = data["right_delta"] - data["left_delta"]
    offsets = {}
    ticks = []
    labels = []
    cursor = 0
    gap = 300_000
    for row in config["chromosomes"]:
        offsets[row["contig"]] = cursor
        ticks.append((cursor + row["length"] / 2) / 1e6)
        label = "Z" if str(row["label"]) == "21" else str(row["label"])
        labels.append(label)
        cursor += row["length"] + gap
    data["x"] = (
        data.apply(lambda row: offsets[row["chrom"]], axis=1)
        + (data["start"] + data["end"]) / 2
    ) / 1e6
    fig, ax = plt.subplots(figsize=(15, 5.1), constrained_layout=True)
    for index, row in enumerate(config["chromosomes"]):
        if index % 2:
            start = offsets[row["contig"]] / 1e6
            end = (offsets[row["contig"]] + row["length"]) / 1e6
            ax.axvspan(start, end, color="#F3F3F3", zorder=0)
    x0 = (offsets[args.contig] + args.start) / 1e6
    x1 = (offsets[args.contig] + args.end) / 1e6
    ax.axvspan(x0, x1, color="#E8D7F0", alpha=0.75, zorder=1,
               label=f"{args.chromosome_label} selected region")
    ax.scatter(
        data["x"], data["right_minus_left_delta"],
        s=8, color=PURPLE, alpha=0.60, linewidths=0,
        label="Right - Left", zorder=3,
    )
    baseline = 0.0
    expected_left_1n_right_3n = 1.0
    expected_left_3n_right_1n = -1.0
    ax.axhline(baseline, color="#222222", linewidth=0.9, zorder=2)
    ax.axhline(
        expected_left_1n_right_3n,
        color=ORANGE, linestyle="--", linewidth=1.2, zorder=2,
    )
    ax.axhline(
        expected_left_3n_right_1n,
        color=BLUE, linestyle="--", linewidth=1.2, zorder=2,
    )
    ax.text(
        cursor / 1e6, expected_left_1n_right_3n,
        " predicted 1n Left / 3n Right (+1.000)", color=ORANGE,
        fontsize=8, va="center", ha="right", backgroundcolor="white",
    )
    ax.text(
        cursor / 1e6, expected_left_3n_right_1n,
        " predicted 3n Left / 1n Right (-1.000)", color=BLUE,
        fontsize=8, va="center", ha="right", backgroundcolor="white",
    )
    ax.text(cursor / 1e6, 0, " equal Left / Right coverage", color="#222222",
            fontsize=8, va="bottom", ha="right", backgroundcolor="white")
    ax.set_xticks(ticks, labels)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Delta normalized depth (Right - Left; autosomal mean = 0)")
    ax.set_title(
        "Genome-wide Right-minus-Left delta coverage",
        fontsize=15, fontweight="bold",
    )
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.6, zorder=0)
    ax.legend(frameon=False, ncol=2, loc="upper right")
    values = data["right_minus_left_delta"].to_numpy()
    lower = min(-1.20, float(np.quantile(values, 0.005) - 0.08))
    upper = max(1.20, float(np.quantile(values, 0.995) + 0.08))
    ax.set_ylim(max(-2, lower), min(2, upper))
    selected = data.loc[
        data["chrom"].eq(args.contig)
        & data["start"].ge(args.start)
        & data["end"].le(args.end)
    ]
    regional_left = float(selected["left_delta"].median())
    regional_right = float(selected["right_delta"].median())
    regional_delta = float(selected["right_minus_left_delta"].median())
    text = (
        f"Selected-region Right - Left delta median: {regional_delta:+.3f}\n"
        f"Raw relative-depth medians: Left {regional_windows.left_depth_relative.median():.3f}; "
        f"Right {regional_windows.right_depth_relative.median():.3f}"
    )
    ax.text(0.012, 0.025, text, transform=ax.transAxes, fontsize=8.5,
            va="bottom", bbox={"facecolor": "white", "alpha": 0.86, "edgecolor": "none"})
    depth = {
        "depth_metric": "(Right relative depth - Right autosomal mean) - (Left relative depth - Left autosomal mean)",
        "left_autosomal_relative_depth_mean": left_center,
        "right_autosomal_relative_depth_mean": right_center,
        "regional_left_delta_depth_median": regional_left,
        "regional_right_delta_depth_median": regional_right,
        "regional_right_minus_left_delta_depth_median": regional_delta,
        "regional_left_relative_depth_median": float(regional_windows.left_depth_relative.median()),
        "regional_right_relative_depth_median": float(regional_windows.right_depth_relative.median()),
        "depth_implied_c_left": float(regional_windows.left_depth_relative.median() - 0.5),
        "depth_implied_c_right": float(1.5 - regional_windows.right_depth_relative.median()),
        "predicted_1n_left_3n_right_delta_depth": expected_left_1n_right_3n,
        "predicted_3n_left_1n_right_delta_depth": expected_left_3n_right_1n,
        "display_y_limits": list(map(float, ax.get_ylim())),
        "display_clipped_points": int(
            np.sum((values < ax.get_ylim()[0]) | (values > ax.get_ylim()[1]))
        ),
    }
    return fig, depth


def main():
    args = parser().parse_args()
    if args.start < 0 or args.end <= args.start:
        raise ValueError("Require 0 <= start < end")
    if not 0 < args.error_rate < 0.1 or args.bins < 10:
        raise ValueError("Invalid error rate or bin count")
    (run, config, chromosome, windows, regional_windows,
     counts, palette, inputs) = load_inputs(args)
    print(
        f"PREFLIGHT_PASS\tcontig={args.contig}\tregion={args.start}-{args.end}"
        f"\tsites={len(counts['pos'])}", flush=True,
    )
    if args.validate_only:
        return
    output = run / "figures" / "chr15_1n3n_followup"
    output.mkdir(parents=True, exist_ok=True)
    model_table = fit_all_models(counts, args.error_rate)
    no_mixing = fit_no_mixing(counts, args.error_rate)
    model_table.to_csv(output / "centroid_model_comparison.tsv", sep="\t", index=False)
    af_fig, target, unfitted, fitted = make_af_figure(
        counts, model_table, no_mixing, args, palette,
    )
    depth_fig, depth = make_depth_figure(
        windows, regional_windows, config, args,
    )
    af_pdf = output / "NDJ_chr15_1n3n_AF_centroid_fit.pdf"
    af_png = output / "NDJ_chr15_1n3n_AF_centroid_fit.png"
    depth_pdf = output / "NDJ_delta_depth_manhattan_1n3n.pdf"
    depth_png = output / "NDJ_delta_depth_manhattan_1n3n.png"
    combined_pdf = output / "NDJ_chr15_1n3n_AF_and_depth_test.pdf"
    af_fig.savefig(af_pdf)
    af_fig.savefig(af_png, dpi=220)
    depth_fig.savefig(depth_pdf)
    depth_fig.savefig(depth_png, dpi=220)
    with PdfPages(combined_pdf) as pdf:
        pdf.savefig(af_fig)
        pdf.savefig(depth_fig)
    plt.close(af_fig)
    plt.close(depth_fig)
    assignments = fitted.copy()
    pd.DataFrame(
        assignments,
        columns=["left_alt_dosage", "right_alt_dosage", "expected_left_AF", "expected_right_AF"],
    ).to_csv(output / "fitted_1n3n_centroids.tsv", sep="\t", index=False)
    audit = {
        "version": VERSION,
        "created": datetime.now(timezone.utc).isoformat(),
        "run": str(run),
        "region": {
            "chromosome": args.chromosome_label,
            "contig": args.contig,
            "start_0based_inclusive": args.start,
            "end_0based_exclusive": args.end,
            "sites": len(counts["pos"]),
        },
        "AF_definition": "FORMAT/AD_ALT / (AD_REF + AD_ALT)",
        "AF_site_selection": "Every eligible discovery-scan site in the region; no candidate or posterior filtering",
        "centroid_model": {
            "likelihood": "Equal-weight mixture over integer dosage pairs with beta-binomial read emissions",
            "warning": "Descriptive centroid comparison only; no pedigree, genotype-frequency, linkage, or founder prior",
            "error_rate": args.error_rate,
            "target_1n3n": target.to_dict(),
            "target_no_mixing": no_mixing,
            "delta_log_likelihood_vs_no_mixing": float(target.log_likelihood - no_mixing["log_likelihood"]),
            "unmixed_centroids": unfitted.tolist(),
            "fitted_centroids": fitted.tolist(),
        },
        "depth": depth,
        "inputs": [stamp(path) for path in inputs],
        "outputs": {},
    }
    for path in [af_pdf, af_png, depth_pdf, depth_png, combined_pdf,
                 output / "centroid_model_comparison.tsv",
                 output / "fitted_1n3n_centroids.tsv"]:
        audit["outputs"][path.name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    atomic_json(
        output / "NDJ_chr15_1n3n_followup.audit.json",
        json_clean(audit),
    )
    print("COMPLETE", flush=True)
    print(f"COMBINED_PDF\t{combined_pdf.resolve()}", flush=True)
    print(f"AF_PDF\t{af_pdf.resolve()}", flush=True)
    print(f"DEPTH_PDF\t{depth_pdf.resolve()}", flush=True)
    print(
        f"FIT_1N3N\tcL={target.c_left:.5f}\tcR={target.c_right:.5f}"
        f"\tBIC_rank={int(target.BIC_rank)}/9\tdelta_BIC={target.delta_BIC:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr, flush=True)
        raise
