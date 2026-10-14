#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm, Normalize
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT / "analysis/autosomal_three_vs_four_founder"
FIGURE_DIR = ANALYSIS / "narrative_figures"
SCRIPT159 = PROJECT / "scripts/159_fit_autosomal_three_vs_four_founder.py"
SCRIPT160 = PROJECT / "scripts/160_simulate_calibrate_autosomal_founder_models.py"
SCRIPT161 = PROJECT / "scripts/161_plot_founder_mixing_narrative.py"

MATCHED_FILE = ANALYSIS / "matched_autosomal_focal_counts.tsv.gz"
OBSERVED_FILE = ANALYSIS / "observed_nested_founder_fit.tsv"
SIMULATION_FILE = ANALYSIS / "founder_simulation_replicates.tsv.gz"
MIXING_FILE = ANALYSIS / "mixing_calibration.tsv"
BOOTSTRAP_FILE = ANALYSIS / "observed_10kb_block_bootstrap.tsv.gz"
SUMMARY_FILE = ANALYSIS / "calibrated_founder_and_mixing_summary.tsv"

METRICS_FILE = ANALYSIS / "founder_simulation_information_metrics.tsv.gz"
REPRESENTATIVE_FILE = ANALYSIS / "representative_founder_visualization_simulations.tsv"
STEP2_PDF = FIGURE_DIR / "step2_observed_and_simulated_AF_hex_comparison.pdf"
STEP2_PNG = FIGURE_DIR / "step2_observed_and_simulated_AF_hex_comparison.png"
STEP3_PDF = FIGURE_DIR / "step3_recovery_information_diagnostics.pdf"
STEP3_PNG = FIGURE_DIR / "step3_recovery_information_diagnostics.png"
STEP3_CAUSE_PNG = FIGURE_DIR / "step3_four_founder_recovery_failure_diagnostic.png"
NARRATIVE_V2 = FIGURE_DIR / "three_vs_four_founder_mixing_narrative_v2.pdf"
AUDIT_FILE = FIGURE_DIR / "founder_mixing_narrative_v2.audit.tsv"

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
BLACK = "#222222"
GREY = "#777777"
CMAP = "viridis"


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def summary_dictionary(path):
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if set(table.columns) != {"statistic", "value"}:
        raise RuntimeError(f"Unexpected summary columns in {path}")
    return dict(zip(table["statistic"], table["value"]))


def boolean_series(values):
    if pd.api.types.is_bool_dtype(values):
        return values.astype(bool)
    normalized = values.astype(str).str.strip().str.lower()
    result = normalized.map({"true": True, "false": False, "1": True, "0": False})
    if result.isna().any():
        bad = sorted(normalized.loc[result.isna()].unique())
        raise RuntimeError(f"Could not interpret Boolean values: {bad}")
    return result.astype(bool)


def regenerate_simulation_information(model159, model160, matched, simulations):
    matrix, donor, included, allele_frequency = model160.build_control_matrix(matched)
    emission = model159.load_emission()
    depths = {
        "left_total": matched["left_AD_total"].to_numpy(float),
        "right_total": matched["right_AD_total"].to_numpy(float),
    }
    counts = simulations.groupby("true_founder_model")["replicate"].nunique()
    if set(counts.index) != {"three_founder", "four_founder"}:
        raise RuntimeError("Simulation table does not contain both founder models.")
    if counts.nunique() != 1:
        raise RuntimeError(f"Unequal simulation counts: {counts.to_dict()}")
    number = int(counts.iloc[0])
    schedule = model160.simulation_schedule(number)
    lookup = simulations.set_index(["true_founder_model", "replicate"])
    donor_lookup = donor.set_index("sample")
    rng = np.random.default_rng(model160.SEED)
    metric_rows = []
    comparison_arrays = {}

    for founder in ["three_founder", "four_founder"]:
        for replicate, (c_left, c_right) in enumerate(schedule, start=1):
            true_left, true_right, donor_a, donor_b, construction = model160.true_genotypes(
                founder, matrix, included, allele_frequency, replicate, rng
            )
            simulated = model160.simulate_arrays(
                model159, true_left, true_right, c_left, c_right,
                depths, emission, rng
            )
            saved = lookup.loc[(founder, replicate)]
            if str(saved["donor_a"]) != donor_a or str(saved["donor_b"]) != donor_b:
                raise RuntimeError(
                    f"Regeneration donor mismatch for {founder} replicate {replicate}."
                )
            if not np.isclose(float(saved["true_c_left"]), c_left) or not np.isclose(
                float(saved["true_c_right"]), c_right
            ):
                raise RuntimeError(
                    f"Regeneration mixing mismatch for {founder} replicate {replicate}."
                )

            dosage_difference = np.abs(true_left - true_right)
            opposing = ((true_left == 0) & (true_right == 2)) | (
                (true_left == 2) & (true_right == 0)
            )
            opposing_fraction = float(opposing.mean())
            if not np.isclose(
                opposing_fraction,
                float(saved["true_opposing_homozygote_fraction"]),
                atol=1e-10,
            ):
                raise RuntimeError(
                    f"Regeneration genotype mismatch for {founder} replicate {replicate}."
                )

            metric_rows.append({
                "true_founder_model": founder,
                "replicate": replicate,
                "genotype_discordance_fraction": float((dosage_difference > 0).mean()),
                "mean_genotype_contrast": float((dosage_difference / 2).mean()),
                "opposing_homozygote_fraction": opposing_fraction,
                "minimum_donor_call_rate": float(min(
                    donor_lookup.loc[donor_a, "call_rate"],
                    donor_lookup.loc[donor_b, "call_rate"],
                )),
            })
            if np.isclose(c_left, 0.20) and np.isclose(c_right, 0.20):
                comparison_arrays[(founder, replicate)] = {
                    "arrays": simulated,
                    "true_left": true_left,
                    "true_right": true_right,
                }

    metrics = pd.DataFrame(metric_rows)
    augmented = simulations.merge(
        metrics,
        on=["true_founder_model", "replicate"],
        how="left",
        validate="one_to_one",
    )
    augmented["either_estimate_at_zero"] = False
    for founder, prefix in [
        ("three_founder", "null"),
        ("four_founder", "alternative"),
    ]:
        mask = augmented["true_founder_model"].eq(founder)
        augmented.loc[mask, "mean_absolute_recovery_error"] = 0.5 * (
            np.abs(augmented.loc[mask, f"{prefix}_c_left"] - augmented.loc[mask, "true_c_left"])
            + np.abs(augmented.loc[mask, f"{prefix}_c_right"] - augmented.loc[mask, "true_c_right"])
        )
        augmented.loc[mask, "either_estimate_at_zero"] = (
            augmented.loc[mask, f"{prefix}_c_left"].le(1e-8)
            | augmented.loc[mask, f"{prefix}_c_right"].le(1e-8)
        )
    augmented["either_estimate_at_zero"] = augmented[
        "either_estimate_at_zero"
    ].astype(bool)
    return augmented, comparison_arrays, matrix, donor, allele_frequency, depths, emission


def representative_simulation(augmented, founder):
    subset = augmented.loc[
        augmented["true_founder_model"].eq(founder)
        & np.isclose(augmented["true_c_left"], 0.20)
        & np.isclose(augmented["true_c_right"], 0.20)
        & boolean_series(augmented["both_donors_high_call"])
    ].copy()
    if subset.empty:
        subset = augmented.loc[
            augmented["true_founder_model"].eq(founder)
            & np.isclose(augmented["true_c_left"], 0.20)
            & np.isclose(augmented["true_c_right"], 0.20)
        ].copy()
    if subset.empty:
        raise RuntimeError(f"No 20%/20% {founder} simulation is available.")
    target = subset["mean_absolute_recovery_error"].median()
    subset["distance_from_median_error"] = np.abs(
        subset["mean_absolute_recovery_error"] - target
    )
    return subset.sort_values(
        ["distance_from_median_error", "replicate"]
    ).iloc[0]


def make_unmixed_control(model159, model160, matrix, donor, allele_frequency, depths, emission):
    candidates = donor.loc[donor["high_call_sensitivity_set"].astype(bool)].copy()
    target = candidates["heterozygosity"].median()
    candidates["distance_from_median_heterozygosity"] = np.abs(
        candidates["heterozygosity"] - target
    )
    selected = candidates.sort_values(
        ["distance_from_median_heterozygosity", "sample"]
    ).iloc[0]
    rng = np.random.default_rng(model160.SEED + 16_200)
    dosage = model160.imputed_dosage(
        matrix, str(selected["sample"]), allele_frequency, rng
    )
    arrays = model160.simulate_arrays(
        model159, dosage, dosage, 0.0, 0.0, depths, emission, rng
    )
    return selected, arrays


def fraction_arrays(arrays):
    return (
        arrays["left_alt"] / arrays["left_total"],
        arrays["right_alt"] / arrays["right_total"],
    )


def expected_centroids(model159, emission, c_left, c_right):
    expected_left = model159.map_fraction(
        (1 - c_left) * model159.GL + c_left * model159.GR,
        emission["error_rate"], emission["bias_log_odds"],
    )
    expected_right = model159.map_fraction(
        c_right * model159.GL + (1 - c_right) * model159.GR,
        emission["error_rate"], emission["bias_log_odds"],
    )
    return expected_left, expected_right


def common_hexbin_max(panels):
    figure, axis = plt.subplots()
    maximum = 1.0
    for panel in panels:
        artist = axis.hexbin(
            panel["x"], panel["y"], gridsize=43, mincnt=1,
            extent=(0, 1, 0, 1), linewidths=0,
        )
        if len(artist.get_array()):
            maximum = max(maximum, float(np.max(artist.get_array())))
        axis.clear()
    plt.close(figure)
    return maximum


def add_centroids(axis, model159, emission, c_left, c_right, model_type):
    left, right = expected_centroids(model159, emission, c_left, c_right)
    compatible = np.array([0, 1, 3, 4, 5, 7, 8])
    opposing = np.array([2, 6])
    diagonal = np.array([0, 4, 8])
    if model_type == "unmixed":
        axis.scatter(
            left[diagonal], right[diagonal], marker="s", s=50,
            facecolor="white", edgecolor=GREEN, linewidth=1.4, zorder=5,
        )
    else:
        axis.scatter(
            left[compatible], right[compatible], marker="s", s=45,
            facecolor="white", edgecolor=BLUE, linewidth=1.2, zorder=5,
        )
        if model_type in {"four_founder", "observed"}:
            axis.scatter(
                left[opposing], right[opposing], marker="D", s=55,
                facecolor=ORANGE, edgecolor="white", linewidth=0.8, zorder=6,
            )


def plot_step2(panels, model159, emission, output_pdf, output_png=None):
    maximum = common_hexbin_max(panels)
    norm = LogNorm(vmin=1, vmax=max(2, maximum))
    figure, axes = plt.subplots(2, 2, figsize=(11.5, 9.2), sharex=True, sharey=True)
    figure.subplots_adjust(
        top=0.86, bottom=0.12, left=0.08, right=0.90,
        hspace=0.30, wspace=0.18,
    )
    figure.suptitle(
        "Observed Left/Right read fractions against matched reference scenarios",
        fontsize=16, weight="bold", y=0.97,
    )
    figure.text(
        0.08, 0.91,
        "Every panel uses the focal depth profile, identical axes and one shared log-density scale.",
        fontsize=9.5, color="#4A4A4A",
    )
    letters = "ABCD"
    artist = None
    for letter, axis, panel in zip(letters, axes.ravel(), panels):
        artist = axis.hexbin(
            panel["x"], panel["y"], gridsize=43, mincnt=1,
            extent=(0, 1, 0, 1), cmap="Greys", norm=norm, linewidths=0,
        )
        add_centroids(
            axis, model159, emission, panel["c_left"], panel["c_right"],
            panel["model_type"],
        )
        axis.set_xlim(-0.02, 1.02)
        axis.set_ylim(-0.02, 1.02)
        axis.set_aspect("equal", adjustable="box")
        axis.set_title(f"{letter}  {panel['title']}", loc="left", weight="bold")
        axis.text(
            0.02, 0.98, panel["subtitle"], transform=axis.transAxes,
            ha="left", va="top", fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#BBBBBB", alpha=0.9),
        )
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
    axes[1, 0].set_xlabel("Observed alternate-read fraction: Left")
    axes[1, 1].set_xlabel("Observed alternate-read fraction: Left")
    axes[0, 0].set_ylabel("Observed alternate-read fraction: Right")
    axes[1, 0].set_ylabel("Observed alternate-read fraction: Right")
    colorbar_axis = figure.add_axes([0.925, 0.19, 0.018, 0.60])
    figure.colorbar(artist, cax=colorbar_axis, label="Sites per hexagon (log scale)")
    figure.text(
        0.08, 0.045,
        "Green squares: same-diploid states. Blue squares: three-founder-compatible states. "
        "Orange diamonds: opposing homozygotes available only without an obligatory shared gamete.",
        fontsize=8.2,
    )
    figure.savefig(output_pdf, bbox_inches="tight")
    if output_png is not None:
        figure.savefig(output_png, dpi=300, bbox_inches="tight")
    return figure


def calibration_panel(axis, data, mixing, founder, side, contrast_norm):
    prefix = "null" if founder == "three_founder" else "alternative"
    color = BLUE if founder == "three_founder" else ORANGE
    subset = data.loc[data["true_founder_model"].eq(founder)].copy()
    high_call_mask = boolean_series(subset["both_donors_high_call"])
    primary = subset.loc[high_call_mask].copy()
    background = subset.loc[~high_call_mask].copy()
    x_column = f"true_c_{side}"
    y_column = f"{prefix}_c_{side}"
    axis.scatter(
        background[x_column], background[y_column], s=12,
        color="#C7C7C7", alpha=0.35, linewidth=0, label="lower-call pairs",
    )
    artist = axis.scatter(
        primary[x_column], primary[y_column], s=22,
        c=primary["mean_genotype_contrast"], cmap=CMAP, norm=contrast_norm,
        alpha=0.75, edgecolor="white", linewidth=0.25, label="high-call pairs",
    )
    grouped = primary.groupby(x_column)[y_column]
    x = np.array(sorted(primary[x_column].unique()))
    median = grouped.median().reindex(x).to_numpy()
    low = grouped.quantile(0.025).reindex(x).to_numpy()
    high = grouped.quantile(0.975).reindex(x).to_numpy()
    axis.fill_between(x, low, high, color=color, alpha=0.13, linewidth=0)
    axis.plot(x, median, color=color, linewidth=2, label="high-call median")
    axis.plot([0, 0.42], [0, 0.42], color=GREY, linestyle="--", linewidth=1)
    boundary = primary.loc[primary[y_column].le(1e-8)]
    axis.scatter(
        boundary[x_column], boundary[y_column], marker="x", s=42,
        color=BLACK, linewidth=1.2, label="estimate at zero",
    )
    row = mixing.loc[
        mixing["founder_model"].eq(founder) & mixing["side"].eq(side)
    ].iloc[0]
    axis.scatter(
        row["calibrated_estimate"], row["observed_raw_estimate"],
        marker="*", s=120, color=BLACK, zorder=7, label="observed projection",
    )
    axis.set_xlim(-0.01, 0.42)
    axis.set_ylim(-0.01, 0.46)
    axis.set_xlabel(f"True mixing into {side.title()}")
    axis.set_ylabel(f"Recovered mixing into {side.title()}")
    axis.set_title(f"{founder.replace('_', ' ').title()}: {side.title()}")
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.grid(axis="y", color="#ECECEC", linewidth=0.6, zorder=0)
    axis.text(
        0.98, 0.04,
        f"High-call: {len(primary)}\nLower-call: {len(background)}\n"
        f"Zero-boundary: {len(boundary)}",
        transform=axis.transAxes, ha="right", va="bottom", fontsize=7,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#BBBBBB", alpha=0.9),
    )
    return artist


def plot_step3_calibration(data, mixing, figure):
    axes = figure.subplots(2, 2)
    figure.subplots_adjust(
        top=0.82, bottom=0.12, left=0.07, right=0.90,
        hspace=0.38, wspace=0.28,
    )
    figure.suptitle(
        "Mixing recovery depends on between-side genotype information",
        fontsize=16, weight="bold", y=0.97,
    )
    figure.text(
        0.07, 0.91,
        "High-call simulations are coloured by mean genotype contrast; lower-call pairs remain as grey context.",
        fontsize=9.5, color="#4A4A4A",
    )
    high_call = data.loc[boolean_series(data["both_donors_high_call"])]
    contrast_norm = Normalize(
        vmin=float(high_call["mean_genotype_contrast"].min()),
        vmax=float(high_call["mean_genotype_contrast"].max()),
    )
    artist = None
    combinations = [
        ("three_founder", "left"), ("three_founder", "right"),
        ("four_founder", "left"), ("four_founder", "right"),
    ]
    for axis, (founder, side) in zip(axes.ravel(), combinations):
        artist = calibration_panel(
            axis, data, mixing, founder, side, contrast_norm
        )
    axes[0, 0].legend(frameon=False, fontsize=6.5, loc="upper left")
    colorbar_axis = figure.add_axes([0.925, 0.20, 0.018, 0.56])
    figure.colorbar(
        artist, cax=colorbar_axis,
        label="Mean genotype contrast between sides",
    )
    figure.text(
        0.07, 0.045,
        "A black x identifies a high-call simulation whose fitted mixing reached the zero boundary. "
        "The black star is the observed inverse-calibration projection, not a simulated replicate.",
        fontsize=8.2,
    )
    return figure


def plot_information_diagnostic(data, figure):
    axes = figure.subplots(1, 2, sharey=True)
    figure.subplots_adjust(top=0.76, bottom=0.20, left=0.08, right=0.89, wspace=0.26)
    figure.suptitle(
        "Why do some four-founder simulations recover mixing poorly?",
        fontsize=16, weight="bold", y=0.96,
    )
    figure.text(
        0.08, 0.875,
        "Two measures of between-side information are tested; colour records the weaker donor's call rate.",
        fontsize=9.5, color="#4A4A4A",
    )
    call_norm = Normalize(
        vmin=float(data["minimum_donor_call_rate"].min()),
        vmax=float(data["minimum_donor_call_rate"].max()),
    )
    artist = None
    subset = data.loc[data["true_founder_model"].eq("four_founder")].copy()
    diagnostics = [
        ("mean_genotype_contrast", "Genotype contrast"),
        ("opposing_homozygote_fraction", "Opposing-homozygote fraction"),
    ]
    for letter, axis, (x_column, x_label) in zip("AB", axes, diagnostics):
        artist = axis.scatter(
            subset[x_column],
            subset["mean_absolute_recovery_error"],
            c=subset["minimum_donor_call_rate"], cmap=CMAP, norm=call_norm,
            s=25, alpha=0.72, edgecolor="white", linewidth=0.25,
        )
        boundary = subset.loc[subset["either_estimate_at_zero"]]
        axis.scatter(
            boundary[x_column],
            boundary["mean_absolute_recovery_error"],
            marker="x", s=45, color=BLACK, linewidth=1.1,
            label="at least one estimate at zero",
        )
        correlation, p_value = spearmanr(
            subset[x_column],
            subset["mean_absolute_recovery_error"],
        )
        axis.set_title(f"{letter}  {x_label}", loc="left", color=ORANGE, weight="bold")
        axis.set_xlabel(
            "Mean genotype contrast between true sides"
            if x_column == "mean_genotype_contrast"
            else "True opposing-homozygote fraction"
        )
        axis.text(
            0.97, 0.62,
            f"Spearman rho = {correlation:.2f}\np = {p_value:.3g}\n"
            f"Zero-boundary pairs = {len(boundary)}",
            transform=axis.transAxes, ha="right", va="top", fontsize=7.5,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#BBBBBB", alpha=0.9),
        )
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(axis="y", color="#ECECEC", linewidth=0.6, zorder=0)
    axes[0].set_ylabel("Mean absolute mixing-recovery error")
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles, labels, frameon=False, fontsize=7,
        loc="lower left", bbox_to_anchor=(0.075, 0.075),
    )
    colorbar_axis = figure.add_axes([0.92, 0.25, 0.018, 0.42])
    figure.colorbar(artist, cax=colorbar_axis, label="Minimum donor call rate")
    return figure


def write_full_narrative(
    narrative, site, block, observed, simulations, mixing, bootstrap,
    summary, emission, model159, panels,
):
    observed_lr = float(summary["observed_nested_LR_statistic"])
    threshold = float(summary["three_founder_simulated_95_percent_threshold"])
    high_threshold = float(
        summary["high_call_three_founder_simulated_95_percent_threshold"]
    )
    intervals = narrative.interval_rows(mixing, bootstrap)
    with PdfPages(NARRATIVE_V2) as pdf, PdfPages(STEP3_PDF) as step3_pdf:
        figure = plt.figure(figsize=(11.7, 8.3))
        narrative.add_page_heading(
            figure, 1, "Define the biological alternatives",
            "Founder sharing and sequencing mixture are distinct parts of the model.",
        )
        axis = figure.add_axes([0.05, 0.08, 0.90, 0.74])
        narrative.model_schematic(axis)
        pdf.savefig(figure, bbox_inches="tight")
        plt.close(figure)

        figure = plot_step2(panels, model159, emission, STEP2_PDF, STEP2_PNG)
        pdf.savefig(figure, bbox_inches="tight")
        plt.close(figure)

        figure = plt.figure(figsize=(11.7, 8.3))
        plot_step3_calibration(simulations, mixing, figure)
        pdf.savefig(figure, bbox_inches="tight")
        step3_pdf.savefig(figure, bbox_inches="tight")
        figure.savefig(STEP3_PNG, dpi=300, bbox_inches="tight")
        plt.close(figure)

        figure = plt.figure(figsize=(11.7, 6.3))
        plot_information_diagnostic(simulations, figure)
        pdf.savefig(figure, bbox_inches="tight")
        step3_pdf.savefig(figure, bbox_inches="tight")
        figure.savefig(STEP3_CAUSE_PNG, dpi=300, bbox_inches="tight")
        plt.close(figure)

        figure, axes = plt.subplots(1, 2, figsize=(11.7, 6.3))
        figure.subplots_adjust(top=0.72, left=0.08, right=0.97, bottom=0.14, wspace=0.25)
        narrative.add_page_heading(
            figure, 4, "Calibrate evidence for four founders",
            "The observed statistic is compared only with simulations near the observed mixing.",
        )
        narrative.plot_founder_panel(axes[0], simulations, observed_lr, threshold, False)
        narrative.plot_founder_panel(axes[1], simulations, observed_lr, high_threshold, True)
        pdf.savefig(figure, bbox_inches="tight")
        plt.close(figure)

        figure, axes = plt.subplots(1, 2, figsize=(11.7, 6.3))
        figure.subplots_adjust(top=0.72, left=0.15, right=0.97, bottom=0.14, wspace=0.32)
        narrative.add_page_heading(
            figure, 5, "Report reciprocal mixing with uncertainty",
            "Thick intervals resample genomic blocks; thin intervals include simulation recovery.",
        )
        narrative.plot_interval_panel(axes[0], intervals)
        narrative.plot_bootstrap_panel(axes[1], bootstrap)
        pdf.savefig(figure, bbox_inches="tight")
        plt.close(figure)

        figure, axes = plt.subplots(2, 1, figsize=(11.7, 7.5), sharex=True)
        figure.subplots_adjust(top=0.76, left=0.10, right=0.97, bottom=0.10, hspace=0.16)
        narrative.add_page_heading(
            figure, 6, "Locate the founder-model signal along Hmr18",
            "Coherent multi-block signal is more persuasive than isolated site-level deviations.",
        )
        narrative.plot_block_panel(axes[0], block, False)
        axes[0].set_xlabel("")
        axes[0].set_title(
            f"Fixed-fit block contributions; total LR = "
            f"{block['nested_LR_contribution'].sum():.2f}"
        )
        narrative.plot_block_panel(axes[1], block, True)
        axes[1].legend(frameon=False, loc="upper right")
        pdf.savefig(figure, bbox_inches="tight")
        plt.close(figure)


def main():
    required = [
        SCRIPT159, SCRIPT160, SCRIPT161, MATCHED_FILE, OBSERVED_FILE,
        SIMULATION_FILE, MIXING_FILE, BOOTSTRAP_FILE, SUMMARY_FILE,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise RuntimeError("Missing required inputs:\n" + "\n".join(missing))
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    model159 = load_module(SCRIPT159, "model159_for_stage162")
    model160 = load_module(SCRIPT160, "model160_for_stage162")
    narrative = load_module(SCRIPT161, "narrative161_for_stage162")
    narrative.set_style()

    matched = pd.read_csv(MATCHED_FILE, sep="\t", low_memory=False)
    observed = pd.read_csv(OBSERVED_FILE, sep="\t")
    simulations = pd.read_csv(SIMULATION_FILE, sep="\t", low_memory=False)
    mixing = pd.read_csv(MIXING_FILE, sep="\t")
    bootstrap = pd.read_csv(BOOTSTRAP_FILE, sep="\t")
    summary = summary_dictionary(SUMMARY_FILE)
    site, block, _ = narrative.calculate_site_evidence(
        matched, observed, model159, model160
    )
    augmented, comparison_arrays, matrix, donor, allele_frequency, depths, emission = (
        regenerate_simulation_information(model159, model160, matched, simulations)
    )
    augmented.to_csv(METRICS_FILE, sep="\t", index=False, compression="gzip")

    m3 = representative_simulation(augmented, "three_founder")
    m4 = representative_simulation(augmented, "four_founder")
    unmixed_donor, unmixed_arrays = make_unmixed_control(
        model159, model160, matrix, donor, allele_frequency, depths, emission
    )
    m3_arrays = comparison_arrays[("three_founder", int(m3["replicate"]))]["arrays"]
    m4_arrays = comparison_arrays[("four_founder", int(m4["replicate"]))]["arrays"]
    observed_alt = observed.loc[
        observed["fit"].eq("opposing_state_weight_alternative")
    ].iloc[0]
    observed_x = matched["left_alt_depth"] / matched["left_AD_total"]
    observed_y = matched["right_alt_depth"] / matched["right_AD_total"]
    unmixed_x, unmixed_y = fraction_arrays(unmixed_arrays)
    m3_x, m3_y = fraction_arrays(m3_arrays)
    m4_x, m4_y = fraction_arrays(m4_arrays)
    panels = [
        {
            "title": "Observed gynandromorph",
            "subtitle": f"Real focal allele depths; fitted cL={float(observed_alt['c_left']):.3f}, "
                        f"cR={float(observed_alt['c_right']):.3f}",
            "x": observed_x, "y": observed_y,
            "c_left": float(observed_alt["c_left"]),
            "c_right": float(observed_alt["c_right"]), "model_type": "observed",
        },
        {
            "title": "Unmixed diploid reference",
            "subtitle": f"Real genotype: {unmixed_donor['sample']}\n"
                        "Independent pseudo-read replicates; c=0",
            "x": unmixed_x, "y": unmixed_y,
            "c_left": 0.0, "c_right": 0.0, "model_type": "unmixed",
        },
        {
            "title": "Simulated three-founder individual",
            "subtitle": f"Representative replicate {int(m3['replicate'])}; "
                        f"{m3['donor_a']} + {m3['donor_b']}\nTrue cL=cR=0.20",
            "x": m3_x, "y": m3_y,
            "c_left": 0.20, "c_right": 0.20, "model_type": "three_founder",
        },
        {
            "title": "Simulated four-founder individual",
            "subtitle": f"Representative replicate {int(m4['replicate'])}; "
                        f"{m4['donor_a']} + {m4['donor_b']}\nTrue cL=cR=0.20",
            "x": m4_x, "y": m4_y,
            "c_left": 0.20, "c_right": 0.20, "model_type": "four_founder",
        },
    ]

    representative = pd.DataFrame([
        {
            "scenario": "unmixed_diploid",
            "sample_or_replicate": str(unmixed_donor["sample"]),
            "selection_rule": "high-call control nearest median heterozygosity",
            "mean_recovery_error": np.nan,
            "mean_genotype_contrast": 0.0,
        },
        {
            "scenario": "three_founder_20_percent",
            "sample_or_replicate": int(m3["replicate"]),
            "selection_rule": "high-call 20%/20% simulation nearest median recovery error",
            "mean_recovery_error": m3["mean_absolute_recovery_error"],
            "mean_genotype_contrast": m3["mean_genotype_contrast"],
        },
        {
            "scenario": "four_founder_20_percent",
            "sample_or_replicate": int(m4["replicate"]),
            "selection_rule": "high-call 20%/20% simulation nearest median recovery error",
            "mean_recovery_error": m4["mean_absolute_recovery_error"],
            "mean_genotype_contrast": m4["mean_genotype_contrast"],
        },
    ])
    representative.to_csv(REPRESENTATIVE_FILE, sep="\t", index=False)
    write_full_narrative(
        narrative, site, block, observed, augmented, mixing, bootstrap,
        summary, emission, model159, panels,
    )

    audit = pd.DataFrame([
        ("matched_sites", len(matched)),
        ("simulations", len(augmented)),
        ("three_founder_simulations", augmented["true_founder_model"].eq("three_founder").sum()),
        ("four_founder_simulations", augmented["true_founder_model"].eq("four_founder").sum()),
        ("high_call_simulations", boolean_series(augmented["both_donors_high_call"]).sum()),
        ("three_founder_zero_boundary", augmented.loc[
            augmented["true_founder_model"].eq("three_founder"), "either_estimate_at_zero"
        ].sum()),
        ("four_founder_zero_boundary", augmented.loc[
            augmented["true_founder_model"].eq("four_founder"), "either_estimate_at_zero"
        ].sum()),
        ("unmixed_reference_sample", unmixed_donor["sample"]),
        ("representative_three_founder_replicate", int(m3["replicate"])),
        ("representative_four_founder_replicate", int(m4["replicate"])),
        ("focal_GT_used", False),
        ("old_T2_T3_used", False),
    ], columns=["item", "value"])
    audit.to_csv(AUDIT_FILE, sep="\t", index=False)
    outputs = [
        METRICS_FILE, REPRESENTATIVE_FILE, STEP2_PDF, STEP2_PNG,
        STEP3_PDF, STEP3_PNG, STEP3_CAUSE_PNG, NARRATIVE_V2, AUDIT_FILE,
    ]
    print("STAGE 162 VISUALIZATION COMPLETE")
    print(audit.to_string(index=False))
    for path in outputs:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty output: {path}")
        print(f"PRESENT\t{path}\t{path.stat().st_size} bytes")


if __name__ == "__main__":
    main()
