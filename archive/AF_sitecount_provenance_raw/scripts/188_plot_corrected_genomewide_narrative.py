import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm

ROOT = Path("analysis/genomewide_mixing_corrected_v1")
BASE = ROOT / "expanded_calibration"
V2 = ROOT / "expanded_calibration_v2"
FIGURE_DIR = V2 / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

BLUE = "#0072B2"
LIGHT_BLUE = "#56B4E9"
GREEN = "#009E73"
ORANGE = "#D55E00"
EDGES = np.linspace(0, 1, 56)


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pilot = load_module(
    "scripts/185_restored_model_pilot.py",
    "pilot185",
)
stage = load_module(
    "scripts/187_expand_simulations_and_window_evidence.py",
    "stage187",
)
pilot.OUT = V2


def density(x, y):
    return np.histogram2d(
        x,
        y,
        bins=[EDGES, EDGES],
    )[0].T


def draw_density(axis, values, norm, title):
    artist = axis.pcolormesh(
        EDGES,
        EDGES,
        np.ma.masked_where(values == 0, values),
        cmap="Greys",
        norm=norm,
        shading="flat",
    )
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.set_aspect("equal")
    axis.set_xlabel("Left ALT-read fraction")
    axis.set_title(
        title,
        loc="left",
        weight="bold",
        fontsize=10,
    )
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    return artist


def add_centroids(axis, mixing, show_opposing=True):
    emission = pilot.load_emission()

    expected_left = (
        (1 - mixing["c_left"]) * pilot.model.GL
        + mixing["c_left"] * pilot.model.GR
    )
    expected_right = (
        mixing["c_right"] * pilot.model.GL
        + (1 - mixing["c_right"]) * pilot.model.GR
    )

    x = pilot.model.map_fraction(
        expected_left,
        emission["error_rate"],
        emission["bias_log_odds"],
    )
    y = pilot.model.map_fraction(
        expected_right,
        emission["error_rate"],
        emission["bias_log_odds"],
    )

    axis.scatter(
        x[pilot.COMPATIBLE],
        y[pilot.COMPATIBLE],
        marker="s",
        s=34,
        facecolor="white",
        edgecolor=BLUE,
        linewidth=1.1,
        zorder=4,
    )

    if show_opposing:
        axis.scatter(
            x[pilot.OPPOSING],
            y[pilot.OPPOSING],
            marker="D",
            s=36,
            facecolor=ORANGE,
            edgecolor="white",
            linewidth=0.7,
            zorder=5,
        )


def full_benchmarks(simulations):
    data = stage.load_data(calibration=False)
    matrix = data["controls"]
    called = matrix >= 0

    allele_frequency = (
        np.where(called, matrix, 0).sum(axis=1) + 0.5
    ) / (
        2 * called.sum(axis=1) + 1
    )

    observed = json.loads(
        (BASE / "observed_expanded.json").read_text()
    )

    panels = [(
        "Observed L/R",
        data["left_alt"] / data["left_total"],
        data["right_alt"] / data["right_total"],
        "observed",
        observed["full"]["alternative"],
    )]

    for founder_index, founder in enumerate([
        "three",
        "four",
    ]):
        subset = simulations.loc[
            (simulations["family"] == "repeat_20")
            & (
                simulations["generating_model"]
                == founder
            )
        ].copy()

        subset["recovery_error"] = (
            (subset["estimated_c_left_alt"] - 0.2).abs()
            + (
                subset["estimated_c_right_alt"] - 0.2
            ).abs()
        )

        subset = subset.sort_values("recovery_error")
        row = subset.iloc[len(subset) // 2]

        identifier = int(row["id"])
        donor_a_index = stage.CONTROLS.index(
            row["donor_a"]
        )
        donor_b_index = stage.CONTROLS.index(
            row["donor_b"]
        )

        rng = np.random.default_rng(
            187000000
            + identifier * 2
            + founder_index
        )

        donor_a = stage.impute(
            matrix,
            donor_a_index,
            allele_frequency,
            rng,
        )
        donor_b = stage.impute(
            matrix,
            donor_b_index,
            allele_frequency,
            rng,
        )

        if founder == "three":
            constructed = (
                stage.gamete(donor_a, rng)
                + stage.gamete(donor_b, rng)
            )

            if identifier % 2 == 0:
                true_left, true_right = (
                    donor_a,
                    constructed,
                )
            else:
                true_left, true_right = (
                    constructed,
                    donor_a,
                )
        else:
            true_left, true_right = donor_a, donor_b

        synthetic = pilot.simulate_counts(
            true_left,
            true_right,
            0.2,
            0.2,
            data,
            pilot.load_emission(),
            rng,
        )

        panels.append((
            f"Simulated {founder}-founder",
            synthetic["left_alt"]
            / synthetic["left_total"],
            synthetic["right_alt"]
            / synthetic["right_total"],
            founder,
            {"c_left": 0.2, "c_right": 0.2},
        ))

    return panels


def plot_equal_scope(pdf, simulations):
    panels = full_benchmarks(simulations)
    densities = [
        density(panel[1], panel[2])
        for panel in panels
    ]

    norm = LogNorm(
        vmin=1,
        vmax=max(values.max() for values in densities),
    )

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(11.7, 5),
        sharex=True,
        sharey=True,
    )
    figure.subplots_adjust(
        left=0.07,
        right=0.90,
        bottom=0.17,
        top=0.77,
        wspace=0.16,
    )

    artist = None

    for axis, panel, values in zip(
        axes,
        panels,
        densities,
    ):
        artist = draw_density(
            axis,
            values,
            norm,
            panel[0],
        )
        add_centroids(
            axis,
            panel[4],
            show_opposing=panel[3] != "three",
        )

    axes[0].set_ylabel("Right ALT-read fraction")

    figure.suptitle(
        "1  Equal-scope allele-fraction patterns",
        x=0.07,
        ha="left",
        fontsize=15,
        weight="bold",
    )
    figure.text(
        0.07,
        0.84,
        "Each panel contains exactly 1,000,000 sites "
        "with identical positions and depth profiles. "
        "One shared log-density scale.",
        fontsize=9,
        color="#444444",
    )

    colorbar_axis = figure.add_axes([
        0.925, 0.23, 0.015, 0.43
    ])
    figure.colorbar(
        artist,
        cax=colorbar_axis,
        label="Sites per hexagon",
    )

    pdf.savefig(figure, bbox_inches="tight")
    plt.close(figure)


def plot_recovery(pdf, simulations):
    frame = simulations.loc[
        simulations["family"].isin([
            "global",
            "focused",
        ])
    ].copy()

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(11.7, 8.1),
        sharex="col",
        sharey="row",
    )
    figure.subplots_adjust(
        left=0.08,
        right=0.97,
        bottom=0.10,
        top=0.84,
        hspace=0.27,
        wspace=0.20,
    )

    for row_index, founder in enumerate([
        "three",
        "four",
    ]):
        subset = frame.loc[
            frame["generating_model"] == founder
        ].copy()

        subset["contrast_bin"] = pd.qcut(
            subset["genotype_contrast"],
            3,
            labels=[
                "low contrast",
                "middle contrast",
                "high contrast",
            ],
            duplicates="drop",
        )

        definitions = [
            (
                "true_c_left",
                "estimated_c_left_alt",
                "Left",
            ),
            (
                "true_c_right",
                "estimated_c_right_alt",
                "Right",
            ),
        ]

        for column_index, (
            true_column,
            estimate_column,
            side,
        ) in enumerate(definitions):
            axis = axes[row_index, column_index]

            axis.scatter(
                subset[true_column],
                subset[estimate_column],
                s=5,
                color="#BBBBBB",
                alpha=0.13,
                rasterized=True,
            )

            for label, color in zip(
                [
                    "low contrast",
                    "middle contrast",
                    "high contrast",
                ],
                [
                    LIGHT_BLUE,
                    GREEN,
                    ORANGE,
                ],
            ):
                selected = subset.loc[
                    subset["contrast_bin"] == label
                ]

                trend = selected.groupby(
                    true_column,
                    as_index=False,
                )[estimate_column].median()

                axis.plot(
                    trend[true_column],
                    trend[estimate_column],
                    marker="o",
                    markersize=3,
                    linewidth=1.5,
                    color=color,
                    label=label,
                )

            axis.plot(
                [0, 0.45],
                [0, 0.45],
                linestyle="--",
                linewidth=1,
                color="#555555",
            )
            axis.set_xlim(-0.01, 0.46)
            axis.set_ylim(-0.01, 0.46)
            axis.set_xlabel(f"True {side} mixing")
            axis.set_ylabel(
                f"Estimated {side} mixing"
            )
            axis.set_title(
                f"{founder.capitalize()}-founder: {side}",
                loc="left",
                weight="bold",
                fontsize=10,
            )
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)

    axes[0, 1].legend(
        frameon=False,
        fontsize=8,
    )

    figure.suptitle(
        "2  Mixing recovery across genotype contrast",
        x=0.08,
        ha="left",
        fontsize=15,
        weight="bold",
    )
    figure.text(
        0.08,
        0.89,
        "Points are simulations; coloured lines are "
        "median recovery trends within genotype-contrast thirds.",
        fontsize=9,
        color="#444444",
    )

    pdf.savefig(figure, bbox_inches="tight")
    plt.close(figure)


def plot_lr_distributions(pdf, simulations):
    observed_lr = json.loads(
        (BASE / "observed_expanded.json").read_text()
    )["calibration"]["LR"]

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(11.7, 4.8),
    )
    figure.subplots_adjust(
        left=0.08,
        right=0.97,
        bottom=0.17,
        top=0.78,
        wspace=0.23,
    )

    for axis, family, title in zip(
        axes,
        ["global", "focused"],
        [
            "Global 0-45% grid",
            "Focused 14-26% grid",
        ],
    ):
        subset = simulations.loc[
            simulations["family"] == family
        ]

        for founder, color in [
            ("three", BLUE),
            ("four", ORANGE),
        ]:
            values = np.log10(
                1
                + subset.loc[
                    subset["generating_model"] == founder,
                    "LR",
                ]
            )
            axis.hist(
                values,
                bins=45,
                density=True,
                histtype="step",
                linewidth=1.8,
                color=color,
                label=f"{founder}-founder",
            )

        axis.axvline(
            np.log10(1 + observed_lr),
            color="black",
            linewidth=1.5,
            label="observed 20k-panel LR",
        )
        axis.set_xlabel("log10(1 + nested LR)")
        axis.set_ylabel("Density")
        axis.set_title(title)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    axes[0].legend(
        frameon=False,
        fontsize=8,
    )

    figure.suptitle(
        "3  Simulation-calibrated founder evidence",
        x=0.08,
        ha="left",
        fontsize=15,
        weight="bold",
    )
    figure.text(
        0.08,
        0.86,
        "Observed and simulated results use the identical "
        "20,000-site chromosome-balanced panel.",
        fontsize=9,
        color="#444444",
    )

    pdf.savefig(figure, bbox_inches="tight")
    plt.close(figure)


def plot_manhattan(pdf):
    windows = pd.read_csv(
        V2 / "autosomal_1Mb_evidence.tsv",
        sep="\t",
    )

    genome_positions = []
    centers = {}
    offset = 0

    for chromosome in range(1, 21):
        selected = windows["chromosome"] == chromosome
        local = (
            windows.loc[selected, "window_start"]
            + windows.loc[selected, "window_end"]
        ) / 2

        genome_positions.extend(local + offset)
        centers[chromosome] = (
            offset + local.max() / 2
        )
        offset += (
            windows.loc[selected, "window_end"].max()
            + 500000
        )

    windows["genome_position"] = genome_positions

    figure, axes = plt.subplots(
        2,
        1,
        figsize=(11.7, 6.4),
        sharex=True,
    )
    figure.subplots_adjust(
        left=0.08,
        right=0.97,
        bottom=0.12,
        top=0.82,
        hspace=0.12,
    )

    colors = [
        BLUE if int(chromosome) % 2 else LIGHT_BLUE
        for chromosome in windows["chromosome"]
    ]

    axes[0].scatter(
        windows["genome_position"],
        windows["LR_per_site"],
        c=colors,
        s=18,
    )
    axes[1].scatter(
        windows["genome_position"],
        windows["opposing_posterior_mean"],
        c=colors,
        s=18,
    )

    axes[0].axhline(
        0,
        color="#777777",
        linewidth=0.8,
    )
    axes[0].set_ylabel("LR contribution / site")
    axes[1].set_ylabel(
        "Mean opposing-state\nposterior"
    )
    axes[1].set_xticks(
        list(centers.values()),
        [str(value) for value in centers],
        fontsize=7,
    )
    axes[1].set_xlabel("Chromosome")

    for axis in axes:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    figure.suptitle(
        "4  Where does opposing-state evidence occur?",
        x=0.08,
        ha="left",
        fontsize=15,
        weight="bold",
    )
    figure.text(
        0.08,
        0.87,
        "One-megabase windows evaluated at the genome-wide "
        "parameters. These are additive contributions, "
        "not independent tests.",
        fontsize=9,
        color="#444444",
    )

    pdf.savefig(figure, bbox_inches="tight")
    plt.close(figure)


def plot_z_diagnostic(pdf):
    autosomes = stage.load_data(calibration=False)
    rng = np.random.default_rng(188)
    selected = rng.choice(
        len(autosomes["pos"]),
        50000,
        replace=False,
    )

    with np.load(
        V2 / "z_diagnostic.npz",
        allow_pickle=False,
    ) as source:
        z = {
            key: source[key]
            for key in source.files
        }

    panels = [
        (
            "Autosomal reference",
            (
                autosomes["left_alt"][selected]
                / autosomes["left_total"][selected]
            ),
            (
                autosomes["right_alt"][selected]
                / autosomes["right_total"][selected]
            ),
        ),
        (
            "Chromosome Z",
            z["left_alt"] / z["left_total"],
            z["right_alt"] / z["right_total"],
        ),
    ]

    densities = [
        density(x, y)
        for title, x, y in panels
    ]
    norm = LogNorm(
        vmin=1,
        vmax=max(values.max() for values in densities),
    )

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(9.4, 5),
        sharex=True,
        sharey=True,
    )
    figure.subplots_adjust(
        left=0.09,
        right=0.89,
        bottom=0.17,
        top=0.76,
        wspace=0.18,
    )

    artist = None

    for axis, panel, values in zip(
        axes,
        panels,
        densities,
    ):
        artist = draw_density(
            axis,
            values,
            norm,
            panel[0],
        )

    axes[0].set_ylabel("Right ALT-read fraction")

    colorbar_axis = figure.add_axes([
        0.92, 0.23, 0.018, 0.42
    ])
    figure.colorbar(
        artist,
        cax=colorbar_axis,
        label="Sites per hexagon",
    )

    audit = json.loads(
        (V2 / "z_diagnostic.audit.json").read_text()
    )

    figure.suptitle(
        "5  Z chromosome sex-ploidy diagnostic",
        x=0.09,
        ha="left",
        fontsize=15,
        weight="bold",
    )
    figure.text(
        0.09,
        0.83,
        "Equal 50,000-site panels. Mean Left:Right Z-depth "
        f"ratio = {audit['left_to_right_depth_ratio']:.3f}. "
        "Z is not fitted with the diploid founder model.",
        fontsize=9,
        color="#444444",
    )

    pdf.savefig(figure, bbox_inches="tight")
    plt.close(figure)


def main():
    simulations = pd.read_csv(
        V2 / "simulations.tsv.gz",
        sep="\t",
    )

    expected_rows = 7490
    if len(simulations) != expected_rows:
        raise RuntimeError(
            f"Expected {expected_rows} simulation rows; "
            f"found {len(simulations)}"
        )

    pdf_path = (
        FIGURE_DIR
        / "genomewide_founder_mixing_narrative_corrected.pdf"
    )

    with PdfPages(pdf_path) as pdf:
        plot_equal_scope(pdf, simulations)
        plot_recovery(pdf, simulations)
        plot_lr_distributions(pdf, simulations)
        plot_manhattan(pdf)
        plot_z_diagnostic(pdf)

    audit = {
        "pages": 5,
        "equal_scope_sites": 1000000,
        "simulation_rows": len(simulations),
        "window_size": 1000000,
        "posterior_model_probability_assigned": False,
        "pdf": str(pdf_path),
    }

    (
        FIGURE_DIR / "narrative.audit.json"
    ).write_text(
        json.dumps(audit, indent=2) + "\n"
    )

    print(json.dumps(audit, indent=2), flush=True)


if __name__ == "__main__":
    main()
