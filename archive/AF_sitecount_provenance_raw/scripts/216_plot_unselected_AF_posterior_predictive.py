#!/usr/bin/env python3

import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm, TwoSlopeNorm
import numpy as np
import pandas as pd


ROOT = Path(
    "analysis/genomewide_mixing_corrected_v1"
)

EMPIRICAL_SCRIPT = Path(
    "scripts/202_204_finish_empirical_sibling_calibration.py"
)

RESULTS_TABLE = (
    ROOT
    / "empirical_full_sib_calibration_v1"
    / "final"
    / "all_empirical_sibling_simulations.tsv.gz"
)

OBSERVED_FIT = (
    ROOT
    / "restored_pilot"
    / "observed.json"
)

OUTPUT_DIRECTORY = (
    ROOT
    / "representative_unselected_AF"
)

OUTPUT_PDF = (
    OUTPUT_DIRECTORY
    / "observed_vs_simulated_unselected_AF.pdf"
)

OUTPUT_PNG = (
    OUTPUT_DIRECTORY
    / "observed_vs_simulated_unselected_AF_page1.png"
)

OUTPUT_AUDIT = (
    OUTPUT_DIRECTORY
    / "observed_vs_simulated_unselected_AF.audit.json"
)

OUTPUT_ROWS = (
    OUTPUT_DIRECTORY
    / "selected_simulation_replicates.tsv"
)


BINS = 81

BLUE = "#0072B2"
ORANGE = "#D55E00"


REQUIRED_COUNTS = {
    "left_alt",
    "left_total",
    "right_alt",
    "right_total",
}


def load_module(path, name):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty module: {path}"
        )

    specification = (
        importlib.util.spec_from_file_location(
            name,
            path,
        )
    )

    module = (
        importlib.util.module_from_spec(
            specification
        )
    )

    specification.loader.exec_module(
        module
    )

    return module


def find_count_mapping(value):
    if isinstance(value, dict):
        if REQUIRED_COUNTS.issubset(value):
            return value

        for nested in value.values():
            try:
                return find_count_mapping(
                    nested
                )
            except RuntimeError:
                pass

    if isinstance(value, (list, tuple)):
        for nested in value:
            try:
                return find_count_mapping(
                    nested
                )
            except RuntimeError:
                pass

    raise RuntimeError(
        "Could not locate read-count arrays "
        "in loaded object"
    )


def validate_counts(data, label):
    missing = REQUIRED_COUNTS - set(data)

    if missing:
        raise RuntimeError(
            f"{label}: missing arrays "
            f"{sorted(missing)}"
        )

    lengths = {
        key: len(
            np.asarray(data[key])
        )
        for key in REQUIRED_COUNTS
    }

    if len(set(lengths.values())) != 1:
        raise RuntimeError(
            f"{label}: inconsistent array lengths "
            f"{lengths}"
        )

    site_count = next(
        iter(lengths.values())
    )

    if site_count == 0:
        raise RuntimeError(
            f"{label}: zero sites"
        )

    return site_count


def allele_fractions(data):
    left_alt = np.asarray(
        data["left_alt"],
        dtype=float,
    )

    left_total = np.asarray(
        data["left_total"],
        dtype=float,
    )

    right_alt = np.asarray(
        data["right_alt"],
        dtype=float,
    )

    right_total = np.asarray(
        data["right_total"],
        dtype=float,
    )

    valid = (
        np.isfinite(left_alt)
        & np.isfinite(left_total)
        & np.isfinite(right_alt)
        & np.isfinite(right_total)
        & (left_total > 0)
        & (right_total > 0)
        & (left_alt >= 0)
        & (right_alt >= 0)
        & (left_alt <= left_total)
        & (right_alt <= right_total)
    )

    left = (
        left_alt[valid]
        / left_total[valid]
    )

    right = (
        right_alt[valid]
        / right_total[valid]
    )

    return left, right, valid


def choose_median_row(
    results,
    scenario,
    generating_model,
):
    subset = results.loc[
        results["scenario"].eq(
            scenario
        )
        & results[
            "generating_model"
        ].eq(generating_model)
    ].copy()

    if subset.empty:
        raise RuntimeError(
            "No simulations for "
            f"{scenario} / {generating_model}"
        )

    subset["LR"] = pd.to_numeric(
        subset["LR"],
        errors="coerce",
    )

    subset = subset.loc[
        np.isfinite(subset["LR"])
    ].copy()

    target = float(
        subset["LR"].median()
    )

    subset["distance_from_median_LR"] = (
        subset["LR"] - target
    ).abs()

    return (
        subset.sort_values(
            [
                "distance_from_median_LR",
                "task",
            ]
        )
        .iloc[0]
        .copy()
    )


def histogram(left, right):
    values, _, _ = np.histogram2d(
        left,
        right,
        bins=BINS,
        range=[
            [0.0, 1.0],
            [0.0, 1.0],
        ],
    )

    return values


def centroid_coordinates(
    c_left,
    c_right,
):
    compatible = []
    opposing = []

    for p_left in [
        0.0,
        0.5,
        1.0,
    ]:
        for p_right in [
            0.0,
            0.5,
            1.0,
        ]:
            q_left = (
                (1.0 - c_left)
                * p_left
                + c_left
                * p_right
            )

            q_right = (
                c_right
                * p_left
                + (1.0 - c_right)
                * p_right
            )

            point = (
                q_left,
                q_right,
            )

            if (
                (p_left == 0.0 and p_right == 1.0)
                or
                (p_left == 1.0 and p_right == 0.0)
            ):
                opposing.append(point)
            else:
                compatible.append(point)

    return compatible, opposing


def add_centroids(
    axis,
    c_left,
    c_right,
):
    compatible, opposing = (
        centroid_coordinates(
            c_left,
            c_right,
        )
    )

    axis.scatter(
        [point[0] for point in compatible],
        [point[1] for point in compatible],
        marker="s",
        s=34,
        facecolor="white",
        edgecolor=BLUE,
        linewidth=1.0,
        alpha=0.80,
        zorder=5,
    )

    axis.scatter(
        [point[0] for point in opposing],
        [point[1] for point in opposing],
        marker="D",
        s=40,
        facecolor=ORANGE,
        edgecolor="white",
        linewidth=0.8,
        alpha=0.90,
        zorder=6,
    )


def configure_axis(axis, title):
    axis.set(
        xlim=(0.0, 1.0),
        ylim=(0.0, 1.0),
        aspect="equal",
        xlabel="Left ALT-read fraction",
        ylabel="Right ALT-read fraction",
        title=title,
    )

    axis.set_title(
        title,
        fontsize=10.5,
        fontweight="bold",
        loc="left",
    )

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def density_page(
    pdf,
    datasets,
    histories,
    overlay_centroids,
):
    maximum = max(
        float(np.max(histogram_value))
        for histogram_value in histories.values()
    )

    norm = LogNorm(
        vmin=1.0,
        vmax=max(
            maximum,
            1.0,
        ),
    )

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(13.2, 8.5),
    )

    axes = axes.ravel()

    image = None

    for index, dataset in enumerate(
        datasets
    ):
        axis = axes[index]

        counts = histories[
            dataset["key"]
        ]

        image = axis.imshow(
            counts.T,
            origin="lower",
            extent=[
                0.0,
                1.0,
                0.0,
                1.0,
            ],
            cmap="Greys",
            norm=norm,
            interpolation="nearest",
            aspect="equal",
        )

        if overlay_centroids:
            add_centroids(
                axis,
                dataset["c_left"],
                dataset["c_right"],
            )

        configure_axis(
            axis,
            (
                f"{dataset['title']}\n"
                f"n={dataset['n']:,}; "
                f"cL={dataset['c_left']:.3f}, "
                f"cR={dataset['c_right']:.3f}"
            ),
        )

    axes[-1].axis("off")

    if overlay_centroids:
        page_title = (
            "All eligible sites with model centroids overlaid"
        )

        subtitle = (
            "Centroids are annotations only: no sites were "
            "selected by genotype state, posterior probability "
            "or distance from a centroid."
        )
    else:
        page_title = (
            "Observed and simulated allele-fraction space: "
            "all eligible sites"
        )

        subtitle = (
            "Every panel uses the same fixed bins and shared "
            "count scale. No post-fit anchor selection is used."
        )

    figure.suptitle(
        page_title,
        fontsize=17,
        fontweight="bold",
        x=0.055,
        y=0.985,
        ha="left",
    )

    figure.text(
        0.055,
        0.940,
        subtitle,
        fontsize=10.2,
        color="#555555",
        ha="left",
    )

    colorbar = figure.colorbar(
        image,
        ax=axes[:5].tolist(),
        fraction=0.025,
        pad=0.025,
    )

    colorbar.set_label(
        "Sites per fixed AF bin",
        fontsize=9.5,
    )

    figure.text(
        0.055,
        0.025,
        (
            "Inclusion rule: fixed Stage-201 analysis sites with "
            "nonzero Left and Right allele-assigned depth. "
            "Simulation rows were selected only by median total LR; "
            "individual sites were not selected by model fit."
        ),
        fontsize=8.7,
        color="#666666",
        ha="left",
    )

    figure.subplots_adjust(
        left=0.07,
        right=0.90,
        top=0.89,
        bottom=0.09,
        wspace=0.27,
        hspace=0.34,
    )

    pdf.savefig(
        figure,
        bbox_inches="tight",
    )

    if not overlay_centroids:
        figure.savefig(
            OUTPUT_PNG,
            dpi=240,
            bbox_inches="tight",
        )

    plt.close(figure)


def residual_page(
    pdf,
    datasets,
    histories,
):
    observed = histories["observed"]

    observed_probability = (
        observed
        / observed.sum()
    )

    residuals = []

    for dataset in datasets[1:]:
        simulation = histories[
            dataset["key"]
        ]

        simulation_probability = (
            simulation
            / simulation.sum()
        )

        pseudocount = 0.5 / max(
            observed.sum(),
            simulation.sum(),
        )

        residual = np.log10(
            (
                observed_probability
                + pseudocount
            )
            /
            (
                simulation_probability
                + pseudocount
            )
        )

        residuals.append(
            (
                dataset,
                residual,
            )
        )

    absolute_values = np.concatenate(
        [
            np.abs(residual).ravel()
            for _, residual in residuals
        ]
    )

    limit = float(
        np.quantile(
            absolute_values,
            0.99,
        )
    )

    limit = max(
        limit,
        0.25,
    )

    norm = TwoSlopeNorm(
        vmin=-limit,
        vcenter=0.0,
        vmax=limit,
    )

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(10.8, 9.0),
    )

    axes = axes.ravel()

    image = None

    for axis, (
        dataset,
        residual,
    ) in zip(
        axes,
        residuals,
    ):
        image = axis.imshow(
            residual.T,
            origin="lower",
            extent=[
                0.0,
                1.0,
                0.0,
                1.0,
            ],
            cmap="RdBu_r",
            norm=norm,
            interpolation="nearest",
            aspect="equal",
        )

        configure_axis(
            axis,
            (
                "Observed minus "
                f"{dataset['title']}"
            ),
        )

    figure.suptitle(
        "Where observed AF density differs from each simulation",
        fontsize=17,
        fontweight="bold",
        x=0.055,
        y=0.985,
        ha="left",
    )

    figure.text(
        0.055,
        0.943,
        (
            "Red bins contain proportionally more observed sites; "
            "blue bins contain proportionally more simulated sites. "
            "All sites and fixed bins are retained."
        ),
        fontsize=10.2,
        color="#555555",
        ha="left",
    )

    colorbar = figure.colorbar(
        image,
        ax=axes.tolist(),
        fraction=0.025,
        pad=0.025,
    )

    colorbar.set_label(
        "log10(observed density / simulated density)",
        fontsize=9.5,
    )

    figure.text(
        0.055,
        0.025,
        (
            "Residual maps compare normalized bin occupancy, "
            "so differences reveal state frequency, dispersion "
            "and between-centroid density rather than sample size."
        ),
        fontsize=8.7,
        color="#666666",
        ha="left",
    )

    figure.subplots_adjust(
        left=0.08,
        right=0.88,
        top=0.89,
        bottom=0.09,
        wspace=0.25,
        hspace=0.32,
    )

    pdf.savefig(
        figure,
        bbox_inches="tight",
    )

    plt.close(figure)


def main():
    empirical = load_module(
        EMPIRICAL_SCRIPT,
        "stage202_for_unselected_AF",
    )

    if not hasattr(empirical, "stage201"):
        raise RuntimeError(
            "Stage 202 did not expose its Stage-201 module"
        )

    loaded_input = (
        empirical.stage201.load_input()
    )

    observed_data = find_count_mapping(
        loaded_input
    )

    observed_total_sites = validate_counts(
        observed_data,
        "observed",
    )

    results = pd.read_csv(
        RESULTS_TABLE,
        sep="\t",
        compression="gzip",
    )

    required_columns = {
        "scenario",
        "generating_model",
        "LR",
        "task",
        "true_c_left",
        "true_c_right",
    }

    missing = (
        required_columns
        - set(results.columns)
    )

    if missing:
        raise RuntimeError(
            "Results table missing columns: "
            f"{sorted(missing)}"
        )

    selected_specifications = [
        (
            "exact_20",
            "three_founder",
            "Three-founder, exact 20%",
        ),
        (
            "exact_20",
            "empirical_full_sib",
            "Full-sib, exact 20%",
        ),
        (
            "observed_matched",
            "three_founder",
            "Three-founder, fitted mixing",
        ),
        (
            "observed_matched",
            "empirical_full_sib",
            "Full-sib, fitted mixing",
        ),
    ]

    selected_rows = []

    for (
        scenario,
        generating_model,
        title,
    ) in selected_specifications:
        row = choose_median_row(
            results,
            scenario,
            generating_model,
        )

        row["display_title"] = title

        selected_rows.append(row)

    selected_table = pd.DataFrame(
        selected_rows
    )

    observed_fit = json.loads(
        OBSERVED_FIT.read_text()
    )

    observed_parameters = (
        observed_fit["alternative"]
    )

    observed_left, observed_right, _ = (
        allele_fractions(
            observed_data
        )
    )

    datasets = [
        {
            "key": "observed",
            "title": "Observed gynandromorph",
            "data": observed_data,
            "left": observed_left,
            "right": observed_right,
            "n": len(observed_left),
            "c_left": float(
                observed_parameters["c_left"]
            ),
            "c_right": float(
                observed_parameters["c_right"]
            ),
            "source": (
                "fixed Stage-201 analysis input"
            ),
        }
    ]

    for row in selected_rows:
        reconstructed = (
            empirical.reconstruct(
                row,
                observed_data,
            )
        )

        reconstructed = find_count_mapping(
            reconstructed
        )

        reconstructed_sites = validate_counts(
            reconstructed,
            row["display_title"],
        )

        if reconstructed_sites != observed_total_sites:
            raise RuntimeError(
                f"{row['display_title']}: "
                f"{reconstructed_sites} sites versus "
                f"{observed_total_sites} observed sites"
            )

        left, right, _ = allele_fractions(
            reconstructed
        )

        key = (
            f"{row['scenario']}_"
            f"{row['generating_model']}"
        )

        datasets.append(
            {
                "key": key,
                "title": row[
                    "display_title"
                ],
                "data": reconstructed,
                "left": left,
                "right": right,
                "n": len(left),
                "c_left": float(
                    row["true_c_left"]
                ),
                "c_right": float(
                    row["true_c_right"]
                ),
                "source": (
                    "median-LR replicate reconstructed "
                    "without site selection"
                ),
            }
        )

    histories = {
        dataset["key"]: histogram(
            dataset["left"],
            dataset["right"],
        )
        for dataset in datasets
    }

    OUTPUT_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    selected_table.to_csv(
        OUTPUT_ROWS,
        sep="\t",
        index=False,
    )

    with PdfPages(
        OUTPUT_PDF
    ) as pdf:
        density_page(
            pdf,
            datasets,
            histories,
            overlay_centroids=False,
        )

        density_page(
            pdf,
            datasets,
            histories,
            overlay_centroids=True,
        )

        residual_page(
            pdf,
            datasets,
            histories,
        )

    audit = {
        "observed_input_sites": int(
            observed_total_sites
        ),
        "AF_bins_per_axis": BINS,
        "site_inclusion": (
            "all fixed Stage-201 analysis sites with "
            "nonzero Left and Right allele-assigned depth"
        ),
        "postfit_anchor_selection_used": False,
        "posterior_filter_used": False,
        "centroid_distance_filter_used": False,
        "LR_contribution_filter_used": False,
        "simulation_replicate_selection": (
            "one median-total-LR replicate per "
            "scenario and generating model"
        ),
        "shared_density_bins": True,
        "shared_density_colour_normalization": True,
        "pages": [
            "all-site raw AF densities",
            "same all-site densities with centroid annotations",
            "observed-versus-simulated normalized density residuals",
        ],
        "datasets": [
            {
                "key": dataset["key"],
                "title": dataset["title"],
                "plotted_sites": int(
                    dataset["n"]
                ),
                "c_left": dataset["c_left"],
                "c_right": dataset["c_right"],
                "source": dataset["source"],
            }
            for dataset in datasets
        ],
        "pdf": str(OUTPUT_PDF),
        "page1_png": str(OUTPUT_PNG),
        "selected_rows": str(
            OUTPUT_ROWS
        ),
    }

    OUTPUT_AUDIT.write_text(
        json.dumps(
            audit,
            indent=2,
        )
        + "\n"
    )

    print("STAGE_216_COMPLETE")

    print(
        "PDF\t"
        f"{OUTPUT_PDF.resolve()}"
    )

    print(
        "PAGE1_PNG\t"
        f"{OUTPUT_PNG.resolve()}"
    )

    print(
        "AUDIT\t"
        f"{OUTPUT_AUDIT.resolve()}"
    )


if __name__ == "__main__":
    main()
