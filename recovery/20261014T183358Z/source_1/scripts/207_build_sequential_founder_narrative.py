#!/usr/bin/env python3
import importlib.util
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import numpy as np
import pandas as pd


ROOT = Path(
    "analysis/genomewide_mixing_corrected_v1"
)
EMP = ROOT / "empirical_full_sib_calibration_v1"
OUT = ROOT / "sequential_founder_narrative_v1"

PDF = (
    OUT
    / "figures"
    / "sequential_founder_mixing_narrative.pdf"
)

AUDIT = (
    OUT
    / "sequential_founder_mixing_narrative.audit.json"
)

SELECTED = (
    OUT
    / "tables"
    / "displayed_sibling_scenarios.tsv"
)

BLUE = "#0072B2"
ORANGE = "#D55E00"

GREY = LinearSegmentedColormap.from_list(
    "grey_density",
    [
        "#FFFFFF",
        "#E5E5E5",
        "#A8A8A8",
        "#555555",
        "#181818",
    ],
)

COUNT_KEYS = {
    "left_alt",
    "left_total",
    "right_alt",
    "right_total",
}


def module(path, name):
    specification = importlib.util.spec_from_file_location(
        name,
        path,
    )

    value = importlib.util.module_from_spec(
        specification
    )

    specification.loader.exec_module(value)

    return value


def find_table():
    required_columns = {
        "generating_model",
        "true_c_left",
        "true_c_right",
        "LR",
    }

    candidates = []

    for root in (
        ROOT / "expanded_calibration_v2",
        ROOT / "expanded_calibration",
    ):
        if not root.is_dir():
            continue

        for path in root.rglob("*.tsv*"):
            try:
                columns = set(
                    pd.read_csv(
                        path,
                        sep="\t",
                        nrows=0,
                    ).columns
                )
            except Exception:
                continue

            if required_columns <= columns:
                candidates.append(path)

    if not candidates:
        raise RuntimeError(
            "Expanded unrelated-donor simulation "
            "table not found"
        )

    return max(
        candidates,
        key=lambda path: path.stat().st_size,
    )


def model_name(value):
    value = (
        str(value)
        .lower()
        .replace("-", "_")
    )

    explicit_models = {
        "three": "three",
        "three_founder": "three",
        "four": "four",
        "four_founder": "four",
        "empirical_full_sib": "four",
    }

    if value in explicit_models:
        return explicit_models[value]

    if "three" in value:
        return "three"

    if "four" in value:
        return "four"

    if (
        "full_sib" in value
        or "sibling" in value
    ):
        return "four"

    return value


def count_sets(value, label="root"):
    answer = []

    if isinstance(value, dict):
        if COUNT_KEYS <= set(value):
            answer.append(
                (
                    label,
                    value,
                )
            )
        else:
            for key, child in value.items():
                answer += count_sets(
                    child,
                    f"{label}.{key}",
                )

    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            answer += count_sets(
                child,
                f"{label}.{index}",
            )

    return answer


def one_count_set(value, label):
    found = count_sets(value)

    if len(found) != 1:
        raise RuntimeError(
            f"{label}: expected one count dataset; "
            f"found {len(found)}"
        )

    return found[0][1]


def benchmark_pair(
    plotter,
    simulations,
    label,
):
    data = plotter.stage.load_data(
        calibration=False
    )

    matrix = data["controls"]
    called = matrix >= 0

    allele_frequency = (
        np.where(
            called,
            matrix,
            0,
        ).sum(axis=1)
        + 0.5
    ) / (
        2 * called.sum(axis=1)
        + 1
    )

    output = {}

    for founder_index, founder in enumerate(
        (
            "three",
            "four",
        )
    ):
        subset = simulations.loc[
            simulations[
                "generating_model"
            ].map(model_name)
            == founder
        ].copy()

        if subset.empty:
            raise RuntimeError(
                f"{label}: no "
                f"{founder}-founder rows"
            )

        # Prefer the high-replication exact-20%
        # family whenever it is available.
        if (
            "family" in subset.columns
            and (
                subset["family"]
                == "repeat_20"
            ).any()
        ):
            subset = subset.loc[
                subset["family"]
                == "repeat_20"
            ].copy()

        target_left = float(
            subset["true_c_left"].median()
        )

        target_right = float(
            subset["true_c_right"].median()
        )

        subset["recovery_error"] = (
            (
                subset[
                    "estimated_c_left_alt"
                ]
                - target_left
            ).abs()
            + (
                subset[
                    "estimated_c_right_alt"
                ]
                - target_right
            ).abs()
        )

        subset = subset.sort_values(
            [
                "recovery_error",
                "id",
            ]
        )

        row = subset.iloc[
            len(subset) // 2
        ]

        identifier = int(
            row["id"]
        )

        donor_a_index = (
            plotter.stage.CONTROLS.index(
                row["donor_a"]
            )
        )

        donor_b_index = (
            plotter.stage.CONTROLS.index(
                row["donor_b"]
            )
        )

        rng = np.random.default_rng(
            207000000
            + identifier * 2
            + founder_index
        )

        donor_a = plotter.stage.impute(
            matrix,
            donor_a_index,
            allele_frequency,
            rng,
        )

        donor_b = plotter.stage.impute(
            matrix,
            donor_b_index,
            allele_frequency,
            rng,
        )

        if founder == "three":
            constructed = (
                plotter.stage.gamete(
                    donor_a,
                    rng,
                )
                + plotter.stage.gamete(
                    donor_b,
                    rng,
                )
            )

            if identifier % 2 == 0:
                true_left = donor_a
                true_right = constructed
            else:
                true_left = constructed
                true_right = donor_a
        else:
            true_left = donor_a
            true_right = donor_b

        output[founder] = (
            plotter.pilot.simulate_counts(
                true_left,
                true_right,
                target_left,
                target_right,
                data,
                plotter.pilot.load_emission(),
                rng,
            )
        )

    return output


def fractions(data):
    left_alt = np.asarray(
        data["left_alt"],
        dtype=float,
    )
    left_depth = np.asarray(
        data["left_total"],
        dtype=float,
    )
    right_alt = np.asarray(
        data["right_alt"],
        dtype=float,
    )
    right_depth = np.asarray(
        data["right_total"],
        dtype=float,
    )

    keep = (
        np.isfinite(
            left_alt
            + left_depth
            + right_alt
            + right_depth
        )
        & (left_depth > 0)
        & (right_depth > 0)
    )

    return (
        left_alt[keep] / left_depth[keep],
        right_alt[keep] / right_depth[keep],
    )


def centroids(
    axis,
    c_left,
    c_right,
    four,
):
    normal = []
    opposing = []

    for left in (
        0.0,
        0.5,
        1.0,
    ):
        for right in (
            0.0,
            0.5,
            1.0,
        ):
            point = (
                (
                    (1.0 - c_left) * left
                    + c_left * right
                ),
                (
                    (1.0 - c_right) * right
                    + c_right * left
                ),
            )

            if (
                left,
                right,
            ) in (
                (0.0, 1.0),
                (1.0, 0.0),
            ):
                opposing.append(point)
            else:
                normal.append(point)

    axis.scatter(
        *zip(*normal),
        s=28,
        facecolors="white",
        edgecolors=BLUE,
        linewidths=1.1,
        zorder=4,
    )

    if four:
        axis.scatter(
            *zip(*opposing),
            s=34,
            marker="D",
            color=ORANGE,
            edgecolors="white",
            linewidths=0.5,
            zorder=5,
        )


def max_hex(data):
    x_values, y_values = fractions(data)

    figure, axis = plt.subplots()

    values = axis.hexbin(
        x_values,
        y_values,
        gridsize=62,
        extent=(
            0.0,
            1.0,
            0.0,
            1.0,
        ),
        mincnt=1,
    ).get_array()

    plt.close(figure)

    return max(
        2.0,
        float(values.max()),
    )


def draw_hex(
    axis,
    data,
    norm,
    title,
    c_left,
    c_right,
    four,
):
    x_values, y_values = fractions(data)

    artist = axis.hexbin(
        x_values,
        y_values,
        gridsize=62,
        extent=(
            0.0,
            1.0,
            0.0,
            1.0,
        ),
        mincnt=1,
        cmap=GREY,
        norm=norm,
        linewidths=0,
    )

    axis.plot(
        [0.0, 1.0],
        [0.0, 1.0],
        color="#D0D0D0",
        linewidth=0.7,
        zorder=0,
    )

    centroids(
        axis,
        c_left,
        c_right,
        four,
    )

    axis.set(
        xlim=(-0.02, 1.02),
        ylim=(-0.02, 1.02),
        title=title,
    )

    axis.set_aspect("equal")
    axis.tick_params(labelsize=8)

    axis.text(
        0.02,
        0.98,
        f"n = {x_values.size:,}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=7,
        color="#666666",
    )

    return artist


def header(
    figure,
    step,
    title,
    subtitle,
):
    figure.text(
        0.055,
        0.965,
        f"STEP {step}",
        va="top",
        fontsize=9,
        fontweight="bold",
        color="#666666",
    )

    figure.text(
        0.055,
        0.925,
        title,
        va="top",
        fontsize=19,
        fontweight="bold",
        color="#222222",
    )

    figure.text(
        0.055,
        0.875,
        subtitle,
        va="top",
        fontsize=9.5,
        color="#666666",
    )


def conceptual_page(pdf):
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(11.0, 7.5),
    )

    figure.subplots_adjust(
        left=0.08,
        right=0.95,
        bottom=0.14,
        top=0.76,
        wspace=0.25,
    )

    header(
        figure,
        1,
        "Define the biological alternatives",
        (
            "Founder sharing and sequencing "
            "mixture are fitted separately."
        ),
    )

    for axis, four in zip(
        axes,
        (
            False,
            True,
        ),
    ):
        color = (
            ORANGE
            if four
            else BLUE
        )

        axis.axis("off")

        axis.set_title(
            (
                "Four-founder hypothesis"
                if four
                else "Three-founder hypothesis"
            ),
            color=color,
            fontsize=14,
            fontweight="bold",
        )

        axis.text(
            0.5,
            0.68,
            "Left = A + B",
            ha="center",
            fontsize=13,
            bbox={
                "boxstyle": "round",
                "fc": "white",
                "ec": color,
            },
        )

        axis.text(
            0.5,
            0.38,
            (
                "Right = C + D"
                if four
                else "Right = A + C"
            ),
            ha="center",
            fontsize=13,
            bbox={
                "boxstyle": "round",
                "fc": "white",
                "ec": color,
            },
        )

        axis.text(
            0.5,
            0.15,
            (
                "All nine genotype pairs possible"
                if four
                else (
                    "Shared gamete A; opposing "
                    "homozygotes impossible"
                )
            ),
            ha="center",
            color=color,
            fontweight="bold",
            fontsize=10,
        )

    figure.text(
        0.5,
        0.06,
        (
            "These are three-founder and "
            "four-founder models, not triploid "
            "and tetraploid states."
        ),
        ha="center",
        fontsize=10,
        color="#555555",
    )

    pdf.savefig(figure)
    plt.close(figure)


def factorial_page(
    pdf,
    pairs,
    step,
    title,
    subtitle,
    model_titles=None,
):
    if model_titles is None:
        model_titles = (
            "Three-founder",
            "Four-founder",
        )
    datasets = [
        pairs[0]["three"],
        pairs[0]["four"],
        pairs[1]["three"],
        pairs[1]["four"],
    ]

    norm = LogNorm(
        vmin=1,
        vmax=max(
            max_hex(value)
            for value in datasets
        ),
    )

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(10.6, 9.0),
        sharex=True,
        sharey=True,
    )

    figure.subplots_adjust(
        left=0.10,
        right=0.89,
        bottom=0.09,
        top=0.78,
        hspace=0.24,
        wspace=0.15,
    )

    header(
        figure,
        step,
        title,
        subtitle,
    )

    artist = None

    for row, mixing in enumerate(
        (
            0.0,
            0.2,
        )
    ):
        for column, model in enumerate(
            (
                "three",
                "four",
            )
        ):
            artist = draw_hex(
                axes[row, column],
                pairs[row][model],
                norm,
                model_titles[column],
                mixing,
                mixing,
                model == "four",
            )

    for row, text in enumerate(
        (
            "No mixing",
            "20% mixing",
        )
    ):
        axes[row, 0].text(
            -0.27,
            0.5,
            text,
            rotation=90,
            transform=axes[row, 0].transAxes,
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
        )

    for axis in axes[1]:
        axis.set_xlabel(
            "Left ALT-read fraction"
        )

    for axis in axes[:, 0]:
        axis.set_ylabel(
            "Right ALT-read fraction"
        )

    color_axis = figure.add_axes(
        [
            0.915,
            0.17,
            0.018,
            0.50,
        ]
    )

    figure.colorbar(
        artist,
        cax=color_axis,
    ).set_label(
        "Sites per hexagon"
    )

    pdf.savefig(figure)
    plt.close(figure)


def nested_arrays(value):
    answer = []

    if isinstance(value, dict):
        for child in value.values():
            answer += nested_arrays(child)

    elif isinstance(value, (list, tuple)):
        for child in value:
            answer += nested_arrays(child)

    elif isinstance(value, np.ndarray):
        answer.append(value)

    return answer


def add_heterozygosity(
    results,
    sibling_input,
    samples,
):
    candidates = []

    for array in nested_arrays(
        sibling_input
    ):
        if (
            array.ndim != 2
            or len(samples) not in array.shape
        ):
            continue

        if array.shape[1] == len(samples):
            matrix = array
        else:
            matrix = array.T

        probe = np.asarray(
            matrix
        ).ravel()[:100000]

        try:
            probe = probe.astype(float)
        except (
            TypeError,
            ValueError,
        ):
            continue

        finite = probe[
            np.isfinite(probe)
        ]

        if (
            finite.size
            and set(
                np.rint(finite).astype(int)
            ) <= {
                -1,
                0,
                1,
                2,
            }
        ):
            candidates.append(matrix)

    if not candidates:
        raise RuntimeError(
            "Could not identify the sibling "
            "genotype matrix"
        )

    matrix = max(
        candidates,
        key=lambda value: value.size,
    ).astype(float)

    rates = (
        (matrix == 1).sum(axis=0)
        / np.maximum(
            1,
            (matrix >= 0).sum(axis=0),
        )
    )

    rates = dict(
        zip(
            samples,
            rates,
        )
    )

    output = results.copy()

    output["left_H"] = (
        output["donor_left"]
        .map(rates)
    )

    output["right_H"] = (
        output["donor_right"]
        .map(rates)
    )

    if (
        output[
            [
                "left_H",
                "right_H",
            ]
        ]
        .isna()
        .any()
        .any()
    ):
        raise RuntimeError(
            "Sibling names did not match the "
            "genotype matrix"
        )

    output["pair_mean_H"] = (
        output["left_H"]
        + output["right_H"]
    ) / 2.0

    return output


def choose_siblings(results):
    exact = results[
        np.isclose(
            results["true_c_left"],
            0.2,
        )
        & np.isclose(
            results["true_c_right"],
            0.2,
        )
    ].copy()

    exact["model"] = (
        exact["generating_model"]
        .map(model_name)
    )

    chosen = []

    for model in (
        "three",
        "four",
    ):
        group = exact[
            exact["model"] == model
        ]

        if group.empty:
            raise RuntimeError(
                "No 20% sibling simulations "
                f"for {model}"
            )

        for label, quantile in (
            (
                "low",
                0.1,
            ),
            (
                "middle",
                0.5,
            ),
            (
                "high",
                0.9,
            ),
        ):
            target = (
                group["pair_mean_H"]
                .quantile(quantile)
            )

            index = (
                group["pair_mean_H"]
                .sub(target)
                .abs()
                .idxmin()
            )

            row = group.loc[index].copy()
            row["H_stratum"] = label
            chosen.append(row)

    return pd.DataFrame(chosen)


def reconstruct(
    finisher,
    sibling_input,
    row,
    mixing,
):
    row = row.copy()

    row["true_c_left"] = mixing
    row["true_c_right"] = mixing

    return one_count_set(
        finisher.reconstruct(
            row,
            sibling_input,
        ),
        f"task {row['task']}",
    )


def sibling_pages(
    pdf,
    finisher,
    sibling_input,
    selected,
):
    middle = {
        model: selected[
            (
                selected["model"] == model
            )
            & (
                selected["H_stratum"]
                == "middle"
            )
        ].iloc[0]
        for model in (
            "three",
            "four",
        )
    }

    pairs = []

    for mixing in (
        0.0,
        0.2,
    ):
        pairs.append(
            {
                model: reconstruct(
                    finisher,
                    sibling_input,
                    row,
                    mixing,
                )
                for model, row in middle.items()
            }
        )

    factorial_page(
        pdf,
        pairs,
        4,
        "Add empirical full-sibling relatedness",
        (
            "Real sibling genotypes at zero "
            "and 20% mixing; no likelihood "
            "model is changed."
        ),
        model_titles=(
            "Three-founder",
            "Empirical full-sib",
        ),
    )

    built = []

    for _, row in selected.iterrows():
        built.append(
            (
                row,
                reconstruct(
                    finisher,
                    sibling_input,
                    row,
                    0.2,
                ),
            )
        )

    norm = LogNorm(
        vmin=1,
        vmax=max(
            max_hex(data)
            for _, data in built
        ),
    )

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(12.3, 7.6),
        sharex=True,
        sharey=True,
    )

    figure.subplots_adjust(
        left=0.08,
        right=0.91,
        bottom=0.10,
        top=0.76,
        hspace=0.24,
        wspace=0.12,
    )

    header(
        figure,
        5,
        (
            "Vary heterozygosity within the "
            "sibling pedigree"
        ),
        (
            "Low, middle and high mean donor "
            "heterozygosity; genotype contrast "
            "and IBD0 are shown."
        ),
    )

    artist = None

    for row, data in built:
        row_index = (
            0
            if row["model"] == "three"
            else 1
        )

        column_index = {
            "low": 0,
            "middle": 1,
            "high": 2,
        }[row["H_stratum"]]

        artist = draw_hex(
            axes[
                row_index,
                column_index,
            ],
            data,
            norm,
            row["H_stratum"].capitalize(),
            0.2,
            0.2,
            row["model"] == "four",
        )

        axes[
            row_index,
            column_index,
        ].text(
            0.98,
            0.02,
            (
                f"mean H="
                f"{row['pair_mean_H']:.3f}\n"
                f"contrast="
                f"{row['genotype_contrast']:.3f}\n"
                f"IBD0="
                f"{row['IBD0_fraction']:.2f}"
            ),
            transform=axes[
                row_index,
                column_index,
            ].transAxes,
            ha="right",
            va="bottom",
            fontsize=7,
            color="#555555",
        )

    axes[0, 0].text(
        -0.27,
        0.5,
        "Three-founder",
        rotation=90,
        transform=axes[0, 0].transAxes,
        ha="center",
        va="center",
        fontweight="bold",
        color=BLUE,
    )

    axes[1, 0].text(
        -0.27,
        0.5,
        "Empirical full-sib",
        rotation=90,
        transform=axes[1, 0].transAxes,
        ha="center",
        va="center",
        fontweight="bold",
        color=ORANGE,
    )

    for axis in axes[1]:
        axis.set_xlabel(
            "Left ALT-read fraction"
        )

    for axis in axes[:, 0]:
        axis.set_ylabel(
            "Right ALT-read fraction"
        )

    color_axis = figure.add_axes(
        [
            0.93,
            0.16,
            0.016,
            0.48,
        ]
    )

    figure.colorbar(
        artist,
        cax=color_axis,
    ).set_label(
        "Sites per hexagon"
    )

    pdf.savefig(figure)
    plt.close(figure)

    return pairs[1]


def observed_page(
    pdf,
    observed,
    sibling_mixed,
    c_left,
    c_right,
):
    datasets = [
        observed,
        sibling_mixed["three"],
        sibling_mixed["four"],
    ]

    norm = LogNorm(
        vmin=1,
        vmax=max(
            max_hex(data)
            for data in datasets
        ),
    )

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(12.5, 5.3),
        sharex=True,
        sharey=True,
    )

    figure.subplots_adjust(
        left=0.075,
        right=0.91,
        bottom=0.14,
        top=0.72,
        wspace=0.14,
    )

    header(
        figure,
        6,
        (
            "Compare the animal with matched "
            "sibling benchmarks"
        ),
        (
            "Observed and simulated panels use "
            "the same genomic scope and density "
            "scale."
        ),
    )

    settings = [
        (
            observed,
            "Observed Left / Right",
            c_left,
            c_right,
            True,
        ),
        (
            sibling_mixed["three"],
            "Sibling three-founder, 20%",
            0.2,
            0.2,
            False,
        ),
        (
            sibling_mixed["four"],
            "Empirical full-sib, 20%",
            0.2,
            0.2,
            True,
        ),
    ]

    artist = None

    for axis, setting in zip(
        axes,
        settings,
    ):
        (
            data,
            title,
            left_mix,
            right_mix,
            four,
        ) = setting

        artist = draw_hex(
            axis,
            data,
            norm,
            title,
            left_mix,
            right_mix,
            four,
        )

        axis.set_xlabel(
            "Left ALT-read fraction"
        )

    axes[0].set_ylabel(
        "Right ALT-read fraction"
    )

    color_axis = figure.add_axes(
        [
            0.93,
            0.19,
            0.016,
            0.40,
        ]
    )

    figure.colorbar(
        artist,
        cax=color_axis,
    ).set_label(
        "Sites per hexagon"
    )

    pdf.savefig(figure)
    plt.close(figure)


def evidence_page(
    pdf,
    results,
    summary,
    observed_lr,
):
    exact = results[
        np.isclose(
            results["true_c_left"],
            0.2,
        )
        & np.isclose(
            results["true_c_right"],
            0.2,
        )
    ].copy()

    exact["model"] = (
        exact["generating_model"]
        .map(model_name)
    )

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(12.0, 5.6),
    )

    figure.subplots_adjust(
        left=0.075,
        right=0.97,
        bottom=0.16,
        top=0.72,
        wspace=0.28,
    )

    header(
        figure,
        7,
        (
            "Calibrate founder evidence with "
            "empirical siblings"
        ),
        (
            "The observed LR is compared with "
            "simulations from the real "
            "full-sibling pedigree."
        ),
    )

    for model, color, label in (
        (
            "three",
            BLUE,
            "Three-founder",
        ),
        (
            "four",
            ORANGE,
            "Empirical full-sib",
        ),
    ):
        values = exact.loc[
            exact["model"] == model,
            "LR",
        ].to_numpy(
            dtype=float
        )

        axes[0].hist(
            np.log10(
                1.0 + values
            ),
            bins=35,
            histtype="step",
            linewidth=2,
            density=True,
            color=color,
            label=(
                f"{label} "
                f"(n={len(values)})"
            ),
        )

    axes[0].axvline(
        np.log10(
            1.0 + observed_lr
        ),
        color="black",
        linewidth=1.8,
        label="Observed",
    )

    axes[0].set(
        xlabel="log10(1 + LR)",
        ylabel="Simulation density",
    )

    axes[0].legend(
        frameon=False,
        fontsize=8,
    )

    table = summary.copy()

    table["model"] = (
        table["generating_model"]
        .map(model_name)
    )

    labels = (
        table["scenario"].astype(str)
        + "\n"
        + table["model"]
    )

    colors = table["model"].map(
        {
            "three": BLUE,
            "four": ORANGE,
        }
    )

    axes[1].scatter(
        range(len(table)),
        table[
            "fraction_simulated_LR_at_or_above_observed"
        ],
        c=colors,
        s=65,
    )

    axes[1].set_xticks(
        range(len(table)),
        labels,
        rotation=25,
        ha="right",
    )

    axes[1].set(
        ylabel=(
            "Fraction of simulated LR "
            "at or above observed"
        ),
        ylim=(-0.03, 1.03),
    )

    for axis in axes:
        axis.spines[
            [
                "top",
                "right",
            ]
        ].set_visible(False)

        axis.grid(
            axis="y",
            color="#E8E8E8",
            linewidth=0.7,
        )

    pdf.savefig(figure)
    plt.close(figure)


def main():
    PDF.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    SELECTED.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    pilot = module(
        "scripts/185_restored_model_pilot.py",
        "pilot207",
    )

    plotter = module(
        (
            "scripts/"
            "188_plot_corrected_genomewide_narrative.py"
        ),
        "plot207",
    )

    sibling = module(
        (
            "scripts/"
            "201_simulate_empirical_full_sib_mixing.py"
        ),
        "sibling207",
    )

    finisher = module(
        (
            "scripts/"
            "202_204_finish_empirical_sibling_calibration.py"
        ),
        "finisher207",
    )

    expanded_path = find_table()

    expanded = pd.read_csv(
        expanded_path,
        sep="\t",
    )

    sibling_results = pd.read_csv(
        (
            EMP
            / "final"
            / "all_empirical_sibling_simulations.tsv.gz"
        ),
        sep="\t",
    )

    summary = pd.read_csv(
        (
            EMP
            / "final"
            / "empirical_sibling_simulation_summary.tsv"
        ),
        sep="\t",
    )

    zero = expanded[
        np.isclose(
            expanded["true_c_left"],
            0.0,
        )
        & np.isclose(
            expanded["true_c_right"],
            0.0,
        )
    ]

    twenty = expanded[
        np.isclose(
            expanded["true_c_left"],
            0.2,
        )
        & np.isclose(
            expanded["true_c_right"],
            0.2,
        )
    ]

    if zero.empty or twenty.empty:
        raise RuntimeError(
            "Both zero- and 20%-mixing "
            "unrelated benchmarks are required"
        )

    unrelated = [
        benchmark_pair(
            plotter,
            zero,
            "zero mixing",
        ),
        benchmark_pair(
            plotter,
            twenty,
            "20% mixing",
        ),
    ]

    sibling_input = sibling.load_input()

    sibling_samples = [
        str(value)
        for value
        in sibling.offspring_samples()
    ]

    sibling_results = add_heterozygosity(
        sibling_results,
        sibling_input,
        sibling_samples,
    )

    selected = choose_siblings(
        sibling_results
    )

    selected.to_csv(
        SELECTED,
        sep="\t",
        index=False,
    )

    observed = (
        plotter.stage.load_data(
            calibration=False
        )
    )

    observed_json = (
        plotter.BASE
        / "observed_expanded.json"
    )

    with open(observed_json) as handle:
        observed_fit = json.load(handle)

    observed_full = observed_fit["full"]
    observed_parameters = (
        observed_full["alternative"]
    )

    if "LR" in observed_full:
        observed_lr = float(
            observed_full["LR"]
        )
    else:
        observed_lr = 2.0 * (
            float(
                observed_full[
                    "alternative"
                ][
                    "logL"
                ]
            )
            - float(
                observed_full[
                    "null"
                ][
                    "logL"
                ]
            )
        )

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titleweight": "bold",
            "axes.edgecolor": "#777777",
            "axes.linewidth": 0.7,
        }
    )

    with PdfPages(PDF) as pdf:
        conceptual_page(pdf)

        factorial_page(
            pdf,
            unrelated,
            2,
            "Add reciprocal sequencing mixture",
            (
                "Unrelated-donor simulations "
                "with identical site and "
                "density scales."
            ),
        )

        plotter.plot_recovery(
            pdf,
            expanded,
        )

        sibling_mixed = sibling_pages(
            pdf,
            finisher,
            sibling_input,
            selected,
        )

        observed_page(
            pdf,
            observed,
            sibling_mixed,
            float(
                observed_parameters[
                    "c_left"
                ]
            ),
            float(
                observed_parameters[
                    "c_right"
                ]
            ),
        )

        evidence_page(
            pdf,
            sibling_results,
            summary,
            observed_lr,
        )

        plotter.plot_manhattan(pdf)
        plotter.plot_z_diagnostic(pdf)

    if (
        not PDF.is_file()
        or PDF.stat().st_size == 0
    ):
        raise RuntimeError(
            "Narrative PDF was not created"
        )

    audit = {
        "status": "PASS",
        "pages_expected": 9,
        "pdf": str(PDF),
        "pdf_bytes": PDF.stat().st_size,
        "expanded_table": str(
            expanded_path
        ),
        "sibling_simulation_rows": len(
            sibling_results
        ),
        "displayed_sibling_scenarios": len(
            selected
        ),
        "zero_mixing_sibling_AF_reconstructed": True,
        "zero_mixing_sibling_model_refitted": False,
        "observed_scope": "expanded_one_million_sites",
        "observed_fit_source": str(observed_json),
        "observed_LR": observed_lr,
        "empirical_full_sib_is_strict_four_founder": False,
        "unrelated_zero_mixing_family": "global",
        "unrelated_twenty_percent_family": "repeat_20",
        "likelihood_model_changed": False,
        "focal_GT_used": False,
        "density_scale": (
            "shared absolute log count "
            "within each comparison page"
        ),
    }

    with open(AUDIT, "w") as handle:
        json.dump(
            audit,
            handle,
            indent=2,
        )
        handle.write("\n")

    print(
        "SEQUENTIAL FOUNDER NARRATIVE COMPLETE"
    )
    print(
        f"pdf\t{PDF.resolve()}"
    )
    print(
        f"audit\t{AUDIT.resolve()}"
    )
    print(
        "likelihood_model_changed\tFalse"
    )


if __name__ == "__main__":
    main()
