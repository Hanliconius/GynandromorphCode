#!/usr/bin/env python3

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(
    "analysis/genomewide_mixing_corrected_v1"
)

SIMULATION_TABLE = (
    ROOT
    / "empirical_full_sib_calibration_v1"
    / "final"
    / "all_empirical_sibling_simulations.tsv.gz"
)

OBSERVED_AUDIT = (
    ROOT
    / "expanded_calibration_v2"
    / "autosomal_1Mb_evidence.audit.json"
)

OUTPUT_DIRECTORY = (
    ROOT
    / "sequential_founder_narrative"
    / "compatibility_panel"
)

OUTPUT_PDF = (
    OUTPUT_DIRECTORY
    / "sequential_model_compatibility.pdf"
)

OUTPUT_PNG = (
    OUTPUT_DIRECTORY
    / "sequential_model_compatibility.png"
)

OUTPUT_TABLE = (
    OUTPUT_DIRECTORY
    / "sequential_model_compatibility.tsv"
)

OUTPUT_AUDIT = (
    OUTPUT_DIRECTORY
    / "sequential_model_compatibility.audit.json"
)


# Three genuinely sequential stages.
MAIN_STAGES = [
    {
        "order": 1,
        "scenario": "exact_20",
        "generating_model": "three_founder",
        "label": (
            "Three-founder baseline\n"
            "20% symmetric mixing"
        ),
        "addition": (
            "Baseline founder structure"
        ),
        "colour": "#2C9BCB",
    },
    {
        "order": 2,
        "scenario": "exact_20",
        "generating_model": "empirical_full_sib",
        "label": (
            "+ empirical full-sibling structure\n"
            "20% symmetric mixing"
        ),
        "addition": (
            "Add observed sibling genotype structure"
        ),
        "colour": "#F28E2B",
    },
    {
        "order": 3,
        "scenario": "observed_matched",
        "generating_model": "empirical_full_sib",
        "label": (
            "+ fitted asymmetric mixing\n"
            "19.7% Left / 20.24% Right"
        ),
        "addition": (
            "Add fitted side-specific mixing"
        ),
        "colour": "#B61F72",
    },
]


# This is scientifically useful but is not part of the causal
# three-stage sequence. It is retained in the output table.
SENSITIVITY_CONTROL = {
    "order": 0,
    "scenario": "observed_matched",
    "generating_model": "three_founder",
    "label": (
        "Three-founder mixing-only control"
    ),
    "addition": (
        "Fitted mixing without sibling structure"
    ),
    "colour": "#777777",
}


def require_file(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            f"Missing or empty input: {path}"
        )


def percentile_text(q):
    percentile = 100.0 * q

    if percentile < 0.1:
        return "<0.1st percentile"

    if percentile > 99.9:
        return ">99.9th percentile"

    if percentile < 1.0:
        return f"{percentile:.1f}th percentile"

    if percentile > 99.0:
        return f"{percentile:.1f}th percentile"

    return f"{percentile:.0f}th percentile"


def compatibility_row(
    simulations,
    observed_lr,
    specification,
    role,
):
    subset = simulations.loc[
        simulations["scenario"].eq(
            specification["scenario"]
        )
        & simulations["generating_model"].eq(
            specification["generating_model"]
        )
    ].copy()

    values = pd.to_numeric(
        subset["LR"],
        errors="coerce",
    ).to_numpy(dtype=float)

    values = values[
        np.isfinite(values)
    ]

    if values.size == 0:
        raise RuntimeError(
            "No finite simulations for "
            f"{specification['scenario']} / "
            f"{specification['generating_model']}"
        )

    less = int(
        np.sum(values < observed_lr)
    )

    equal = int(
        np.sum(
            np.isclose(
                values,
                observed_lr,
                rtol=1e-12,
                atol=1e-9,
            )
        )
    )

    greater = int(
        values.size - less - equal
    )

    # Mid-rank empirical percentile with a half-observation
    # continuity correction. This prevents finite simulations
    # from producing literal zero or one.
    q = (
        less
        + 0.5 * equal
        + 0.5
    ) / (
        values.size + 1.0
    )

    compatibility = 2.0 * min(
        q,
        1.0 - q,
    )

    return {
        "role": role,
        "order": int(
            specification["order"]
        ),
        "scenario": specification["scenario"],
        "generating_model": (
            specification["generating_model"]
        ),
        "label": specification["label"],
        "addition": specification["addition"],
        "colour": specification["colour"],
        "simulations": int(values.size),
        "observed_LR": float(observed_lr),
        "simulated_LR_min": float(
            np.min(values)
        ),
        "simulated_LR_median": float(
            np.median(values)
        ),
        "simulated_LR_max": float(
            np.max(values)
        ),
        "simulated_LR_below_observed": less,
        "simulated_LR_equal_observed": equal,
        "simulated_LR_above_observed": greater,
        "observed_empirical_percentile": float(q),
        "two_sided_compatibility": float(
            compatibility
        ),
        "display_coordinate": float(
            math.sqrt(compatibility)
        ),
        "percentile_label": percentile_text(q),
    }


def load_results():
    require_file(SIMULATION_TABLE)
    require_file(OBSERVED_AUDIT)

    simulations = pd.read_csv(
        SIMULATION_TABLE,
        sep="\t",
        compression="gzip",
    )

    required = {
        "scenario",
        "generating_model",
        "LR",
        "sites",
    }

    missing = required - set(
        simulations.columns
    )

    if missing:
        raise RuntimeError(
            "Simulation table is missing columns: "
            f"{sorted(missing)}"
        )

    audit = json.loads(
        OBSERVED_AUDIT.read_text()
    )

    if "global_LR" not in audit:
        raise RuntimeError(
            "Observed audit lacks global_LR"
        )

    observed_lr = float(
        audit["global_LR"]
    )

    observed_sites = int(
        audit.get("sites", 0)
    )

    if not np.isfinite(observed_lr):
        raise RuntimeError(
            "Observed LR is not finite"
        )

    simulation_sites = sorted(
        pd.to_numeric(
            simulations["sites"],
            errors="coerce",
        )
        .dropna()
        .astype(int)
        .unique()
        .tolist()
    )

    if (
        observed_sites > 0
        and simulation_sites != [
            observed_sites
        ]
    ):
        raise RuntimeError(
            "Observed and simulated scopes differ: "
            f"observed={observed_sites}, "
            f"simulated={simulation_sites}"
        )

    rows = []

    for specification in MAIN_STAGES:
        rows.append(
            compatibility_row(
                simulations,
                observed_lr,
                specification,
                "main_sequence",
            )
        )

    rows.append(
        compatibility_row(
            simulations,
            observed_lr,
            SENSITIVITY_CONTROL,
            "sensitivity_control",
        )
    )

    results = pd.DataFrame(rows)

    return (
        results,
        observed_lr,
        observed_sites,
        simulation_sites,
    )


def make_plot(results):
    main = (
        results.loc[
            results["role"].eq(
                "main_sequence"
            )
        ]
        .sort_values("order")
        .reset_index(drop=True)
    )

    if len(main) != 3:
        raise RuntimeError(
            "Expected exactly three main stages"
        )

    y = np.arange(
        len(main) - 1,
        -1,
        -1,
        dtype=float,
    )

    x = main[
        "display_coordinate"
    ].to_numpy(dtype=float)

    colours = main[
        "colour"
    ].tolist()

    figure, axis = plt.subplots(
        figsize=(11.4, 5.4)
    )

    figure.subplots_adjust(
        left=0.36,
        right=0.96,
        top=0.77,
        bottom=0.22,
    )

    axis.plot(
        x,
        y,
        color="#B8B8B8",
        linewidth=2.0,
        zorder=1,
    )

    for index, row in main.iterrows():
        y_value = y[index]
        x_value = float(
            row["display_coordinate"]
        )

        axis.hlines(
            y_value,
            0.0,
            x_value,
            color=row["colour"],
            linewidth=3.0,
            alpha=0.55,
            zorder=2,
        )

        axis.scatter(
            [x_value],
            [y_value],
            s=150 if index < 2 else 205,
            color=row["colour"],
            edgecolor="white",
            linewidth=1.6,
            zorder=4,
        )

        if index == 2:
            axis.scatter(
                [x_value],
                [y_value],
                s=295,
                facecolor="none",
                edgecolor="#222222",
                linewidth=1.5,
                zorder=3,
            )

        compatibility = float(
            row["two_sided_compatibility"]
        )

        annotation = (
            f"C = {compatibility:.4f}"
            f"   |   observed at "
            f"{row['percentile_label']}"
            f"   |   n = "
            f"{int(row['simulations'])}"
        )

        text_x = min(
            x_value + 0.025,
            0.66,
        )

        axis.text(
            text_x,
            y_value,
            annotation,
            va="center",
            ha="left",
            fontsize=10.2,
            color="#303030",
        )

    tick_probabilities = np.array(
        [
            0.0,
            0.01,
            0.05,
            0.25,
            0.50,
            1.0,
        ]
    )

    axis.set_xticks(
        np.sqrt(tick_probabilities)
    )

    axis.set_xticklabels(
        [
            "0",
            "0.01",
            "0.05",
            "0.25",
            "0.50",
            "1",
        ]
    )

    axis.set_xlim(
        -0.015,
        1.02,
    )

    axis.set_ylim(
        -0.55,
        len(main) - 0.45,
    )

    axis.set_yticks(y)

    axis.set_yticklabels(
        main["label"].tolist(),
        fontsize=10.5,
    )

    axis.set_xlabel(
        "Empirical two-sided simulation compatibility, C"
        "\n(square-root display scale; tick labels show C)",
        fontsize=11,
    )

    axis.axvline(
        math.sqrt(0.5),
        color="#D0D0D0",
        linestyle="--",
        linewidth=1.0,
        zorder=0,
    )

    axis.axvline(
        1.0,
        color="#D0D0D0",
        linestyle=":",
        linewidth=1.0,
        zorder=0,
    )

    axis.grid(
        axis="x",
        color="#E8E8E8",
        linewidth=0.8,
        zorder=0,
    )

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)

    axis.tick_params(
        axis="y",
        length=0,
    )

    figure.text(
        0.055,
        0.925,
        "Do successive biological layers reproduce the observed founder evidence?",
        fontsize=18,
        fontweight="bold",
        ha="left",
        va="top",
        color="#202020",
    )

    figure.text(
        0.055,
        0.855,
        (
            "Compatibility approaches 1 when the observed LR lies "
            "near the centre of the simulated distribution; "
            "values near 0 indicate a tail mismatch."
        ),
        fontsize=10.8,
        ha="left",
        va="top",
        color="#555555",
    )

    figure.text(
        0.055,
        0.055,
        (
            "C = 2 × min(q, 1 − q), where q is the empirical "
            "percentile of the observed LR. This is a simulation-"
            "compatibility diagnostic, not a posterior probability."
        ),
        fontsize=9.2,
        ha="left",
        va="bottom",
        color="#666666",
    )

    OUTPUT_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure.savefig(
        OUTPUT_PDF,
        bbox_inches="tight",
    )

    figure.savefig(
        OUTPUT_PNG,
        dpi=240,
        bbox_inches="tight",
    )

    plt.close(figure)


def main(validate_only=False):
    (
        results,
        observed_lr,
        observed_sites,
        simulation_sites,
    ) = load_results()

    print(
        f"observed_LR={observed_lr}"
    )

    print(
        f"observed_sites={observed_sites}"
    )

    print(
        f"simulation_sites={simulation_sites}"
    )

    print(
        results[
            [
                "role",
                "order",
                "scenario",
                "generating_model",
                "simulations",
                "observed_empirical_percentile",
                "two_sided_compatibility",
            ]
        ].to_string(index=False)
    )

    if validate_only:
        print("VALIDATION_PASS")
        return

    OUTPUT_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    results.to_csv(
        OUTPUT_TABLE,
        sep="\t",
        index=False,
    )

    make_plot(results)

    main_results = results.loc[
        results["role"].eq(
            "main_sequence"
        )
    ].sort_values("order")

    audit = {
        "observed_LR": observed_lr,
        "observed_sites": observed_sites,
        "simulation_sites": simulation_sites,
        "compatibility_definition": (
            "2 * min(q, 1-q), with continuity-corrected "
            "empirical percentile q"
        ),
        "display_scale": (
            "square root; axis labels report untransformed "
            "compatibility"
        ),
        "main_stage_count": int(
            len(main_results)
        ),
        "best_main_stage": str(
            main_results.loc[
                main_results[
                    "two_sided_compatibility"
                ].idxmax(),
                "label",
            ]
        ),
        "final_stage_is_best": bool(
            main_results[
                "two_sided_compatibility"
            ].iloc[-1]
            == main_results[
                "two_sided_compatibility"
            ].max()
        ),
        "posterior_probability_assigned": False,
        "pdf": str(OUTPUT_PDF),
        "png": str(OUTPUT_PNG),
        "table": str(OUTPUT_TABLE),
    }

    OUTPUT_AUDIT.write_text(
        json.dumps(
            audit,
            indent=2,
        )
        + "\n"
    )

    print("STAGE_213_COMPLETE")

    print(
        "PDF\t"
        f"{OUTPUT_PDF.resolve()}"
    )

    print(
        "PNG\t"
        f"{OUTPUT_PNG.resolve()}"
    )

    print(
        "AUDIT\t"
        f"{OUTPUT_AUDIT.resolve()}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--validate-only",
        action="store_true",
    )

    arguments = parser.parse_args()

    main(
        validate_only=arguments.validate_only
    )
