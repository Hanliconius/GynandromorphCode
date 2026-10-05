import argparse
import csv
import gzip
import importlib.util
import json
import math
import os
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


ROOT = Path("analysis/genomewide_mixing_corrected_v1")
OUT = ROOT / "restored_pilot"

COMPATIBLE = np.array([0, 1, 3, 4, 5, 7, 8])
OPPOSING = np.array([2, 6])


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


qc = load_module(
    "scripts/184_control_qc_and_emission.py",
    "control_qc_184",
)

model = load_module(
    "scripts/159_fit_autosomal_three_vs_four_founder.py",
    "founder_model_159",
)


def load_emission():
    path = (
        ROOT
        / "control_qc"
        / "emission_calibration.diagnostic.json"
    )
    values = json.loads(path.read_text())
    if len(values) != 2:
        raise ValueError("Expected two emission diagnostic scopes")

    all_biallelic, snps_only = values

    for key in [
        "error_rate",
        "bias_log_odds",
        "concentration",
    ]:
        if not np.isclose(
            all_biallelic[key],
            snps_only[key],
            rtol=1e-12,
            atol=1e-12,
        ):
            raise ValueError(
                "SNP-only and all-biallelic calibration differ; "
                "inspect before selecting one"
            )

    return {
        key: float(all_biallelic[key])
        for key in [
            "error_rate",
            "bias_log_odds",
            "concentration",
        ]
    }


def save_npz(path, **data):
    temporary = str(path) + ".new"

    with open(temporary, "wb") as handle:
        np.savez_compressed(handle, **data)

    os.replace(temporary, path)


def parse_focal_counts(cell):
    fields = cell.split(":")
    if len(fields) != 5:
        raise ValueError("Expected GT:AD:DP:PL:SP")

    gt, ad, dp, pl, sp = fields

    try:
        depths = [int(value) for value in ad.split(",")]
        format_dp = int(dp)
    except ValueError as error:
        raise ValueError(
            f"Invalid focal AD or DP: {cell}"
        ) from error

    if len(depths) != 2 or min(depths) < 0:
        raise ValueError(f"Invalid focal AD vector: {cell}")

    return depths[1], sum(depths), format_dp


def prepare_chromosome(chromosome):
    if not 1 <= chromosome <= 20:
        raise ValueError("Pilot preparation uses autosomes 1-20")

    old_path = Path(
        "analysis/genomewide_mixing/observed_fit/"
        "autosomal_fit_input.npz"
    )

    with np.load(old_path, allow_pickle=False) as old:
        indices = np.flatnonzero(old["chrom"] == chromosome)

        wanted = {
            int(old["pos"][index]): tuple(
                int(old[key][index])
                for key in [
                    "left_alt",
                    "left_dp",
                    "right_alt",
                    "right_dp",
                ]
            )
            for index in indices
        }

    if not wanted or len(wanted) != len(indices):
        raise ValueError(
            "Missing or duplicated positions in previous fit input"
        )

    contig = f"NC_{140388 + chromosome}.1"
    source = (
        ROOT
        / "inputs"
        / f"chr{chromosome}_{contig}.qc.tsv.gz"
    )

    rows = []
    seen = set()
    rejected = 0

    with gzip.open(source, "rt") as handle:
        header = next(handle).rstrip("\n").split("\t")

        control_indices = [
            header.index(sample)
            for sample in qc.CONTROLS
        ]

        left_index = header.index("TN_Left_merged")
        right_index = header.index("TN_Right_merged")

        for line in handle:
            row = line.rstrip("\n").split("\t")
            position = int(row[1])

            if position not in wanted:
                continue

            if position in seen:
                raise ValueError(
                    f"Ambiguous position: chr{chromosome}:{position}"
                )

            seen.add(position)

            left_alt, left_total, left_dp = (
                parse_focal_counts(row[left_index])
            )
            right_alt, right_total, right_dp = (
                parse_focal_counts(row[right_index])
            )

            old_counts = wanted[position]

            if (
                left_alt,
                left_dp,
                right_alt,
                right_dp,
            ) != old_counts:
                raise ValueError(
                    "Previous and re-extracted counts disagree at "
                    f"chr{chromosome}:{position}"
                )

            if min(
                left_total,
                right_total,
                left_dp,
                right_dp,
            ) < 10:
                rejected += 1
                continue

            control_values = [
                qc.classify(row[index])[1]
                for index in control_indices
            ]

            control_dosage = [
                value[0] if value is not None else -1
                for value in control_values
            ]

            rows.append(
                (
                    position,
                    row[2],
                    row[3],
                    left_alt,
                    left_total,
                    right_alt,
                    right_total,
                    control_dosage,
                )
            )

    if seen != set(wanted):
        missing = len(set(wanted) - seen)
        raise ValueError(
            f"{missing} previous positions were absent from extraction"
        )

    if not rows:
        raise ValueError(
            f"No sites retained for chromosome {chromosome}"
        )

    output = OUT / f"chr{chromosome}.npz"

    save_npz(
        output,
        chrom=np.full(len(rows), chromosome),
        pos=np.array([row[0] for row in rows]),
        ref=np.array([row[1] for row in rows]),
        alt=np.array([row[2] for row in rows]),
        left_alt=np.array([row[3] for row in rows]),
        left_total=np.array([row[4] for row in rows]),
        right_alt=np.array([row[5] for row in rows]),
        right_total=np.array([row[6] for row in rows]),
        controls=np.array(
            [row[7] for row in rows],
            dtype=np.int8,
        ),
    )

    audit = {
        "chromosome": chromosome,
        "source": str(source),
        "previous_fit_input": str(old_path),
        "old_positions": len(wanted),
        "retained": len(rows),
        "rejected_AD_or_DP_below_10": rejected,
        "inherited_previous_site_selection": True,
        "focal_depth_used_for_likelihood": "AD_REF_plus_AD_ALT",
        "focal_depth_cap": None,
        "control_confidence_rule": (
            "DP>=2; sumAD>=2; GT at unique PL minimum; "
            "next-best minus best PL>=20"
        ),
    }

    (
        OUT / f"chr{chromosome}.audit.json"
    ).write_text(json.dumps(audit, indent=2) + "\n")

    print(
        f"PREPARED chromosome={chromosome} "
        f"retained={len(rows)} rejected={rejected}",
        flush=True,
    )


def load_dataset():
    pieces = []

    for chromosome in range(1, 21):
        path = OUT / f"chr{chromosome}.npz"

        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing prepared input: {path}")

        with np.load(path, allow_pickle=False) as data:
            pieces.append(
                {
                    key: data[key]
                    for key in data.files
                }
            )

    keys = pieces[0].keys()

    return {
        key: np.concatenate(
            [piece[key] for piece in pieces]
        )
        for key in keys
    }


def model_evidence(state_likelihood, epsilon):
    compatible = (
        logsumexp(
            state_likelihood[:, COMPATIBLE],
            axis=1,
        )
        - math.log(7)
    )

    if epsilon == 0:
        return compatible

    opposing = (
        logsumexp(
            state_likelihood[:, OPPOSING],
            axis=1,
        )
        - math.log(2)
    )

    return np.logaddexp(
        compatible + math.log1p(-epsilon),
        opposing + math.log(epsilon),
    )


def fit_pair(data, emission, label):
    optimizer_trace = []

    def run_model(alternative, starts):
        def objective(parameters):
            state_likelihood = model.state_log_likelihoods(
                data,
                float(parameters[0]),
                float(parameters[1]),
                emission,
            )

            epsilon = (
                float(parameters[2])
                if alternative
                else 0.0
            )

            return -float(
                model_evidence(
                    state_likelihood,
                    epsilon,
                ).sum()
            )

        bounds = [
            (0.0, 0.45),
            (0.0, 0.45),
        ]

        if alternative:
            bounds.append((0.0, 0.5))

        successful = []

        for start in starts:
            result = minimize(
                objective,
                start,
                method="L-BFGS-B",
                bounds=bounds,
                options={
                    "maxiter": 350,
                    "ftol": 1e-10,
                },
            )

            optimizer_trace.append(
                {
                    "alternative": alternative,
                    "method": "L-BFGS-B",
                    "success": bool(result.success),
                    "message": str(result.message),
                    "objective": float(result.fun),
                    "parameters": result.x.tolist(),
                }
            )

            if result.success and np.isfinite(result.fun):
                successful.append(result)

        if not successful:
            fallback = minimize(
                objective,
                starts[0],
                method="Powell",
                bounds=bounds,
                options={
                    "maxiter": 750,
                    "ftol": 1e-9,
                    "xtol": 1e-6,
                },
            )

            optimizer_trace.append(
                {
                    "alternative": alternative,
                    "method": "Powell",
                    "success": bool(fallback.success),
                    "message": str(fallback.message),
                    "objective": float(fallback.fun),
                    "parameters": fallback.x.tolist(),
                }
            )

            if (
                fallback.success
                and np.isfinite(fallback.fun)
            ):
                successful.append(fallback)

        if not successful:
            (
                OUT / f"{label}.optimizer.json"
            ).write_text(
                json.dumps(
                    optimizer_trace,
                    indent=2,
                )
                + "\n"
            )

            raise RuntimeError(
                f"No converged optimizer for {label}; "
                "inspect optimizer trace"
            )

        best = min(
            successful,
            key=lambda result: result.fun,
        )

        return {
            "c_left": float(best.x[0]),
            "c_right": float(best.x[1]),
            "epsilon": (
                float(best.x[2])
                if alternative
                else 0.0
            ),
            "logL": -float(best.fun),
        }

    null = run_model(
        False,
        [
            (0.05, 0.05),
            (0.18, 0.18),
            (0.35, 0.35),
            (0.10, 0.30),
            (0.30, 0.10),
        ],
    )

    alternative = run_model(
        True,
        [
            (
                null["c_left"],
                null["c_right"],
                epsilon,
            )
            for epsilon in [
                0.001,
                0.05,
                0.20,
            ]
        ],
    )

    (
        OUT / f"{label}.optimizer.json"
    ).write_text(
        json.dumps(
            optimizer_trace,
            indent=2,
        )
        + "\n"
    )

    # The alternative parameter space contains the null at epsilon=0.
    # Retain the exact boundary solution if numerical optimization is
    # microscopically worse.
    if alternative["logL"] < null["logL"]:
        alternative = dict(null)

    result = {
        "label": label,
        "sites": len(data["left_alt"]),
        "null": null,
        "alternative": alternative,
        "LR": 2.0 * (
            alternative["logL"]
            - null["logL"]
        ),
        "emission": emission,
    }

    print(json.dumps(result), flush=True)
    return result


def fit_observed():
    data = load_dataset()
    emission = load_emission()

    result = fit_pair(
        data,
        emission,
        "observed",
    )

    result["scope"] = (
        "bridge_on_previous_selected_sites; "
        "not calibrated founder support"
    )

    (
        OUT / "observed.json"
    ).write_text(
        json.dumps(result, indent=2) + "\n"
    )


def simulate_counts(
    true_left,
    true_right,
    c_left,
    c_right,
    data,
    emission,
    rng,
):
    output = {}

    expected = {
        "left": (
            (1.0 - c_left) * true_left / 2.0
            + c_left * true_right / 2.0
        ),
        "right": (
            c_right * true_left / 2.0
            + (1.0 - c_right) * true_right / 2.0
        ),
    }

    for side in ["left", "right"]:
        mean = model.map_fraction(
            expected[side],
            emission["error_rate"],
            emission["bias_log_odds"],
        )

        concentration = emission["concentration"]

        latent = rng.beta(
            mean * concentration,
            (1.0 - mean) * concentration,
        )

        depth = data[f"{side}_total"]

        output[f"{side}_total"] = depth
        output[f"{side}_alt"] = rng.binomial(
            depth,
            latent,
        )

    return output


def run_pilot(task):
    if not 0 <= task <= 8:
        raise ValueError("Pilot task must be 0-8")

    data = load_dataset()
    rng = np.random.default_rng(185000 + task)

    # Fixed chromosome-balanced site set shared by all pilot tasks.
    selector = np.random.default_rng(185)

    keep = np.concatenate(
        [
            selector.choice(
                np.flatnonzero(
                    data["chrom"] == chromosome
                ),
                size=150,
                replace=False,
            )
            for chromosome in range(1, 21)
        ]
    )

    data = {
        key: value[keep]
        for key, value in data.items()
    }

    matrix = data["controls"]
    called = matrix >= 0

    allele_frequency = (
        np.where(called, matrix, 0).sum(axis=1)
        + 0.5
    ) / (
        2.0 * called.sum(axis=1)
        + 1.0
    )

    call_rate = called.mean(axis=0)

    heterozygosity = (
        (matrix == 1).sum(axis=0)
        / np.maximum(1, called.sum(axis=0))
    )

    order = np.argsort(
        heterozygosity,
        kind="stable",
    )

    donor_a_index = int(order[task])
    donor_b_index = int(order[task + 9])

    def donor(index):
        genotype = matrix[:, index].copy()
        missing = genotype < 0
        genotype[missing] = rng.binomial(
            2,
            allele_frequency[missing],
        )
        return genotype

    donor_a = donor(donor_a_index)
    donor_b = donor(donor_b_index)

    def gamete(genotype):
        return rng.binomial(
            1,
            genotype / 2.0,
        )

    constructed = (
        gamete(donor_a)
        + gamete(donor_b)
    )

    emission = load_emission()
    results = []

    panels = {
        key: data[key]
        for key in [
            "left_alt",
            "left_total",
            "right_alt",
            "right_total",
        ]
    }

    for founder in [
        "three",
        "four",
    ]:
        if founder == "three":
            true_left = donor_a
            true_right = constructed

            if task % 2:
                true_left, true_right = (
                    true_right,
                    true_left,
                )
        else:
            true_left = donor_a
            true_right = donor_b

        if (
            founder == "three"
            and np.any(
                np.abs(
                    true_left - true_right
                )
                == 2
            )
        ):
            raise AssertionError(
                "Shared-allele construction was violated"
            )

        synthetic = simulate_counts(
            true_left,
            true_right,
            0.20,
            0.20,
            data,
            emission,
            rng,
        )

        result = fit_pair(
            synthetic,
            emission,
            f"pilot{task}_{founder}",
        )

        result.update(
            generating_model=founder,
            true_c_left=0.20,
            true_c_right=0.20,
            donor_a=qc.CONTROLS[
                donor_a_index
            ],
            donor_b=qc.CONTROLS[
                donor_b_index
            ],
            donor_a_call_rate=float(
                call_rate[donor_a_index]
            ),
            donor_b_call_rate=float(
                call_rate[donor_b_index]
            ),
            mean_genotype_contrast=float(
                np.mean(
                    np.abs(
                        true_left
                        - true_right
                    )
                    / 2.0
                )
            ),
            imputation=(
                "original per-site binomial fallback; "
                "audited call rates"
            ),
        )

        results.append(result)

        for key, value in synthetic.items():
            panels[
                f"{founder}_{key}"
            ] = value

    save_npz(
        OUT / f"pilot{task}.panels.npz",
        **panels,
    )

    (
        OUT / f"pilot{task}.json"
    ).write_text(
        json.dumps(results, indent=2) + "\n"
    )


def make_plots():
    import matplotlib
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LogNorm

    observed_result = json.loads(
        (OUT / "observed.json").read_text()
    )

    summary = []

    pdf_path = (
        OUT
        / "restored_pilot_AF_comparison.pdf"
    )

    with PdfPages(pdf_path) as pdf:
        for task in range(9):
            results = json.loads(
                (
                    OUT
                    / f"pilot{task}.json"
                ).read_text()
            )

            for result in results:
                summary.append(
                    {
                        "task": task,
                        "generating_model": result[
                            "generating_model"
                        ],
                        "donor_a": result["donor_a"],
                        "donor_b": result["donor_b"],
                        "call_rate_a": result[
                            "donor_a_call_rate"
                        ],
                        "call_rate_b": result[
                            "donor_b_call_rate"
                        ],
                        "contrast": result[
                            "mean_genotype_contrast"
                        ],
                        "true_mixing": 0.20,
                        "null_left": result[
                            "null"
                        ]["c_left"],
                        "null_right": result[
                            "null"
                        ]["c_right"],
                        "alternative_left": result[
                            "alternative"
                        ]["c_left"],
                        "alternative_right": result[
                            "alternative"
                        ]["c_right"],
                        "epsilon": result[
                            "alternative"
                        ]["epsilon"],
                        "LR": result["LR"],
                    }
                )

            with np.load(
                OUT
                / f"pilot{task}.panels.npz",
                allow_pickle=False,
            ) as source:
                data = {
                    key: source[key]
                    for key in source.files
                }

            figure, axes = plt.subplots(
                1,
                3,
                figsize=(11.7, 4.8),
                sharex=True,
                sharey=True,
            )

            figure.subplots_adjust(
                left=0.07,
                right=0.89,
                bottom=0.18,
                top=0.72,
                wspace=0.15,
            )

            artists = []

            for (
                axis,
                prefix,
                title,
            ) in zip(
                axes,
                [
                    "",
                    "three_",
                    "four_",
                ],
                [
                    "Observed L/R",
                    "Simulated three-founder",
                    "Simulated four-founder",
                ],
            ):
                artist = axis.hexbin(
                    (
                        data[
                            f"{prefix}left_alt"
                        ]
                        / data[
                            f"{prefix}left_total"
                        ]
                    ),
                    (
                        data[
                            f"{prefix}right_alt"
                        ]
                        / data[
                            f"{prefix}right_total"
                        ]
                    ),
                    gridsize=43,
                    extent=(0, 1, 0, 1),
                    mincnt=1,
                    cmap="Greys",
                    linewidths=0,
                )

                artists.append(artist)

                axis.set_xlim(-0.02, 1.02)
                axis.set_ylim(-0.02, 1.02)
                axis.set_aspect(
                    "equal",
                    adjustable="box",
                )
                axis.set_xlabel(
                    "Left ALT-read fraction"
                )
                axis.set_title(
                    title,
                    loc="left",
                    fontweight="bold",
                    fontsize=10,
                )
                axis.spines["top"].set_visible(
                    False
                )
                axis.spines["right"].set_visible(
                    False
                )

                if prefix:
                    mixing = {
                        "c_left": 0.20,
                        "c_right": 0.20,
                    }
                else:
                    mixing = observed_result[
                        "alternative"
                    ]

                emission = load_emission()

                expected_left = (
                    (1.0 - mixing["c_left"])
                    * model.GL
                    + mixing["c_left"]
                    * model.GR
                )

                expected_right = (
                    mixing["c_right"]
                    * model.GL
                    + (
                        1.0
                        - mixing["c_right"]
                    )
                    * model.GR
                )

                x = model.map_fraction(
                    expected_left,
                    emission["error_rate"],
                    emission["bias_log_odds"],
                )

                y = model.map_fraction(
                    expected_right,
                    emission["error_rate"],
                    emission["bias_log_odds"],
                )

                axis.scatter(
                    x[COMPATIBLE],
                    y[COMPATIBLE],
                    s=28,
                    marker="s",
                    facecolor="none",
                    edgecolor="#0072B2",
                    linewidth=1.1,
                    zorder=4,
                )

                if prefix != "three_":
                    axis.scatter(
                        x[OPPOSING],
                        y[OPPOSING],
                        s=30,
                        marker="D",
                        facecolor="none",
                        edgecolor="#D55E00",
                        linewidth=1.1,
                        zorder=4,
                    )

            maximum = max(
                float(
                    artist.get_array().max()
                )
                for artist in artists
            )

            shared_norm = LogNorm(
                vmin=1,
                vmax=max(2, maximum),
            )

            for artist in artists:
                artist.set_norm(shared_norm)

            colorbar_axis = figure.add_axes(
                [
                    0.92,
                    0.23,
                    0.015,
                    0.42,
                ]
            )

            figure.colorbar(
                artists[0],
                cax=colorbar_axis,
                label="Sites per hexagon",
            )

            axes[0].set_ylabel(
                "Right ALT-read fraction"
            )

            representative = results[0]

            figure.suptitle(
                (
                    f"Matched-depth benchmark "
                    f"{task + 1}: "
                    f"{representative['donor_a']} / "
                    f"{representative['donor_b']}"
                ),
                fontsize=13,
                fontweight="bold",
                y=0.96,
            )

            figure.text(
                0.07,
                0.83,
                (
                    "3,000 sites; uncapped AD depths; "
                    "simulations at 20% in each direction.\n"
                    "Donor call rates: "
                    f"{representative['donor_a_call_rate']:.1%}, "
                    f"{representative['donor_b_call_rate']:.1%}; "
                    "missing donor GTs imputed."
                ),
                fontsize=9,
            )

            figure.text(
                0.07,
                0.05,
                (
                    "Blue squares: shared-allele-compatible states. "
                    "Orange diamonds: opposing homozygotes.\n"
                    "Pilot only: no calibrated founder probability "
                    "or mixing confidence intervals."
                ),
                fontsize=8,
            )

            pdf.savefig(
                figure,
                bbox_inches="tight",
            )
            plt.close(figure)

    recovery_path = (
        OUT
        / "pilot_recovery.tsv"
    )

    with recovery_path.open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(summary[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary)

    audit = {
        "observed_scope": (
            "previous selected sites"
        ),
        "observed_site_selection_reused": True,
        "focal_depth_denominator": (
            "AD_REF_plus_AD_ALT"
        ),
        "focal_depth_cap": None,
        "pilot_tasks": 9,
        "pilot_results": 18,
        "sites_per_pilot": 3000,
        "true_mixing": 0.20,
        "founder_support_calibrated": False,
        "confidence_intervals_calculated": False,
        "pdf": str(pdf_path),
        "recovery_table": str(
            recovery_path
        ),
    }

    (
        OUT
        / "restored_pilot.audit.json"
    ).write_text(
        json.dumps(audit, indent=2)
        + "\n"
    )

    print(
        "ENDPOINT=RESTORED_FIT_AND_PILOT_READY",
        flush=True,
    )


def run_checks():
    emission = load_emission()

    toy = {
        "left_alt": np.array(
            [0, 5, 10]
        ),
        "left_total": np.array(
            [10, 10, 10]
        ),
        "right_alt": np.array(
            [10, 5, 0]
        ),
        "right_total": np.array(
            [10, 10, 10]
        ),
    }

    state_likelihood = (
        model.state_log_likelihoods(
            toy,
            0.20,
            0.20,
            emission,
        )
    )

    if (
        state_likelihood.shape != (3, 9)
        or not np.isfinite(
            state_likelihood
        ).all()
    ):
        raise RuntimeError(
            "Model likelihood API check failed"
        )

    explicit_null = logsumexp(
        state_likelihood[:, COMPATIBLE]
        - math.log(7),
        axis=1,
    )

    nested_null = model_evidence(
        state_likelihood,
        0.0,
    )

    if not np.allclose(
        explicit_null,
        nested_null,
        rtol=1e-12,
        atol=1e-12,
    ):
        raise RuntimeError(
            "Seven-state null boundary check failed"
        )

    print(
        "PASS: model API, finite likelihood, "
        "exact seven-state null boundary",
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mode",
        choices=[
            "check",
            "prepare",
            "observed",
            "pilot",
            "plot",
        ],
    )
    parser.add_argument(
        "--task",
        type=int,
        default=0,
    )
    arguments = parser.parse_args()

    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    started = time.monotonic()

    print(
        f"START mode={arguments.mode} "
        f"task={arguments.task}",
        flush=True,
    )

    if arguments.mode == "check":
        run_checks()
    elif arguments.mode == "prepare":
        prepare_chromosome(
            arguments.task + 1
        )
    elif arguments.mode == "observed":
        fit_observed()
    elif arguments.mode == "pilot":
        run_pilot(arguments.task)
    else:
        make_plots()

    print(
        "END elapsed_seconds="
        f"{time.monotonic() - started:.1f}",
        flush=True,
    )
