import argparse
import csv
import gzip
import importlib.util
import itertools
import json
import math
import os
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

ROOT = Path("analysis/genomewide_mixing_corrected_v1")
BASE = ROOT / "expanded_calibration"
OUT = ROOT / "expanded_calibration_v2"

CONTROLS = """ERR3790122 ERR3790123 ERR2298257 ERR2298260
ERR3790118 ERR2298201 ERR2298265 ERR2298235 ERR2298236
ERR2298216 ERR1307113 ERR1307112 ERR1307111 ERR3790124
ERR2298208 ERR2298263 ERR2298221 SRR11243339""".split()

TASKS = 150
PER_TASK = 25
WINDOW = 1_000_000


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pilot = load_module(
    "scripts/185_restored_model_pilot.py",
    "pilot185",
)
pilot.OUT = OUT


def load_data(calibration=False):
    suffix = ".calibration.npz" if calibration else ".npz"
    pieces = []

    for chromosome in range(1, 21):
        with np.load(
            BASE / f"chr{chromosome}{suffix}",
            allow_pickle=False,
        ) as source:
            pieces.append({
                key: source[key]
                for key in source.files
            })

    return {
        key: np.concatenate([
            piece[key]
            for piece in pieces
        ])
        for key in pieces[0]
    }


def simulation_design():
    pairs = list(itertools.combinations(range(18), 2))
    rows = []
    identifier = 0

    global_grid = np.arange(0, 0.451, 0.05)
    focused_grid = np.arange(0.14, 0.261, 0.02)

    for family, grid, repetitions in [
        ("global", global_grid, 20),
        ("focused", focused_grid, 20),
    ]:
        for cell, (c_left, c_right) in enumerate(
            itertools.product(grid, grid)
        ):
            for repetition in range(repetitions):
                pair_index = (
                    cell * 37
                    + repetition * 11
                    + (0 if family == "global" else 71)
                ) % len(pairs)

                rows.append((
                    identifier,
                    family,
                    float(c_left),
                    float(c_right),
                    pair_index,
                    repetition,
                ))
                identifier += 1

    for pair_index in range(len(pairs)):
        for repetition in range(5):
            rows.append((
                identifier,
                "repeat_20",
                0.20,
                0.20,
                pair_index,
                repetition,
            ))
            identifier += 1

    if len(rows) != 3745:
        raise AssertionError(
            f"Expected 3745 configurations; found {len(rows)}"
        )

    return rows, pairs


def impute(matrix, index, allele_frequency, rng):
    genotype = matrix[:, index].copy()
    missing = genotype < 0
    genotype[missing] = rng.binomial(
        2,
        allele_frequency[missing],
    )
    return genotype


def gamete(genotype, rng):
    return rng.binomial(1, genotype / 2)


def run_simulations(task):
    data = load_data(calibration=True)
    matrix = data["controls"]
    called = matrix >= 0

    allele_frequency = (
        np.where(called, matrix, 0).sum(axis=1) + 0.5
    ) / (
        2 * called.sum(axis=1) + 1
    )

    schedule, pairs = simulation_design()
    start = task * PER_TASK
    stop = min(len(schedule), start + PER_TASK)

    if start >= len(schedule):
        raise ValueError(f"Task {task} lies outside schedule")

    rows = []
    emission = pilot.load_emission()

    for index in range(start, stop):
        (
            identifier,
            family,
            c_left,
            c_right,
            pair_index,
            draw,
        ) = schedule[index]

        donor_a_index, donor_b_index = pairs[pair_index]

        for founder_index, founder in enumerate([
            "three",
            "four",
        ]):
            rng = np.random.default_rng(
                187000000
                + identifier * 2
                + founder_index
            )

            donor_a = impute(
                matrix,
                donor_a_index,
                allele_frequency,
                rng,
            )
            donor_b = impute(
                matrix,
                donor_b_index,
                allele_frequency,
                rng,
            )

            if founder == "three":
                constructed = (
                    gamete(donor_a, rng)
                    + gamete(donor_b, rng)
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

                if np.any(
                    np.abs(true_left - true_right) == 2
                ):
                    raise AssertionError(
                        "Shared-allele construction violated"
                    )
            else:
                true_left, true_right = donor_a, donor_b

            synthetic = pilot.simulate_counts(
                true_left,
                true_right,
                c_left,
                c_right,
                data,
                emission,
                rng,
            )

            fit = pilot.fit_pair(
                synthetic,
                emission,
                f"v2_t{task}_i{identifier}_{founder}",
            )

            rows.append({
                "id": identifier,
                "family": family,
                "draw": draw,
                "generating_model": founder,
                "donor_a": CONTROLS[donor_a_index],
                "donor_b": CONTROLS[donor_b_index],
                "true_c_left": c_left,
                "true_c_right": c_right,
                "genotype_contrast": float(
                    np.mean(
                        np.abs(true_left - true_right) / 2
                    )
                ),
                "estimated_c_left_null":
                    fit["null"]["c_left"],
                "estimated_c_right_null":
                    fit["null"]["c_right"],
                "estimated_c_left_alt":
                    fit["alternative"]["c_left"],
                "estimated_c_right_alt":
                    fit["alternative"]["c_right"],
                "epsilon":
                    fit["alternative"]["epsilon"],
                "LR": fit["LR"],
                "sites": len(true_left),
            })

        if (index - start + 1) % 5 == 0:
            print(
                f"PROGRESS task={task} "
                f"configurations={index-start+1}",
                flush=True,
            )

    output = OUT / "tasks" / f"task_{task:03d}.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"COMPLETE task={task} "
        f"configurations={stop-start} rows={len(rows)}",
        flush=True,
    )


def collect_simulations():
    rows = []

    for task in range(TASKS):
        path = OUT / "tasks" / f"task_{task:03d}.tsv"
        with path.open() as handle:
            rows.extend(
                csv.DictReader(handle, delimiter="\t")
            )

    if len(rows) != 7490:
        raise RuntimeError(
            f"Expected 7490 rows; found {len(rows)}"
        )

    output = OUT / "simulations.tsv.gz"
    temporary = str(output) + ".new"

    with gzip.open(temporary, "wt") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    os.replace(temporary, output)

    counts = defaultdict(int)
    for row in rows:
        counts[
            row["family"]
            + "_"
            + row["generating_model"]
        ] += 1

    audit = {
        "rows": len(rows),
        "configurations": 3745,
        "counts": dict(counts),
        "global_grid":
            "0.00-0.45 by 0.05; 20 donor pairs per cell",
        "focused_grid":
            "0.14-0.26 by 0.02; 20 donor pairs per cell",
        "repeat_20":
            "all 153 donor pairs; five stochastic draws",
        "sites_per_simulation": 20000,
        "posterior_assigned": False,
    }

    (
        OUT / "simulation_design.audit.json"
    ).write_text(
        json.dumps(audit, indent=2) + "\n"
    )

    print(json.dumps(audit, indent=2), flush=True)


def calculate_autosomal_evidence():
    data = load_data(calibration=False)

    fit = json.loads(
        (
            BASE / "observed_expanded.json"
        ).read_text()
    )["full"]

    emission = pilot.load_emission()

    null_likelihood = (
        pilot.model.state_log_likelihoods(
            data,
            fit["null"]["c_left"],
            fit["null"]["c_right"],
            emission,
        )
    )

    null_evidence = pilot.model_evidence(
        null_likelihood,
        0,
    )

    epsilon = fit["alternative"]["epsilon"]

    alternative_likelihood = (
        pilot.model.state_log_likelihoods(
            data,
            fit["alternative"]["c_left"],
            fit["alternative"]["c_right"],
            emission,
        )
    )

    alternative_evidence = pilot.model_evidence(
        alternative_likelihood,
        epsilon,
    )

    compatible = (
        logsumexp(
            alternative_likelihood[
                :, pilot.COMPATIBLE
            ],
            axis=1,
        )
        - math.log(7)
        + math.log1p(-epsilon)
    )

    opposing = (
        logsumexp(
            alternative_likelihood[
                :, pilot.OPPOSING
            ],
            axis=1,
        )
        - math.log(2)
        + math.log(epsilon)
    )

    opposing_posterior = np.exp(
        opposing
        - np.logaddexp(
            compatible,
            opposing,
        )
    )

    site_lr = 2 * (
        alternative_evidence
        - null_evidence
    )

    groups = defaultdict(list)

    for index, (chromosome, position) in enumerate(
        zip(data["chrom"], data["pos"])
    ):
        start = (
            (int(position) - 1)
            // WINDOW
            * WINDOW
        )
        groups[
            (int(chromosome), start)
        ].append(index)

    rows = []

    for (
        chromosome,
        start,
    ), indices in sorted(groups.items()):
        selected = np.asarray(indices)

        rows.append({
            "chromosome": chromosome,
            "window_start": start,
            "window_end": start + WINDOW,
            "sites": len(selected),
            "LR_contribution": float(
                site_lr[selected].sum()
            ),
            "LR_per_site": float(
                site_lr[selected].mean()
            ),
            "opposing_posterior_sum": float(
                opposing_posterior[selected].sum()
            ),
            "opposing_posterior_mean": float(
                opposing_posterior[selected].mean()
            ),
            "mean_left_depth": float(
                data["left_total"][selected].mean()
            ),
            "mean_right_depth": float(
                data["right_total"][selected].mean()
            ),
        })

    output = OUT / "autosomal_1Mb_evidence.tsv"

    with output.open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    audit = {
        "windows": len(rows),
        "sites": len(site_lr),
        "global_LR": fit["LR"],
        "window_LR_sum": float(site_lr.sum()),
        "difference": float(
            site_lr.sum() - fit["LR"]
        ),
        "interpretation": (
            "Fixed-global-parameter site contributions; "
            "not independent window tests"
        ),
    }

    (
        OUT / "autosomal_1Mb_evidence.audit.json"
    ).write_text(
        json.dumps(audit, indent=2) + "\n"
    )

    print(json.dumps(audit, indent=2), flush=True)


def prepare_z_diagnostic():
    source = (
        ROOT
        / "inputs"
        / "chr21_NC_140409.1.qc.tsv.gz"
    )

    rng = np.random.default_rng(18721)
    reservoir = []
    eligible = 0

    with gzip.open(source, "rt") as handle:
        header = next(handle).rstrip().split("\t")
        left_index = header.index("TN_Left_merged")
        right_index = header.index("TN_Right_merged")
        control_indices = [
            header.index(sample)
            for sample in CONTROLS
        ]

        for line in handle:
            row = line.rstrip().split("\t")

            try:
                left_alt, left_total, left_dp = (
                    pilot.parse_focal_counts(
                        row[left_index]
                    )
                )
                right_alt, right_total, right_dp = (
                    pilot.parse_focal_counts(
                        row[right_index]
                    )
                )
            except ValueError:
                continue

            if min(
                left_total,
                right_total,
                left_dp,
                right_dp,
            ) < 10:
                continue

            calls = [
                pilot.qc.classify(row[index])[1]
                for index in control_indices
            ]

            dosage = np.array([
                value[0] if value is not None else -1
                for value in calls
            ])
            called = dosage >= 0

            if called.sum() < 13:
                continue

            alternate = dosage[called].sum()

            if (
                alternate == 0
                or alternate == 2 * called.sum()
            ):
                continue

            record = (
                int(row[1]),
                left_alt,
                left_total,
                right_alt,
                right_total,
            )
            eligible += 1

            if len(reservoir) < 50000:
                reservoir.append(record)
            else:
                replacement = int(
                    rng.integers(0, eligible)
                )
                if replacement < 50000:
                    reservoir[replacement] = record

    reservoir.sort()
    array = np.asarray(reservoir)

    output = OUT / "z_diagnostic.npz"
    temporary = str(output) + ".new"

    with open(temporary, "wb") as handle:
        np.savez_compressed(
            handle,
            pos=array[:, 0],
            left_alt=array[:, 1],
            left_total=array[:, 2],
            right_alt=array[:, 3],
            right_total=array[:, 4],
        )

    os.replace(temporary, output)

    audit = {
        "eligible": eligible,
        "retained": len(array),
        "left_mean_depth": float(
            array[:, 2].mean()
        ),
        "right_mean_depth": float(
            array[:, 4].mean()
        ),
        "left_to_right_depth_ratio": float(
            array[:, 2].mean()
            / array[:, 4].mean()
        ),
        "usage": (
            "sex-ploidy diagnostic only; "
            "not diploid founder fitting"
        ),
    }

    (
        OUT / "z_diagnostic.audit.json"
    ).write_text(
        json.dumps(audit, indent=2) + "\n"
    )

    print(json.dumps(audit, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "mode",
        choices=[
            "simulate",
            "collect",
            "evidence",
            "zdiag",
        ],
    )
    parser.add_argument("--task", type=int, default=0)
    arguments = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()

    print(
        f"START mode={arguments.mode} "
        f"task={arguments.task}",
        flush=True,
    )

    if arguments.mode == "simulate":
        run_simulations(arguments.task)
    elif arguments.mode == "collect":
        collect_simulations()
    elif arguments.mode == "evidence":
        calculate_autosomal_evidence()
    else:
        prepare_z_diagnostic()

    print(
        f"END elapsed_seconds="
        f"{time.monotonic()-started:.1f}",
        flush=True,
    )
