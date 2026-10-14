#!/usr/bin/env python3

import argparse
import csv
import importlib.util
import itertools
import json
import os
from pathlib import Path

import numpy as np


ROOT = Path("analysis/genomewide_mixing_corrected_v1")
BASE = ROOT / "expanded_calibration"
FAMILY = ROOT / "full_sib_validation"
OUT = ROOT / "empirical_full_sib_calibration_v1"
INPUT = OUT / "inputs"
TASK_DIRECTORY = OUT / "tasks"

SAMPLE_TABLE = (
    FAMILY
    / "metadata/bam_set/full_sib_14_samples.tsv"
)

IBD_TABLE = (
    FAMILY
    / "ibd_transmission/genomewide/"
    "sibling_pair_IBD_summary.tsv"
)

OBSERVED = (
    BASE
    / "observed_expanded.json"
)

REPETITIONS = 5


def load_module(path, name):
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


pilot = load_module(
    "scripts/185_restored_model_pilot.py",
    "pilot185_empirical_siblings",
)


def offspring_samples():
    with SAMPLE_TABLE.open() as handle:
        rows = list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )

    offspring = [
        row["sample"]
        for row in rows
        if row["role"] == "F2_offspring"
    ]

    if len(offspring) != 12:
        raise RuntimeError(
            f"Expected 12 offspring; "
            f"found {len(offspring)}"
        )

    return offspring


def simulation_schedule():
    offspring = offspring_samples()

    pairs = list(
        itertools.combinations(
            range(12),
            2,
        )
    )

    if len(pairs) != 66:
        raise RuntimeError(
            f"Expected 66 sibling pairs; "
            f"found {len(pairs)}"
        )

    observed = json.loads(
        OBSERVED.read_text()
    )["full"]["alternative"]

    scenarios = [
        (
            "exact_20",
            0.20,
            0.20,
        ),
        (
            "observed_matched",
            float(observed["c_left"]),
            float(observed["c_right"]),
        ),
    ]

    rows = []
    identifier = 0

    for (
        scenario,
        c_left,
        c_right,
    ) in scenarios:
        for pair_index, (
            offspring_a,
            offspring_b,
        ) in enumerate(pairs):
            for repetition in range(
                REPETITIONS
            ):
                rows.append(
                    {
                        "id": identifier,
                        "scenario": scenario,
                        "c_left": c_left,
                        "c_right": c_right,
                        "pair_index": pair_index,
                        "offspring_a_index":
                            offspring_a,
                        "offspring_b_index":
                            offspring_b,
                        "offspring_a":
                            offspring[offspring_a],
                        "offspring_b":
                            offspring[offspring_b],
                        "repetition": repetition,
                    }
                )

                identifier += 1

    if len(rows) != 660:
        raise RuntimeError(
            f"Expected 660 configurations; "
            f"found {len(rows)}"
        )

    return rows


def load_input():
    pieces = []

    for chromosome in range(1, 21):
        path = (
            INPUT
            / f"chr{chromosome}.npz"
        )

        with np.load(
            path,
            allow_pickle=False,
        ) as data:
            pieces.append(
                {
                    key: data[key]
                    for key in data.files
                    if key != "offspring"
                }
            )

    combined = {
        key: np.concatenate(
            [
                piece[key]
                for piece in pieces
            ]
        )
        for key in pieces[0]
    }

    if len(combined["pos"]) != 1_000_000:
        raise RuntimeError(
            "Expected 1,000,000 combined sites; "
            f"found {len(combined['pos'])}"
        )

    if combined["siblings"].shape != (
        1_000_000,
        12,
    ):
        raise RuntimeError(
            "Unexpected sibling genotype-matrix shape: "
            f"{combined['siblings'].shape}"
        )

    return combined


def sibling_ibd():
    values = {}

    with IBD_TABLE.open() as handle:
        for row in csv.DictReader(
            handle,
            delimiter="\t",
        ):
            key = tuple(
                sorted(
                    (
                        row["offspring_a"],
                        row["offspring_b"],
                    )
                )
            )

            values[key] = {
                "IBD0_fraction":
                    float(row["IBD0_fraction"]),
                "IBD1_fraction":
                    float(row["IBD1_fraction"]),
                "IBD2_fraction":
                    float(row["IBD2_fraction"]),
            }

    if len(values) != 66:
        raise RuntimeError(
            f"Expected 66 sibling IBD rows; "
            f"found {len(values)}"
        )

    return values


def impute(
    genotype,
    allele_frequency,
    rng,
):
    output = genotype.copy()
    missing = output < 0

    output[missing] = rng.binomial(
        2,
        allele_frequency[missing],
    )

    return output


def gamete(genotype, rng):
    return rng.binomial(
        1,
        genotype / 2,
    )


def main(task):
    TASK_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    output = (
        TASK_DIRECTORY
        / f"task_{task:03d}.tsv"
    )

    completion = (
        TASK_DIRECTORY
        / f"task_{task:03d}.complete.tsv"
    )

    if (
        output.is_file()
        and output.stat().st_size > 0
        and completion.is_file()
        and completion.stat().st_size > 0
    ):
        print(
            f"REUSE task={task} output={output}"
        )

        return

    schedule = simulation_schedule()

    if not 0 <= task < len(schedule):
        raise RuntimeError(
            f"Task {task} is outside the "
            f"0-{len(schedule) - 1} schedule"
        )

    configuration = schedule[task]
    data = load_input()

    controls = data["controls"]
    called = controls >= 0
    called_count = called.sum(axis=1)

    if np.any(called_count == 0):
        raise RuntimeError(
            "A model site has no called control genotype"
        )

    allele_frequency = (
        np.where(
            called,
            controls,
            0,
        ).sum(axis=1)
        + 0.5
    ) / (
        2 * called_count
        + 1
    )

    left_index = (
        configuration["offspring_a_index"]
    )

    right_index = (
        configuration["offspring_b_index"]
    )

    left_sample = (
        configuration["offspring_a"]
    )

    right_sample = (
        configuration["offspring_b"]
    )

    # Reverse odd replicates so both sibling orientations
    # occur on the left and right.
    if configuration["repetition"] % 2 == 1:
        left_index, right_index = (
            right_index,
            left_index,
        )

        left_sample, right_sample = (
            right_sample,
            left_sample,
        )

    ibd_key = tuple(
        sorted(
            (
                configuration["offspring_a"],
                configuration["offspring_b"],
            )
        )
    )

    ibd = sibling_ibd()[ibd_key]
    emission = pilot.load_emission()

    optimizer_directory = (
        TASK_DIRECTORY
        / f"optimizer_{task:03d}"
    )

    optimizer_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    pilot.OUT = optimizer_directory
    result_rows = []

    for founder_index, generating_model in enumerate(
        (
            "three_founder",
            "empirical_full_sib",
        )
    ):
        rng = np.random.default_rng(
            201000000
            + task * 2
            + founder_index
        )

        sibling_left = impute(
            data["siblings"][:, left_index],
            allele_frequency,
            rng,
        )

        sibling_right = impute(
            data["siblings"][:, right_index],
            allele_frequency,
            rng,
        )

        if generating_model == "three_founder":
            constructed = (
                gamete(
                    sibling_left,
                    rng,
                )
                + gamete(
                    sibling_right,
                    rng,
                )
            )

            if task % 2 == 0:
                true_left = sibling_left
                true_right = constructed
            else:
                true_left = constructed
                true_right = sibling_left

            if np.any(
                np.abs(
                    true_left
                    - true_right
                )
                == 2
            ):
                raise RuntimeError(
                    "The restored three-founder "
                    "shared-allele construction failed"
                )
        else:
            true_left = sibling_left
            true_right = sibling_right

        synthetic = pilot.simulate_counts(
            true_left,
            true_right,
            configuration["c_left"],
            configuration["c_right"],
            data,
            emission,
            rng,
        )

        fit = pilot.fit_pair(
            synthetic,
            emission,
            (
                f"sibling_t{task}_"
                f"{configuration['scenario']}_"
                f"{generating_model}"
            ),
        )

        result_rows.append(
            {
                "id": configuration["id"],
                "task": task,
                "scenario":
                    configuration["scenario"],
                "repetition":
                    configuration["repetition"],
                "pair_index":
                    configuration["pair_index"],
                "generating_model":
                    generating_model,
                "donor_left": left_sample,
                "donor_right": right_sample,
                "true_c_left":
                    configuration["c_left"],
                "true_c_right":
                    configuration["c_right"],
                "genotype_contrast":
                    float(
                        np.mean(
                            np.abs(
                                true_left
                                - true_right
                            )
                            / 2
                        )
                    ),
                "IBD0_fraction":
                    ibd["IBD0_fraction"],
                "IBD1_fraction":
                    ibd["IBD1_fraction"],
                "IBD2_fraction":
                    ibd["IBD2_fraction"],
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
            }
        )

    temporary = Path(
        str(output) + ".new"
    )

    with temporary.open("w") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=list(
                result_rows[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(result_rows)

    os.replace(
        temporary,
        output,
    )

    completion.write_text(
        "status\ttask\trows\tsites\n"
        f"PASS\t{task}\t"
        f"{len(result_rows)}\t"
        f"{len(data['pos'])}\n"
    )

    print(
        f"COMPLETE task={task} "
        f"scenario={configuration['scenario']} "
        f"pair={left_sample},{right_sample} "
        f"rows={len(result_rows)} "
        f"sites={len(data['pos'])}",
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--task",
        required=True,
        type=int,
    )

    arguments = parser.parse_args()
    main(arguments.task)
