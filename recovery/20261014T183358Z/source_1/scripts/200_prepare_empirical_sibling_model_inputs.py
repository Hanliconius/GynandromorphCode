#!/usr/bin/env python3

import csv
import importlib.util
import json
import os
import subprocess
from pathlib import Path

import numpy as np


ROOT = Path("analysis/genomewide_mixing_corrected_v1")
BASE = ROOT / "expanded_calibration"
FAMILY = ROOT / "full_sib_validation"
OUT = ROOT / "empirical_full_sib_calibration_v1"
INPUT = OUT / "inputs"

BCFTOOLS = Path(
    "/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/"
    ".local/share/mamba/envs/colias_w_popgen/bin/bcftools"
)

SAMPLE_TABLE = (
    FAMILY
    / "metadata/bam_set/full_sib_14_samples.tsv"
)


def sample_information():
    with SAMPLE_TABLE.open() as handle:
        rows = list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )

    if len(rows) != 14:
        raise RuntimeError(
            f"Expected 14 pedigree samples; found {len(rows)}"
        )

    if [row["role"] for row in rows[:2]] != [
        "F1_mother",
        "F1_father",
    ]:
        raise RuntimeError(
            "The first two pedigree samples are not "
            "the mother and father"
        )

    offspring = [
        row["sample"]
        for row in rows
        if row["role"] == "F2_offspring"
    ]

    if len(offspring) != 12:
        raise RuntimeError(
            f"Expected 12 offspring; found {len(offspring)}"
        )

    return offspring


def one_directory(root, chromosome):
    matches = sorted(
        root.glob(f"chr{chromosome}_*")
    )

    if len(matches) != 1:
        raise RuntimeError(
            f"Chromosome {chromosome}: expected one directory "
            f"under {root}; found {len(matches)}"
        )

    return matches[0]


def high_confidence_genotype(cell):
    fields = cell.split(":", 2)

    if len(fields) != 3:
        return -1

    gt, depth_text, pl_text = fields

    try:
        depth = int(depth_text)
        pls = [
            int(value)
            for value in pl_text.split(",")
        ]
    except ValueError:
        return -1

    alleles = (
        gt.replace("|", "/")
        .split("/")
    )

    if (
        len(alleles) != 2
        or any(
            allele not in {"0", "1"}
            for allele in alleles
        )
    ):
        return -1

    dosage = sum(
        int(allele)
        for allele in alleles
    )

    if depth < 5 or len(pls) < 3:
        return -1

    ordered = sorted(pls[:3])

    if ordered[1] - ordered[0] < 20:
        return -1

    if pls[dosage] != ordered[0]:
        return -1

    return dosage


def valid_existing(path, expected_sites):
    if not path.is_file() or path.stat().st_size == 0:
        return False

    try:
        with np.load(
            path,
            allow_pickle=False,
        ) as data:
            return (
                data["pos"].shape == (expected_sites,)
                and data["siblings"].shape
                == (expected_sites, 12)
                and data["offspring"].shape == (12,)
            )
    except Exception:
        return False


def main():
    offspring = sample_information()

    call_root = (
        FAMILY
        / "variant_calls/by_chromosome"
    )

    INPUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    audit_rows = []

    for chromosome in range(1, 21):
        source = (
            BASE
            / f"chr{chromosome}.npz"
        )

        output = (
            INPUT
            / f"chr{chromosome}.npz"
        )

        if not source.is_file():
            raise RuntimeError(
                f"Missing model input: {source}"
            )

        with np.load(
            source,
            allow_pickle=False,
        ) as data:
            arrays = {
                key: data[key]
                for key in data.files
            }

        sites = len(arrays["pos"])

        if valid_existing(output, sites):
            with np.load(
                output,
                allow_pickle=False,
            ) as existing:
                sibling_matrix = (
                    existing["siblings"].copy()
                )

            print(
                f"REUSE chromosome={chromosome} sites={sites}",
                flush=True,
            )
        else:
            call_directory = one_directory(
                call_root,
                chromosome,
            )

            bcf = (
                call_directory
                / "full_sib_14."
                "joint_at_existing_control_sites.bcf"
            )

            if (
                not bcf.is_file()
                or bcf.stat().st_size == 0
            ):
                raise RuntimeError(
                    f"Missing family BCF: {bcf}"
                )

            lookup = {
                (
                    int(position),
                    str(ref),
                    str(alt),
                ): index
                for index, (
                    position,
                    ref,
                    alt,
                ) in enumerate(
                    zip(
                        arrays["pos"],
                        arrays["ref"],
                        arrays["alt"],
                    )
                )
            }

            sibling_matrix = np.full(
                (sites, 12),
                -1,
                dtype=np.int8,
            )

            matched_sites = 0

            process = subprocess.Popen(
                [
                    str(BCFTOOLS),
                    "query",
                    "-f",
                    (
                        r"%POS\t%REF\t%ALT\t"
                        r"[%GT:%DP:%PL;]\n"
                    ),
                    str(bcf),
                ],
                stdout=subprocess.PIPE,
                text=True,
            )

            if process.stdout is None:
                raise RuntimeError(
                    "Could not read bcftools output"
                )

            for line_number, line in enumerate(
                process.stdout,
                start=1,
            ):
                fields = (
                    line.rstrip("\n")
                    .split("\t", 3)
                )

                if len(fields) != 4:
                    process.kill()

                    raise RuntimeError(
                        f"Chromosome {chromosome}, line "
                        f"{line_number}: expected four fields; "
                        f"found {len(fields)}"
                    )

                key = (
                    int(fields[0]),
                    fields[1],
                    fields[2],
                )

                model_index = lookup.get(key)

                if model_index is None:
                    continue

                cells = fields[3].split(";")

                if cells and cells[-1] == "":
                    cells.pop()

                if len(cells) != 14:
                    process.kill()

                    raise RuntimeError(
                        f"Chromosome {chromosome}, position "
                        f"{fields[0]}: expected 14 sample "
                        f"cells; found {len(cells)}"
                    )

                sibling_matrix[model_index] = [
                    high_confidence_genotype(cell)
                    for cell in cells[2:]
                ]

                matched_sites += 1

            query_status = process.wait()

            if query_status != 0:
                raise RuntimeError(
                    f"bcftools query failed with "
                    f"status {query_status}"
                )

            temporary = Path(
                str(output) + ".new.npz"
            )

            np.savez_compressed(
                temporary,
                **arrays,
                siblings=sibling_matrix,
                offspring=np.asarray(offspring),
            )

            os.replace(
                temporary,
                output,
            )

            print(
                f"COMPLETE chromosome={chromosome} "
                f"sites={sites} "
                f"matched_sites={matched_sites}",
                flush=True,
            )

        called = sibling_matrix >= 0
        per_sample_call_rate = called.mean(axis=0)

        audit_rows.append(
            {
                "chromosome": chromosome,
                "sites": sites,
                "mean_sibling_call_rate":
                    float(called.mean()),
                "minimum_sibling_call_rate":
                    float(per_sample_call_rate.min()),
                "maximum_sibling_call_rate":
                    float(per_sample_call_rate.max()),
                "output": str(output),
                "bytes": output.stat().st_size,
            }
        )

    total_sites = sum(
        row["sites"]
        for row in audit_rows
    )

    if total_sites != 1_000_000:
        raise RuntimeError(
            f"Expected 1,000,000 sites; "
            f"found {total_sites}"
        )

    audit_path = (
        OUT
        / "input_preparation.audit.tsv"
    )

    with audit_path.open("w") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=list(
                audit_rows[0].keys()
            ),
        )

        writer.writeheader()
        writer.writerows(audit_rows)

    overall_call_rate = float(
        np.average(
            [
                row["mean_sibling_call_rate"]
                for row in audit_rows
            ],
            weights=[
                row["sites"]
                for row in audit_rows
            ],
        )
    )

    diagnostic = {
        "status": "PASS",
        "chromosomes": 20,
        "offspring": 12,
        "sites": total_sites,
        "genotype_filters": {
            "minimum_DP": 5,
            "minimum_PL_gap": 20,
            "GT_must_match_best_PL": True,
        },
        "mean_sibling_call_rate":
            overall_call_rate,
        "missing_genotypes_in_simulations":
            "imputed using the existing control allele-frequency rule",
        "scientific_model_changed": False,
    }

    (
        OUT
        / "input_preparation.audit.json"
    ).write_text(
        json.dumps(
            diagnostic,
            indent=2,
        )
        + "\n"
    )

    (
        OUT
        / "input_preparation.complete.tsv"
    ).write_text(
        "status\tchromosomes\toffspring\tsites\n"
        "PASS\t20\t12\t1000000\n"
    )

    print(
        "ENDPOINT="
        "EMPIRICAL_SIBLING_MODEL_INPUTS_READY"
    )
    print(f"sites={total_sites}")
    print(
        f"mean_sibling_call_rate="
        f"{overall_call_rate:.6f}"
    )


if __name__ == "__main__":
    main()
