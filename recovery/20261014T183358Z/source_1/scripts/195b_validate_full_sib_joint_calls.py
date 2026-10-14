#!/usr/bin/env python3

import csv
import gzip
import json
import os
import re
import subprocess
import time
from pathlib import Path


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/"
    "GynandroCalling"
)

SOURCE_ROOT = (
    PROJECT
    / "analysis/genomewide_mixing_input/"
    "variant_calls/joint_24/by_chromosome"
)

CALL_ROOT = (
    PROJECT
    / "analysis/genomewide_mixing_corrected_v1/"
    "full_sib_validation/variant_calls/by_chromosome"
)

METADATA_ROOT = (
    PROJECT
    / "analysis/genomewide_mixing_corrected_v1/"
    "full_sib_validation/metadata"
)

EXPECTED_SAMPLES = (
    METADATA_ROOT
    / "bam_set/full_sib_14_samples.samples"
)

BCFTOOLS = (
    "/gpfs/automountdir/gpfs/homes/CCAS/home/"
    "jjhanly/.local/share/mamba/envs/"
    "colias_w_popgen/bin/bcftools"
)


def count_target_sites(path):
    count = 0

    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                count += 1

    return count


def read_completion(path):
    with path.open(newline="") as handle:
        rows = list(
            csv.DictReader(handle, delimiter="\t")
        )

    if len(rows) != 1:
        raise RuntimeError(
            f"Malformed completion marker: {path}"
        )

    return rows[0]


def write_tsv(path, fieldnames, rows):
    temporary = Path(str(path) + ".tmp")

    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    os.replace(temporary, path)


def main():
    started = time.time()

    expected_samples = [
        line.strip()
        for line in EXPECTED_SAMPLES.read_text().splitlines()
        if line.strip()
    ]

    if len(expected_samples) != 14:
        raise RuntimeError(
            f"Expected 14 samples; found {len(expected_samples)}"
        )

    manifest_rows = []

    for chromosome in range(1, 22):
        call_directories = sorted(
            CALL_ROOT.glob(f"chr{chromosome}_*")
        )
        source_directories = sorted(
            SOURCE_ROOT.glob(f"chr{chromosome}_*")
        )

        if len(call_directories) != 1:
            raise RuntimeError(
                f"Chromosome {chromosome}: expected one "
                f"call directory; found {call_directories}"
            )

        if len(source_directories) != 1:
            raise RuntimeError(
                f"Chromosome {chromosome}: expected one "
                f"source directory; found {source_directories}"
            )

        call_directory = call_directories[0]
        source_directory = source_directories[0]
        contig = re.sub(
            rf"^chr{chromosome}_",
            "",
            call_directory.name,
        )

        bcf = (
            call_directory
            / "full_sib_14.joint_at_existing_control_sites.bcf"
        )
        index = Path(str(bcf) + ".csi")
        completion = call_directory / "call.complete.tsv"
        alleles = (
            source_directory
            / "controls.discovery.alleles.tsv.gz"
        )

        for required in (
            bcf,
            index,
            completion,
            alleles,
        ):
            if (
                not required.is_file()
                or required.stat().st_size == 0
            ):
                raise RuntimeError(
                    f"Missing or empty required file: {required}"
                )

        samples = subprocess.run(
            [
                BCFTOOLS,
                "query",
                "-l",
                str(bcf),
            ],
            check=True,
            text=True,
            capture_output=True,
        ).stdout.splitlines()

        if samples != expected_samples:
            raise RuntimeError(
                f"Chromosome {chromosome}: sample order mismatch"
            )

        emitted_sites = int(
            subprocess.run(
                [
                    BCFTOOLS,
                    "index",
                    "-n",
                    str(bcf),
                ],
                check=True,
                text=True,
                capture_output=True,
            ).stdout.strip()
        )

        target_sites = count_target_sites(alleles)

        if (
            emitted_sites <= 0
            or emitted_sites > target_sites
        ):
            raise RuntimeError(
                f"Chromosome {chromosome}: invalid site counts "
                f"emitted={emitted_sites}, target={target_sites}"
            )

        completion_row = read_completion(completion)

        if int(completion_row["extra_sites"]) != 0:
            raise RuntimeError(
                f"Chromosome {chromosome}: non-target sites reported"
            )

        if (
            int(completion_row["emitted_sites"])
            != emitted_sites
        ):
            raise RuntimeError(
                f"Chromosome {chromosome}: completion count mismatch"
            )

        manifest_rows.append(
            {
                "chromosome": chromosome,
                "contig": contig,
                "samples": len(samples),
                "target_sites": target_sites,
                "emitted_sites": emitted_sites,
                "retained_fraction":
                    emitted_sites / target_sites,
                "missing_target_sites":
                    target_sites - emitted_sites,
                "extra_sites": 0,
                "bcf": bcf,
                "bcf_bytes": bcf.stat().st_size,
                "index": index,
                "completion": completion,
            }
        )

        print(
            "PASS"
            f"\tchromosome={chromosome}"
            f"\ttarget={target_sites}"
            f"\temitted={emitted_sites}"
            f"\tretained={emitted_sites / target_sites:.6f}"
        )

    manifest_path = (
        METADATA_ROOT
        / "full_sib_21chromosome_call_manifest.tsv"
    )

    write_tsv(
        manifest_path,
        [
            "chromosome",
            "contig",
            "samples",
            "target_sites",
            "emitted_sites",
            "retained_fraction",
            "missing_target_sites",
            "extra_sites",
            "bcf",
            "bcf_bytes",
            "index",
            "completion",
        ],
        manifest_rows,
    )

    total_target = sum(
        row["target_sites"]
        for row in manifest_rows
    )
    total_emitted = sum(
        row["emitted_sites"]
        for row in manifest_rows
    )

    audit = {
        "chromosomes": len(manifest_rows),
        "samples": 14,
        "total_target_sites": total_target,
        "total_emitted_sites": total_emitted,
        "overall_retained_fraction":
            total_emitted / total_target,
        "total_missing_target_sites":
            total_target - total_emitted,
        "total_extra_sites": 0,
        "site_relationship":
            "family calls are a validated subset of existing control sites",
        "reason_missing_sites_are_allowed":
            "bcftools -T restricts candidate positions but does not "
            "force output where all family samples lack usable pileup data",
        "scientific_filters_changed": False,
        "elapsed_seconds":
            round(time.time() - started, 3),
        "pipeline_status": "PASS",
    }

    audit_path = (
        METADATA_ROOT
        / "full_sib_joint_calls.audit.json"
    )

    temporary_audit = Path(str(audit_path) + ".tmp")
    temporary_audit.write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    os.replace(temporary_audit, audit_path)

    print("\nFULL-SIB JOINT-CALL VALIDATION PASSED")
    print(json.dumps(audit, indent=2))
    print("ENDPOINT=FULL_SIB_JOINT_CALLS_READY")


if __name__ == "__main__":
    main()
