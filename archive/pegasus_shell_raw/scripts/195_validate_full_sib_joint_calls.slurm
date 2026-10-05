#!/usr/bin/env bash
#SBATCH --job-name=gyn195sibqc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/195_validate_full_sib_joint_calls.%j.out
#SBATCH --error=logs/195_validate_full_sib_joint_calls.%j.err

set -euo pipefail

project_root=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

bcftools=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/bcftools

source_root=${project_root}/analysis/genomewide_mixing_input/variant_calls/joint_24/by_chromosome

call_root=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/variant_calls/by_chromosome

set_dir=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata/bam_set

metadata_dir=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata

cd "${project_root}"

mkdir -p \
    "${metadata_dir}" \
    logs

printf 'STAGE_START\tstage=validate_full_sib_joint_calls\tjob=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(date --iso-8601=seconds)"

"${python_bin}" - \
    "${source_root}" \
    "${call_root}" \
    "${set_dir}" \
    "${metadata_dir}" \
    "${bcftools}" <<'PYEOF'
import csv
import gzip
import os
import subprocess
import sys
from pathlib import Path


source_root = Path(sys.argv[1])
call_root = Path(sys.argv[2])
set_dir = Path(sys.argv[3])
metadata_dir = Path(sys.argv[4])
bcftools = sys.argv[5]

expected_samples_path = (
    set_dir
    / "full_sib_14_samples.samples"
)

expected_samples = (
    expected_samples_path
    .read_text()
    .splitlines()
)

if (
    len(expected_samples) != 14
    or len(set(expected_samples)) != 14
):
    raise RuntimeError(
        "Expected 14 unique samples in "
        "the BAM-set sample file"
    )

results = []

for chromosome in range(1, 22):
    source_matches = sorted(
        path
        for path in source_root.glob(
            f"chr{chromosome}_*"
        )
        if path.is_dir()
    )

    call_matches = sorted(
        path
        for path in call_root.glob(
            f"chr{chromosome}_*"
        )
        if path.is_dir()
    )

    if (
        len(source_matches) != 1
        or len(call_matches) != 1
    ):
        raise RuntimeError(
            f"Chromosome {chromosome}: "
            "source directories="
            f"{len(source_matches)}, "
            "call directories="
            f"{len(call_matches)}"
        )

    source_directory = (
        source_matches[0]
    )
    call_directory = (
        call_matches[0]
    )

    contig = (
        source_directory.name
        .removeprefix(
            f"chr{chromosome}_"
        )
    )

    expected_call_name = (
        f"chr{chromosome}_{contig}"
    )

    if (
        call_directory.name
        != expected_call_name
    ):
        raise RuntimeError(
            f"Chromosome {chromosome}: "
            "call-directory contig mismatch"
        )

    alleles = (
        source_directory
        / "controls.discovery.alleles.tsv.gz"
    )

    bcf = (
        call_directory
        / (
            "full_sib_14."
            "joint_at_existing_control_sites.bcf"
        )
    )

    index = Path(
        str(bcf) + ".csi"
    )

    completion = (
        call_directory
        / "call.complete.tsv"
    )

    for path in (
        alleles,
        bcf,
        index,
        completion,
    ):
        if (
            not path.is_file()
            or path.stat().st_size == 0
        ):
            raise RuntimeError(
                "Missing or empty required "
                f"file: {path}"
            )

    samples = subprocess.run(
        [
            bcftools,
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
            f"Chromosome {chromosome}: "
            "sample order mismatch"
        )

    observed_sites = int(
        subprocess.run(
            [
                bcftools,
                "index",
                "-n",
                str(bcf),
            ],
            check=True,
            text=True,
            capture_output=True,
        ).stdout.strip()
    )

    with gzip.open(
        alleles,
        "rt",
    ) as handle:
        expected_sites = sum(
            1
            for line in handle
            if (
                line
                and not line.startswith("#")
            )
        )

    if observed_sites != expected_sites:
        raise RuntimeError(
            f"Chromosome {chromosome}: "
            f"expected {expected_sites} "
            f"sites; found {observed_sites}"
        )

    results.append(
        {
            "chromosome":
                chromosome,
            "contig":
                contig,
            "samples":
                len(samples),
            "sites":
                observed_sites,
            "bcf":
                str(bcf),
            "index":
                str(index),
            "bcf_bytes":
                bcf.stat().st_size,
        }
    )

manifest = (
    metadata_dir
    / (
        "full_sib_21chromosome_"
        "call_manifest.tsv"
    )
)

temporary = Path(
    str(manifest) + ".new"
)

with temporary.open(
    "w",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=list(
            results[0]
        ),
    )
    writer.writeheader()
    writer.writerows(results)

os.replace(
    temporary,
    manifest,
)

audit = (
    metadata_dir
    / "full_sib_joint_calls.audit.tsv"
)

audit_rows = [
    (
        "chromosomes",
        len(results),
    ),
    (
        "autosomes",
        sum(
            row["chromosome"] <= 20
            for row in results
        ),
    ),
    (
        "Z_chromosomes",
        sum(
            row["chromosome"] == 21
            for row in results
        ),
    ),
    (
        "samples_per_chromosome",
        len(expected_samples),
    ),
    (
        "total_sites",
        sum(
            row["sites"]
            for row in results
        ),
    ),
    (
        "site_ascertainment",
        "existing_18-control discovery alleles",
    ),
    (
        "family_calling_HWE_prior",
        "disabled with bcftools call -G -",
    ),
    (
        "concatenated_copy_created",
        False,
    ),
    (
        "next_stage",
        (
            "Mendelian QC and empirical "
            "IBD/state-space audit"
        ),
    ),
]

temporary = Path(
    str(audit) + ".new"
)

with temporary.open("w") as handle:
    handle.write(
        "item\tvalue\n"
    )

    for item, value in audit_rows:
        handle.write(
            f"{item}\t{value}\n"
        )

os.replace(
    temporary,
    audit,
)

print(
    "FULL-SIB JOINT-CALL "
    "VALIDATION PASSED"
)

print(
    audit.read_text(),
    end="",
)
PYEOF

printf 'ENDPOINT=FULL_SIB_JOINT_CALLS_READY_FOR_MENDELIAN_QC\n'

printf 'STAGE_END\tstage=validate_full_sib_joint_calls\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
