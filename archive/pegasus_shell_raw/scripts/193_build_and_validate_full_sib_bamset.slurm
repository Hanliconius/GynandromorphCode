#!/usr/bin/env bash
#SBATCH --job-name=gyn193sibbs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/193_build_and_validate_full_sib_bamset.%j.out
#SBATCH --error=logs/193_build_and_validate_full_sib_bamset.%j.err

set -euo pipefail

project_root=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

samtools=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/samtools

manifest=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata/initial_full_sib_download_manifest.tsv

alignment_root=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/alignment

set_dir=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata/bam_set

cd "${project_root}"

mkdir -p \
    "${set_dir}" \
    logs

printf 'STAGE_START\tstage=build_full_sib_bamset\tjob=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(date --iso-8601=seconds)"

"${python_bin}" - \
    "${manifest}" \
    "${alignment_root}" \
    "${set_dir}" \
    "${samtools}" <<'PYEOF'
import csv
import os
import subprocess
import sys
from pathlib import Path


manifest_path = Path(sys.argv[1])
alignment_root = Path(sys.argv[2])
set_dir = Path(sys.argv[3])
samtools = sys.argv[4]

with manifest_path.open(
    newline="",
) as handle:
    rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

if len(rows) != 14:
    raise RuntimeError(
        "Expected 14 pedigree samples; "
        f"found {len(rows)}"
    )

samples = [
    row["sample_accession"]
    for row in rows
]

if len(set(samples)) != 14:
    raise RuntimeError(
        "Pedigree sample accessions "
        "are not unique"
    )

output_rows = []

for row in rows:
    sample = row["sample_accession"]

    bam = (
        alignment_root
        / f"{sample}.dedup.bam"
    )

    bai = Path(
        str(bam) + ".bai"
    )

    completion = (
        alignment_root
        / f"{sample}.complete.tsv"
    )

    for path in (
        bam,
        bai,
        completion,
    ):
        if (
            not path.is_file()
            or path.stat().st_size == 0
        ):
            raise RuntimeError(
                "Missing or empty "
                f"required file: {path}"
            )

    subprocess.run(
        [
            samtools,
            "quickcheck",
            "-v",
            str(bam),
        ],
        check=True,
    )

    header = subprocess.run(
        [
            samtools,
            "view",
            "-H",
            str(bam),
        ],
        check=True,
        text=True,
        capture_output=True,
    ).stdout

    header_samples = set()

    for line in header.splitlines():
        if not line.startswith("@RG\t"):
            continue

        for field in line.split("\t"):
            if field.startswith("SM:"):
                header_samples.add(
                    field[3:]
                )

    if header_samples != {sample}:
        raise RuntimeError(
            f"BAM {bam} has read-group "
            f"samples {sorted(header_samples)}; "
            f"expected {sample}"
        )

    output_rows.append(
        {
            "sample": sample,
            "role": row["role"],
            "progeny_number": (
                row["progeny_number"]
                or "NA"
            ),
            "bam": str(bam),
            "bam_bytes":
                bam.stat().st_size,
        }
    )

table_path = (
    set_dir
    / "full_sib_14_samples.tsv"
)

bamlist_path = (
    set_dir
    / "full_sib_14_samples.bamlist"
)

samples_path = (
    set_dir
    / "full_sib_14_samples.samples"
)

audit_path = (
    set_dir
    / "full_sib_14_samples.audit.tsv"
)

temporary = Path(
    str(table_path) + ".new"
)

with temporary.open(
    "w",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=list(
            output_rows[0]
        ),
    )
    writer.writeheader()
    writer.writerows(output_rows)

os.replace(
    temporary,
    table_path,
)

temporary = Path(
    str(bamlist_path) + ".new"
)

temporary.write_text(
    "".join(
        f"{row['bam']}\n"
        for row in output_rows
    )
)

os.replace(
    temporary,
    bamlist_path,
)

temporary = Path(
    str(samples_path) + ".new"
)

temporary.write_text(
    "".join(
        f"{row['sample']}\n"
        for row in output_rows
    )
)

os.replace(
    temporary,
    samples_path,
)

role_counts = {}

for row in output_rows:
    role = row["role"]

    role_counts[role] = (
        role_counts.get(role, 0)
        + 1
    )

offspring = role_counts.get(
    "F2_offspring",
    0,
)

audit_rows = [
    (
        "samples",
        len(output_rows),
    ),
    (
        "F1_mother",
        role_counts.get(
            "F1_mother",
            0,
        ),
    ),
    (
        "F1_father",
        role_counts.get(
            "F1_father",
            0,
        ),
    ),
    (
        "F2_offspring",
        offspring,
    ),
    (
        "possible_sibling_pairs",
        offspring
        * (offspring - 1)
        // 2,
    ),
    (
        "bamlist",
        str(bamlist_path),
    ),
]

temporary = Path(
    str(audit_path) + ".new"
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
    audit_path,
)

print(
    "FULL-SIB BAM SET COMPLETE"
)

print(
    audit_path.read_text(),
    end="",
)
PYEOF

printf 'ENDPOINT=FULL_SIB_BAMSET_READY\n'

printf 'STAGE_END\tstage=build_full_sib_bamset\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
