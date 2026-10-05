#!/usr/bin/env bash
#SBATCH --job-name=gyn191sibqc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=logs/191_validate_full_sib_downloads.%j.out
#SBATCH --error=logs/191_validate_full_sib_downloads.%j.err

set -euo pipefail

project_root=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

manifest=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata/initial_full_sib_download_manifest.tsv

fastq_root=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/fastq

audit=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata/full_sib_download_validation.tsv

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

cd "${project_root}"

mkdir -p \
    "$(dirname "${audit}")" \
    logs

printf 'STAGE_START\tstage=validate_full_sib_downloads\tjob=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(date --iso-8601=seconds)"

"${python_bin}" - \
    "${manifest}" \
    "${fastq_root}" \
    "${audit}" <<'PYEOF'
import csv
import os
import sys
from pathlib import Path


manifest_path = Path(sys.argv[1])
fastq_root = Path(sys.argv[2])
audit_path = Path(sys.argv[3])

with manifest_path.open(
    newline="",
) as handle:
    rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

results = []
failed = 0

for row in rows:
    sample = row["sample_accession"]
    directory = fastq_root / sample

    mate1 = (
        directory
        / f"{sample}_1.fastq.gz"
    )

    mate2 = (
        directory
        / f"{sample}_2.fastq.gz"
    )

    completion = (
        directory
        / "download.complete.tsv"
    )

    valid = all(
        path.is_file()
        and path.stat().st_size > 0
        for path
        in (
            mate1,
            mate2,
            completion,
        )
    )

    failed += int(not valid)

    results.append(
        {
            "array_id":
                row["array_id"],
            "sample_accession":
                sample,
            "role":
                row["role"],
            "progeny_number":
                (
                    row["progeny_number"]
                    or "NA"
                ),
            "status":
                (
                    "PASS"
                    if valid
                    else "MISSING_OR_EMPTY"
                ),
            "mate1_bytes":
                (
                    mate1.stat().st_size
                    if mate1.is_file()
                    else 0
                ),
            "mate2_bytes":
                (
                    mate2.stat().st_size
                    if mate2.is_file()
                    else 0
                ),
        }
    )

fields = list(results[0])

temporary = Path(
    str(audit_path) + ".new"
)

with temporary.open(
    "w",
    newline="",
) as handle:
    writer = csv.DictWriter(
        handle,
        delimiter="\t",
        fieldnames=fields,
    )
    writer.writeheader()
    writer.writerows(results)

os.replace(
    temporary,
    audit_path,
)

with audit_path.open() as handle:
    print(
        handle.read(),
        end="",
    )

print(
    f"samples={len(results)}"
)
print(
    f"failed_samples={failed}"
)

if failed:
    raise RuntimeError(
        f"{failed} full-sib FASTQ "
        "samples failed validation"
    )
PYEOF

printf 'ENDPOINT=FULL_SIB_FASTQS_READY\n'

printf 'STAGE_END\tstage=validate_full_sib_downloads\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
