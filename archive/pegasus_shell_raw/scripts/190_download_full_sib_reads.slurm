#!/usr/bin/env bash
#SBATCH --job-name=gyn190sibdl
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=logs/190_download_full_sib_reads.%A_%a.out
#SBATCH --error=logs/190_download_full_sib_reads.%A_%a.err

set -euo pipefail

project_root=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

manifest=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata/initial_full_sib_download_manifest.tsv

output_root=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/fastq

temporary_root=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/tmp/downloads

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python
curl_bin=/usr/bin/curl
gzip_bin=/usr/bin/gzip

cd "${project_root}"

task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}

IFS=$'\t' read -r \
    sample_accession \
    sample_role \
    progeny_number \
    run_accessions \
    mate1_urls \
    mate2_urls \
    expected_fastq_bytes \
    < <(
        "${python_bin}" - \
            "${manifest}" \
            "${task_id}" <<'PYEOF'
import csv
import sys

manifest = sys.argv[1]
task = int(sys.argv[2])

with open(
    manifest,
    newline="",
) as handle:
    rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

matching = [
    row
    for row in rows
    if int(row["array_id"]) == task
]

if len(matching) != 1:
    raise RuntimeError(
        "Expected one manifest row "
        f"for task {task}; found "
        f"{len(matching)}"
    )

row = matching[0]

values = [
    row["sample_accession"],
    row["role"],
    row["progeny_number"] or "NA",
    row["run_accessions"],
    row["mate1_urls"],
    row["mate2_urls"],
    row["fastq_bytes"],
]

if any(
    "\t" in value
    or "\n" in value
    for value in values
):
    raise RuntimeError(
        "Manifest value contains a "
        "tab or newline"
    )

print("\t".join(values))
PYEOF
    )

sample_dir=${output_root}/${sample_accession}

task_tmp=${temporary_root}/${sample_accession}.${SLURM_ARRAY_JOB_ID}.${task_id}

final_mate1=${sample_dir}/${sample_accession}_1.fastq.gz
final_mate2=${sample_dir}/${sample_accession}_2.fastq.gz
completion=${sample_dir}/download.complete.tsv

mkdir -p \
    "${sample_dir}" \
    "${task_tmp}" \
    logs

printf 'STAGE_START\tstage=download_full_sib_reads\tjob=%s\ttask=%s\tsample=%s\ttime=%s\n' \
    "${SLURM_ARRAY_JOB_ID}" \
    "${task_id}" \
    "${sample_accession}" \
    "$(date --iso-8601=seconds)"

existing_valid=0

if [[
    -s "${completion}" &&
    -s "${final_mate1}" &&
    -s "${final_mate2}"
]]
then
    if "${gzip_bin}" -t \
        "${final_mate1}" \
        "${final_mate2}"
    then
        existing_valid=1
    fi
fi

if [[ "${existing_valid}" -eq 1 ]]
then
    printf 'REUSE\tsample=%s\tmate1_bytes=%s\tmate2_bytes=%s\n' \
        "${sample_accession}" \
        "$(wc -c < "${final_mate1}")" \
        "$(wc -c < "${final_mate2}")"
else
    IFS=',' read -r -a mate1_array \
        <<< "${mate1_urls}"

    IFS=',' read -r -a mate2_array \
        <<< "${mate2_urls}"

    mate1_count=${#mate1_array[@]}
    mate2_count=${#mate2_array[@]}

    if [[ "${mate1_count}" -ne "${mate2_count}" || "${mate1_count}" -lt 1 ]]
    then
        printf 'ERROR: unequal or empty mate URL lists for %s\n' \
            "${sample_accession}" >&2
        false
    fi

    mate1_parts=()
    mate2_parts=()

    for index in "${!mate1_array[@]}"
    do
        part1=${task_tmp}/run_${index}_1.fastq.gz
        part2=${task_tmp}/run_${index}_2.fastq.gz

        printf 'DOWNLOAD\tsample=%s\trun_index=%s\tmate=1\turl=%s\n' \
            "${sample_accession}" \
            "${index}" \
            "${mate1_array[$index]}"

        "${curl_bin}" \
            --location \
            --fail \
            --show-error \
            --retry 5 \
            --retry-delay 5 \
            --continue-at - \
            --output "${part1}.part" \
            "${mate1_array[$index]}"

        mv \
            "${part1}.part" \
            "${part1}"

        printf 'DOWNLOAD\tsample=%s\trun_index=%s\tmate=2\turl=%s\n' \
            "${sample_accession}" \
            "${index}" \
            "${mate2_array[$index]}"

        "${curl_bin}" \
            --location \
            --fail \
            --show-error \
            --retry 5 \
            --retry-delay 5 \
            --continue-at - \
            --output "${part2}.part" \
            "${mate2_array[$index]}"

        mv \
            "${part2}.part" \
            "${part2}"

        "${gzip_bin}" -t \
            "${part1}" \
            "${part2}"

        mate1_parts+=(
            "${part1}"
        )

        mate2_parts+=(
            "${part2}"
        )
    done

    cat "${mate1_parts[@]}" \
        > "${final_mate1}.new"

    cat "${mate2_parts[@]}" \
        > "${final_mate2}.new"

    "${gzip_bin}" -t \
        "${final_mate1}.new" \
        "${final_mate2}.new"

    mv \
        "${final_mate1}.new" \
        "${final_mate1}"

    mv \
        "${final_mate2}.new" \
        "${final_mate2}"

    observed_fastq_bytes=$((
        $(wc -c < "${final_mate1}") +
        $(wc -c < "${final_mate2}")
    ))

    {
        printf 'sample_accession\trole\tprogeny_number\trun_accessions\tmate1_bytes\tmate2_bytes\ttotal_bytes\texpected_ENA_bytes\tcompleted_at\n'

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${sample_accession}" \
            "${sample_role}" \
            "${progeny_number}" \
            "${run_accessions}" \
            "$(wc -c < "${final_mate1}")" \
            "$(wc -c < "${final_mate2}")" \
            "${observed_fastq_bytes}" \
            "${expected_fastq_bytes}" \
            "$(date --iso-8601=seconds)"
    } > "${completion}.new"

    mv \
        "${completion}.new" \
        "${completion}"

    rm -f \
        "${task_tmp}"/run_*.fastq.gz

    printf 'COMPLETE\tsample=%s\truns=%s\ttotal_bytes=%s\n' \
        "${sample_accession}" \
        "${run_accessions}" \
        "${observed_fastq_bytes}"
fi

printf 'STAGE_END\tstage=download_full_sib_reads\tstatus=0\tsample=%s\ttime=%s\n' \
    "${sample_accession}" \
    "$(date --iso-8601=seconds)"
