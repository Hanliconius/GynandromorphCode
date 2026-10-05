#!/usr/bin/env bash
#SBATCH --job-name=gyn164cr
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/164c_repair_corrupt_fastq_array.%A_%a.out
#SBATCH --error=logs/164c_repair_corrupt_fastq_array.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
root="${project}/analysis/genomewide_mixing_input"
manifest="${root}/metadata/download_manifest.tsv"
task_id=${SLURM_ARRAY_TASK_ID}

tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
pigz_bin="${tool_bin}/pigz"
fasterq_bin="${tool_bin}/fasterq-dump"

cd "${project}"

field() {
    local column=$1

    awk -F '\t' \
        -v task="${task_id}" \
        -v column="${column}" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                index_by_name[$i] = i
            }
            next
        }

        $1 == task {
            print $index_by_name[column]
            found++
        }

        END {
            if (found != 1) {
                exit 2
            }
        }
    ' "${manifest}"
}

accession=$(field accession)
library_id=$(field library_id)
role=$(field role)
side=$(field side)
leg=$(field leg)
merge_group=$(field merge_group)

final_dir="${root}/fastq/${library_id}"
final_r1="${final_dir}/${library_id}_1.fastq.gz"
final_r2="${final_dir}/${library_id}_2.fastq.gz"

mapping_dir="${root}/metadata/resolved_runs"
mapping_file="${mapping_dir}/$(printf '%03d' "${task_id}").${library_id}.tsv"
complete_file="${final_dir}/download.complete"

if [[ \
    -s "${complete_file}" \
    && -s "${mapping_file}" \
    && -s "${final_r1}" \
    && -s "${final_r2}" \
]]
then
    "${pigz_bin}" -t "${final_r1}" "${final_r2}"

    printf 'Already complete and valid: task=%s library=%s\n' \
        "${task_id}" \
        "${library_id}"

    exit 0
fi

mapfile -t old_temp_dirs < <(
    find "${root}/tmp" \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        -name "fasterq_${accession}.*.${task_id}" \
        -print
)

if [[ "${#old_temp_dirs[@]}" -ne 1 ]]
then
    printf 'ERROR: task %s (%s) has %s old temporary directories; expected 1.\n' \
        "${task_id}" \
        "${accession}" \
        "${#old_temp_dirs[@]}" >&2

    printf '%s\n' "${old_temp_dirs[@]:-}" >&2
    false
fi

old_temp_dir=${old_temp_dirs[0]}
work_dir="${old_temp_dir}"

raw_r1="${work_dir}/${accession}_1.fastq"
raw_r2="${work_dir}/${accession}_2.fastq"
gz_r1="${raw_r1}.gz"
gz_r2="${raw_r2}.gz"

need_regeneration=0

inspect_mate() {
    local mate=$1
    local raw=$2
    local compressed=$3

    if [[ -s "${compressed}" ]]
    then
        if "${pigz_bin}" -t "${compressed}"
        then
            printf 'VALID_EXISTING_GZIP\tmate=%s\t%s\n' \
                "${mate}" \
                "${compressed}"
        elif [[ -s "${raw}" ]]
        then
            printf 'RECOMPRESSING_CORRUPT_GZIP\tmate=%s\t%s\n' \
                "${mate}" \
                "${compressed}"

            rm -f -- "${compressed}"

            "${pigz_bin}" -p 2 "${raw}"
            "${pigz_bin}" -t "${compressed}"
        else
            printf 'CORRUPT_GZIP_WITHOUT_RAW\tmate=%s\t%s\n' \
                "${mate}" \
                "${compressed}"

            need_regeneration=1
        fi
    elif [[ -s "${raw}" ]]
    then
        printf 'COMPRESSING_EXISTING_RAW\tmate=%s\t%s\n' \
            "${mate}" \
            "${raw}"

        "${pigz_bin}" -p 2 "${raw}"
        "${pigz_bin}" -t "${compressed}"
    else
        printf 'MATE_MISSING\tmate=%s\n' "${mate}"
        need_regeneration=1
    fi
}

inspect_mate 1 "${raw_r1}" "${gz_r1}"
inspect_mate 2 "${raw_r2}" "${gz_r2}"

if [[ "${need_regeneration}" -eq 1 ]]
then
    mapfile -t sra_files < <(
        find "${root}/sra/${accession}" \
            -type f \
            -name "${accession}.sra" \
            -print
    )

    if [[ "${#sra_files[@]}" -ne 1 ]]
    then
        printf 'ERROR: expected one local SRA file for %s; found %s.\n' \
            "${accession}" \
            "${#sra_files[@]}" >&2

        printf '%s\n' "${sra_files[@]:-}" >&2
        false
    fi

    recovery_dir="${root}/tmp/recovery_${accession}.${SLURM_JOB_ID}.${task_id}"

    mkdir -p "${recovery_dir}"

    printf 'Regenerating FASTQ from local SRA: %s\n' \
        "${sra_files[0]}"

    "${fasterq_bin}" \
        --split-files \
        --threads 2 \
        --temp "${recovery_dir}" \
        --outdir "${recovery_dir}" \
        "${sra_files[0]}"

    recovery_r1="${recovery_dir}/${accession}_1.fastq"
    recovery_r2="${recovery_dir}/${accession}_2.fastq"

    if [[ ! -s "${recovery_r1}" || ! -s "${recovery_r2}" ]]
    then
        printf 'ERROR: fasterq regeneration did not produce both mates.\n' >&2
        false
    fi

    "${pigz_bin}" -p 2 \
        "${recovery_r1}" \
        "${recovery_r2}"

    gz_r1="${recovery_r1}.gz"
    gz_r2="${recovery_r2}.gz"
    work_dir="${recovery_dir}"

    "${pigz_bin}" -t "${gz_r1}" "${gz_r2}"
fi

"${pigz_bin}" -t "${gz_r1}" "${gz_r2}"

if [[ -e "${final_r1}" || -e "${final_r2}" ]]
then
    printf 'ERROR: final FASTQ path exists without complete metadata; refusing overwrite.\n' >&2
    printf '%s\n%s\n' "${final_r1}" "${final_r2}" >&2
    false
fi

mkdir -p "${final_dir}" "${mapping_dir}"

mv -- "${gz_r1}" "${final_r1}"
mv -- "${gz_r2}" "${final_r2}"

"${pigz_bin}" -t "${final_r1}" "${final_r2}"

mapping_tmp="${mapping_file}.tmp.${SLURM_JOB_ID}"

printf '%s\n' \
    $'array_id\tsource_accession\tresolved_run\tlibrary_id\trole\tside\tleg\tmerge_group\tlayout\tfastq1\tfastq2' \
    > "${mapping_tmp}"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tPAIRED\t%s\t%s\n' \
    "${task_id}" \
    "${accession}" \
    "${accession}" \
    "${library_id}" \
    "${role}" \
    "${side}" \
    "${leg}" \
    "${merge_group}" \
    "${final_r1}" \
    "${final_r2}" \
    >> "${mapping_tmp}"

mv -- "${mapping_tmp}" "${mapping_file}"

printf 'task_id\t%s\nlibrary_id\t%s\nsource_accession\t%s\nresolved_run\t%s\nlayout\tPAIRED\ncompleted_at\t%s\n' \
    "${task_id}" \
    "${library_id}" \
    "${accession}" \
    "${accession}" \
    "$(date --iso-8601=seconds)" \
    > "${complete_file}"

for cleanup_dir in "${old_temp_dir}" "${work_dir}"
do
    if [[ \
        -d "${cleanup_dir}" \
        && "${cleanup_dir}" == "${root}/tmp/"* \
    ]]
    then
        find "${cleanup_dir}" \
            -mindepth 1 \
            -type f \
            -delete

        find "${cleanup_dir}" \
            -depth \
            -type d \
            -empty \
            -delete
    fi
done

printf 'REPAIRED\ttask=%s\taccession=%s\tlibrary=%s\n' \
    "${task_id}" \
    "${accession}" \
    "${library_id}"

printf 'FASTQ1\t%s\nFASTQ2\t%s\nMAPPING\t%s\n' \
    "${final_r1}" \
    "${final_r2}" \
    "${mapping_file}"
