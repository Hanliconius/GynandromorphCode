#!/usr/bin/env bash
#SBATCH --job-name=gyn164fx
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=logs/164b_finalize_existing_downloads.%A_%a.out
#SBATCH --error=logs/164b_finalize_existing_downloads.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
root="${project}/analysis/genomewide_mixing_input"
manifest="${root}/metadata/download_manifest.tsv"
task_id=${SLURM_ARRAY_TASK_ID}

tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
pigz_bin="${tool_bin}/pigz"

cd "${project}"

if [[ ! -s "${manifest}" ]]
then
    printf 'ERROR: missing manifest: %s\n' "${manifest}" >&2
    false
fi

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

mapfile -t temp_dirs < <(
    find "${root}/tmp" \
        -mindepth 1 \
        -maxdepth 1 \
        -type d \
        -name "fasterq_${accession}.*.${task_id}" \
        -print
)

if [[ "${#temp_dirs[@]}" -ne 1 ]]
then
    printf 'ERROR: task %s (%s) has %s matching temporary directories; expected 1.\n' \
        "${task_id}" \
        "${accession}" \
        "${#temp_dirs[@]}" >&2

    printf '%s\n' "${temp_dirs[@]:-}" >&2
    false
fi

temp_dir=${temp_dirs[0]}

temp_r1="${temp_dir}/${accession}_1.fastq"
temp_r2="${temp_dir}/${accession}_2.fastq"
temp_r1_gz="${temp_r1}.gz"
temp_r2_gz="${temp_r2}.gz"

for mate in 1 2
do
    if [[ "${mate}" -eq 1 ]]
    then
        raw="${temp_r1}"
        compressed="${temp_r1_gz}"
    else
        raw="${temp_r2}"
        compressed="${temp_r2_gz}"
    fi

    if [[ -s "${compressed}" ]]
    then
        printf 'Validating existing compressed mate %s: %s\n' \
            "${mate}" \
            "${compressed}"

        "${pigz_bin}" -t "${compressed}"
    elif [[ -s "${raw}" ]]
    then
        printf 'Compressing raw mate %s: %s\n' \
            "${mate}" \
            "${raw}"

        "${pigz_bin}" -p 1 "${raw}"
        "${pigz_bin}" -t "${compressed}"
    else
        printf 'ERROR: neither raw nor compressed mate %s exists for %s in %s\n' \
            "${mate}" \
            "${accession}" \
            "${temp_dir}" >&2

        false
    fi
done

if [[ -e "${final_r1}" || -e "${final_r2}" ]]
then
    printf 'ERROR: final FASTQ path exists without a complete manifest; refusing overwrite.\n' >&2
    printf '%s\n%s\n' "${final_r1}" "${final_r2}" >&2
    false
fi

mkdir -p "${final_dir}" "${mapping_dir}"

mv -- "${temp_r1_gz}" "${final_r1}"
mv -- "${temp_r2_gz}" "${final_r2}"

"${pigz_bin}" -t "${final_r1}" "${final_r2}"

mapping_tmp="${mapping_file}.tmp.${SLURM_JOB_ID}"

printf '%s\n' \
    $'array_id\tsource_accession\tresolved_run\tlibrary_id\trole\tside\tleg\tmerge_group\tlayout\tfastq1\tfastq2' \
    > "${mapping_tmp}"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${task_id}" \
    "${accession}" \
    "${accession}" \
    "${library_id}" \
    "${role}" \
    "${side}" \
    "${leg}" \
    "${merge_group}" \
    PAIRED \
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

find "${temp_dir}" \
    -mindepth 1 \
    -type f \
    -delete

find "${temp_dir}" \
    -depth \
    -type d \
    -empty \
    -delete

printf 'FINALIZED\ttask=%s\taccession=%s\tlibrary=%s\n' \
    "${task_id}" \
    "${accession}" \
    "${library_id}"

printf 'FASTQ1\t%s\nFASTQ2\t%s\nMAPPING\t%s\n' \
    "${final_r1}" \
    "${final_r2}" \
    "${mapping_file}"
