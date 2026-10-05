#!/usr/bin/env bash
#SBATCH --job-name=gyn167rt
#SBATCH --partition=cpu
#SBATCH --array=0-21
#SBATCH --cpus-per-task=3
#SBATCH --mem=20G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/167b_align_genomewide_libraries_safe_tsv.%A_%a.out
#SBATCH --error=logs/167b_align_genomewide_libraries_safe_tsv.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
manifest="${project}/analysis/genomewide_mixing_input/metadata/resolved_download_manifest.tsv"
reference_path_file="${project}/analysis/genomewide_mixing_input/metadata/reference_fasta.path"
output_dir="${project}/analysis/genomewide_mixing_input/alignment/libraries"
tmp_root="${project}/analysis/genomewide_mixing_input/tmp/alignment"

bwa="${tool_bin}/bwa"
samtools="${tool_bin}/samtools"

mkdir -p \
    "${project}/logs" \
    "${output_dir}" \
    "${tmp_root}"

for required in \
    "${manifest}" \
    "${reference_path_file}" \
    "${bwa}" \
    "${samtools}"
do
    if [[ ! -s "${required}" ]]
    then
        printf 'ERROR: required input is missing or empty: %s\n' \
            "${required}" >&2
        false
    fi
done

reference=$(head -n 1 "${reference_path_file}")

if [[
    "${reference}" != *'/GCF_056956725.1/'* ||
    ! -s "${reference}" ||
    ! -s "${reference}.fai"
]]
then
    printf 'ERROR: active reference is not indexed GCF_056956725.1: %s\n' \
        "${reference}" >&2
    false
fi

mapfile -t selected_fields < <(
    awk \
        -F '\t' \
        -v id="${SLURM_ARRAY_TASK_ID}" \
        '
            NR > 1 && $1 == id {
                print $1
                print $2
                print $3
                print $4
                print $5
                print $9
                print $10
                print $11
                found++
            }

            END {
                if (found != 1) {
                    exit 1
                }
            }
        ' \
        "${manifest}"
)

if [[ "${#selected_fields[@]}" -ne 8 ]]
then
    printf 'ERROR: expected eight manifest values for task %s; found %s.
' \
        "${SLURM_ARRAY_TASK_ID}" \
        "${#selected_fields[@]}" >&2
    false
fi

array_id="${selected_fields[0]}"
source_accession="${selected_fields[1]}"
resolved_run="${selected_fields[2]}"
library_id="${selected_fields[3]}"
role="${selected_fields[4]}"
layout="${selected_fields[5]}"
fastq1="${selected_fields[6]}"
fastq2="${selected_fields[7]}"

if [[
    "${layout}" != "PAIRED" ||
    ! -s "${fastq1}" ||
    ! -s "${fastq2}"
]]
then
    printf 'ERROR: invalid paired FASTQ inputs for task %s (%s)\n' \
        "${array_id}" \
        "${library_id}" >&2
    false
fi

final_bam="${output_dir}/${library_id}.dedup.bam"
final_bai="${final_bam}.bai"
flagstat="${output_dir}/${library_id}.dedup.flagstat.txt"
stats="${output_dir}/${library_id}.dedup.stats.txt"
markdup_stats="${output_dir}/${library_id}.markdup.stats.txt"
done_file="${output_dir}/${library_id}.complete.tsv"

if [[
    -s "${final_bam}" &&
    -s "${final_bai}" &&
    -s "${done_file}"
]] && "${samtools}" quickcheck -q "${final_bam}"
then
    printf 'Already complete and valid: task=%s library=%s\n' \
        "${array_id}" \
        "${library_id}"
else
    work_dir="${tmp_root}/${library_id}.${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}"

    mkdir -p "${work_dir}"

    name_bam="${work_dir}/${library_id}.name.bam"
    fixmate_bam="${work_dir}/${library_id}.fixmate.bam"
    position_bam="${work_dir}/${library_id}.position.bam"
    new_bam="${work_dir}/${library_id}.dedup.bam"

    rg="@RG\\tID:${library_id}\\tSM:${library_id}\\tLB:${library_id}\\tPL:ILLUMINA\\tPU:${resolved_run}"

    "${bwa}" mem \
        -M \
        -t 2 \
        -R "${rg}" \
        "${reference}" \
        "${fastq1}" \
        "${fastq2}" \
        | "${samtools}" sort \
            -n \
            -@ 1 \
            -m 3G \
            -o "${name_bam}" \
            -

    "${samtools}" fixmate \
        -m \
        -@ 1 \
        "${name_bam}" \
        "${fixmate_bam}"

    "${samtools}" sort \
        -@ 2 \
        -m 3G \
        -o "${position_bam}" \
        "${fixmate_bam}"

    "${samtools}" markdup \
        -r \
        -@ 2 \
        -s \
        -f "${markdup_stats}.new" \
        "${position_bam}" \
        "${new_bam}"

    "${samtools}" index \
        -@ 2 \
        "${new_bam}"

    "${samtools}" quickcheck \
        -v \
        "${new_bam}"

    "${samtools}" flagstat \
        -@ 2 \
        "${new_bam}" \
        > "${flagstat}.new"

    "${samtools}" stats \
        -@ 2 \
        "${new_bam}" \
        > "${stats}.new"

    mv "${new_bam}" "${final_bam}"
    mv "${new_bam}.bai" "${final_bai}"
    mv "${flagstat}.new" "${flagstat}"
    mv "${stats}.new" "${stats}"
    mv "${markdup_stats}.new" "${markdup_stats}"

    printf 'array_id\tlibrary_id\trole\tsource_accession\tresolved_run\tbam\treference\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${array_id}" \
        "${library_id}" \
        "${role}" \
        "${source_accession}" \
        "${resolved_run}" \
        "${final_bam}" \
        "${reference}" \
        > "${done_file}.new"

    mv "${done_file}.new" "${done_file}"

    rm -f \
        "${name_bam}" \
        "${fixmate_bam}" \
        "${position_bam}"

    rmdir "${work_dir}" 2>/dev/null || true
fi

printf 'ALIGNMENT_COMPLETE\ttask=%s\tlibrary=%s\tbam=%s\n' \
    "${array_id}" \
    "${library_id}" \
    "${final_bam}"
