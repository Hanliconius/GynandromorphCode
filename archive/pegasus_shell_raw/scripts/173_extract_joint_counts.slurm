#!/usr/bin/env bash
#SBATCH --job-name=gyn173ex
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=06:00:00
#SBATCH --output=logs/173_extract_joint_counts.%A_%a.out
#SBATCH --error=logs/173_extract_joint_counts.%A_%a.err

set -euo pipefail

stage_name=extract_joint_counts
start_seconds=${SECONDS}
printf 'STAGE_START\tstage=%s\tjob=%s\ttask=%s\thost=%s\ttime=%s\n' \
    "${stage_name}" "${SLURM_JOB_ID:-NA}" "${SLURM_ARRAY_TASK_ID:-NA}" \
    "$(hostname)" "$(date --iso-8601=seconds)"
trap 'status=$?; printf "STAGE_END\tstage=%s\tstatus=%s\telapsed_seconds=%s\thost=%s\ttime=%s\n" "${stage_name}" "${status}" "$((SECONDS - start_seconds))" "$(hostname)" "$(date --iso-8601=seconds)"' EXIT

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
environment_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
bcftools="${environment_bin}/bcftools"
manifest="${project}/analysis/genomewide_mixing_input/metadata/chromosome_calling_manifest.tsv"
output_dir="${project}/analysis/genomewide_mixing/inputs"

mkdir -p "${output_dir}"

task_id="${SLURM_ARRAY_TASK_ID}"

row=$(
    awk -F '\t' -v task_id="${task_id}" \
        'NR > 1 && $1 == task_id {print; found++}
         END {if (found != 1) exit 1}' \
        "${manifest}"
)

IFS=$'\t' read -r \
    array_id molecule contig length chromosome_type \
    <<< "${row}"

chrom_dir="${project}/analysis/genomewide_mixing_input/variant_calls/joint_24/by_chromosome/chr${molecule}_${contig}"
bcf="${chrom_dir}/all_24_samples.joint_at_control_sites.bcf"
completion="${chrom_dir}/calling.complete.tsv"
output="${output_dir}/chr${molecule}_${contig}.counts.tsv.gz"

if [[ ! -s "${bcf}" || ! -s "${bcf}.csi" || ! -s "${completion}" ]]
then
    printf 'ERROR: chromosome %s lacks a validated final joint BCF.\n' \
        "${molecule}" >&2
    false
fi

if [[ -s "${output}" && "${output}" -nt "${bcf}" ]] \
    && gzip -t "${output}" 2>/dev/null
then
    printf 'ALREADY_COMPLETE\tchromosome=%s\tfile=%s\n' \
        "${molecule}" "${output}"
else
    sample_header=$(
        "${bcftools}" query -l "${bcf}" \
            | paste -sd $'\t' -
    )

    sample_count=$(
        "${bcftools}" query -l "${bcf}" \
            | awk 'END {print NR + 0}'
    )

    if [[ "${sample_count}" -ne 24 ]]
    then
        printf 'ERROR: chromosome %s contains %s samples, not 24.\n' \
            "${molecule}" "${sample_count}" >&2
        false
    fi

    {
        printf '#CHROM\tPOS\tREF\tALT\t%s\n' "${sample_header}"

        "${bcftools}" query \
            -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT:%AD:%DP]\n' \
            "${bcf}"
    } \
        | gzip -c \
        > "${output}.new"

    gzip -t "${output}.new"

    row_count=$(
        gzip -dc "${output}.new" \
            | awk 'NR > 1 {count++} END {print count + 0}'
    )

    expected_rows=$(
        "${bcftools}" index -n "${bcf}"
    )

    if [[ "${row_count}" -ne "${expected_rows}" ]]
    then
        printf 'ERROR: chromosome %s extraction has %s rows; expected %s.\n' \
            "${molecule}" "${row_count}" "${expected_rows}" >&2
        false
    fi

    mv "${output}.new" "${output}"

    printf 'COUNT_EXTRACTION_COMPLETE\tchromosome=%s\tcontig=%s\trows=%s\tfile=%s\n' \
        "${molecule}" "${contig}" "${row_count}" "${output}"
fi
