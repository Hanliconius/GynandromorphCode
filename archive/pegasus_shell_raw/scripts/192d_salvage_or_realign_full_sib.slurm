#!/usr/bin/env bash
#SBATCH --job-name=gyn192sibsv
#SBATCH --partition=cpu
#SBATCH --array=0-13%14
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=logs/192d_salvage_or_realign_full_sib.%A_%a.out
#SBATCH --error=logs/192d_salvage_or_realign_full_sib.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
samtools=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/samtools
alignment_root="${project}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/alignment"
original_job=138256779
original_log="${project}/logs/192b_align_full_sib_pedigree_fixed_rg.${original_job}_${SLURM_ARRAY_TASK_ID}.out"
corrected_core="${project}/scripts/192c_align_full_sib_pedigree_corrected.slurm"

cd "${project}"
mkdir -p "${alignment_root}"

started_epoch=$(date +%s)

printf 'STAGE_START\tstage=salvage_or_realign_full_sib\tjob=%s\ttask=%s\thost=%s\ttime=%s\n' \
    "${SLURM_ARRAY_JOB_ID}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

sample=""

if [[ -s "${original_log}" ]]
then
    sample=$(
        awk -F '\t' '
            NR == 1 {
                for (field = 1; field <= NF; field++) {
                    if ($field ~ /^sample=/) {
                        sub(/^sample=/, "", $field)
                        print $field
                    }
                }
            }
        ' "${original_log}"
    )
fi

if [[ -z "${sample}" ]]
then
    printf 'ERROR: could not recover sample name from %s\n' \
        "${original_log}" >&2
    false
fi

bam="${alignment_root}/${sample}.dedup.bam"
bai="${bam}.bai"
flagstat="${alignment_root}/${sample}.dedup.flagstat.txt"
stats="${alignment_root}/${sample}.dedup.stats.txt"
completion="${alignment_root}/${sample}.complete.tsv"

printf 'TASK_SAMPLE\ttask=%s\tsample=%s\n' \
    "${SLURM_ARRAY_TASK_ID}" \
    "${sample}"

bam_valid=0

if [[ -s "${bam}" ]]
then
    if "${samtools}" quickcheck -v "${bam}"
    then
        observed_sample=$(
            "${samtools}" view -H "${bam}" \
                | awk -F '\t' '/^@RG/ {for (field=1; field<=NF; field++) {if ($field ~ /^SM:/) {sub(/^SM:/, "", $field); print $field}}}' \
                | sort -u
        )

        if [[ "${observed_sample}" == "${sample}" ]]
        then
            bam_valid=1
        else
            printf 'INVALID_EXISTING_BAM_RG\texpected=%s\tobserved=%s\n' \
                "${sample}" \
                "${observed_sample}"
        fi
    else
        printf 'INVALID_EXISTING_BAM_QUICKCHECK\t%s\n' \
            "${bam}"
    fi
fi

if [[ "${bam_valid}" -eq 1 ]]
then
    printf 'SALVAGING_VALID_BAM\t%s\n' "${bam}"

    if [[ ! -s "${bai}" ]]
    then
        "${samtools}" index \
            -@ "${SLURM_CPUS_PER_TASK}" \
            "${bam}"
    fi

    if [[ ! -s "${flagstat}" ]]
    then
        "${samtools}" flagstat \
            -@ "${SLURM_CPUS_PER_TASK}" \
            "${bam}" \
            > "${flagstat}.tmp"

        mv "${flagstat}.tmp" "${flagstat}"
    fi

    if [[ ! -s "${stats}" ]]
    then
        "${samtools}" stats \
            -@ "${SLURM_CPUS_PER_TASK}" \
            "${bam}" \
            > "${stats}.tmp"

        mv "${stats}.tmp" "${stats}"
    fi

    completion_tmp="${completion}.tmp.${SLURM_JOB_ID}"

    printf 'sample\tbam\tbam_bytes\tstatus\tcompleted_at\n' \
        > "${completion_tmp}"

    printf '%s\t%s\t%s\tPASS_SALVAGED\t%s\n' \
        "${sample}" \
        "${bam}" \
        "$(wc -c < "${bam}")" \
        "$(date --iso-8601=seconds)" \
        >> "${completion_tmp}"

    mv "${completion_tmp}" "${completion}"

    action=salvaged
else
    printf 'REALIGN_REQUIRED\ttask=%s\tsample=%s\n' \
        "${SLURM_ARRAY_TASK_ID}" \
        "${sample}"

    if [[ -f "${bam}" ]]
    then
        rm -f "${bam}"
    fi

    if [[ -f "${bai}" ]]
    then
        rm -f "${bai}"
    fi

    bash "${corrected_core}"

    action=realigned
fi

if [[ ! -s "${bam}" ]]
then
    printf 'ERROR: final BAM missing or empty: %s\n' \
        "${bam}" >&2
    false
fi

if [[ ! -s "${bai}" ]]
then
    printf 'ERROR: BAM index missing or empty: %s\n' \
        "${bai}" >&2
    false
fi

if [[ ! -s "${completion}" ]]
then
    printf 'ERROR: completion marker missing: %s\n' \
        "${completion}" >&2
    false
fi

"${samtools}" quickcheck -v "${bam}"

ended_epoch=$(date +%s)

printf 'SALVAGE_COMPLETE\ttask=%s\tsample=%s\taction=%s\tbam_bytes=%s\n' \
    "${SLURM_ARRAY_TASK_ID}" \
    "${sample}" \
    "${action}" \
    "$(wc -c < "${bam}")"

printf 'STAGE_END\tstage=salvage_or_realign_full_sib\tstatus=0\telapsed_seconds=%s\thost=%s\ttime=%s\n' \
    "$((ended_epoch - started_epoch))" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"
