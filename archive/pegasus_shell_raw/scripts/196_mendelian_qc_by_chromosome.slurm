#!/usr/bin/env bash
#SBATCH --job-name=gyn196sibmq
#SBATCH --partition=cpu
#SBATCH --array=0-19%20
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=04:00:00
#SBATCH --output=logs/196_mendelian_qc.%A_%a.out
#SBATCH --error=logs/196_mendelian_qc.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python
bcftools=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/bcftools

call_root="${project}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/variant_calls/by_chromosome"
sample_file="${project}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata/bam_set/full_sib_14_samples.samples"
output_root="${project}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/mendelian_qc/by_chromosome"

chromosome=$((SLURM_ARRAY_TASK_ID + 1))
started_epoch=$(date +%s)

cd "${project}"
mkdir -p "${output_root}"

printf 'STAGE_START\tstage=mendelian_qc\tchromosome=%s\tjob=%s\ttask=%s\thost=%s\ttime=%s\n' \
    "${chromosome}" \
    "${SLURM_ARRAY_JOB_ID}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/196_mendelian_qc_by_chromosome.py \
    --chromosome "${chromosome}" \
    --bcftools "${bcftools}" \
    --call-root "${call_root}" \
    --sample-file "${sample_file}" \
    --output-root "${output_root}" \
    --minimum-depth 5 \
    --minimum-pl-gap 20

ended_epoch=$(date +%s)

printf 'STAGE_END\tstage=mendelian_qc\tchromosome=%s\tstatus=0\telapsed_seconds=%s\thost=%s\ttime=%s\n' \
    "${chromosome}" \
    "$((ended_epoch - started_epoch))" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"
