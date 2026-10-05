#!/usr/bin/env bash
#SBATCH --job-name=gyn198sibib
#SBATCH --partition=nano
#SBATCH --array=0-19%20
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=00:30:00
#SBATCH --output=logs/198_infer_sibling_ibd.%A_%a.out
#SBATCH --error=logs/198_infer_sibling_ibd.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python
bcftools=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/bcftools

call_root="${project}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/variant_calls/by_chromosome"
sample_file="${project}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata/bam_set/full_sib_14_samples.samples"
output_root="${project}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/ibd_transmission/by_chromosome"

chromosome=$((SLURM_ARRAY_TASK_ID + 1))
started_epoch=$(date +%s)

cd "${project}"
mkdir -p "${output_root}"

printf 'STAGE_START\tstage=infer_sibling_ibd\tchromosome=%s\tjob=%s\ttask=%s\thost=%s\ttime=%s\n' \
    "${chromosome}" \
    "${SLURM_ARRAY_JOB_ID}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/198_infer_sibling_ibd_by_chromosome.py \
    --chromosome "${chromosome}" \
    --bcftools "${bcftools}" \
    --call-root "${call_root}" \
    --sample-file "${sample_file}" \
    --output-root "${output_root}" \
    --window-size 1000000 \
    --minimum-depth 5 \
    --minimum-pl-gap 20 \
    --classification-error 0.01 \
    --minimum-informative-sites 10 \
    --minimum-absolute-llr 5

ended_epoch=$(date +%s)

printf 'STAGE_END\tstage=infer_sibling_ibd\tchromosome=%s\tstatus=0\telapsed_seconds=%s\thost=%s\ttime=%s\n' \
    "${chromosome}" \
    "$((ended_epoch - started_epoch))" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"
