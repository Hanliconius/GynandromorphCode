#!/usr/bin/env bash
#SBATCH --job-name=gyn176jk
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --output=logs/176_jackknife_genomewide_mixing.%A_%a.out
#SBATCH --error=logs/176_jackknife_genomewide_mixing.%A_%a.err

set -euo pipefail

stage_name=jackknife_genomewide_mixing
start_seconds=${SECONDS}
printf 'STAGE_START\tstage=%s\tjob=%s\ttask=%s\thost=%s\ttime=%s\n' \
    "${stage_name}" "${SLURM_JOB_ID:-NA}" "${SLURM_ARRAY_TASK_ID:-NA}" \
    "$(hostname)" "$(date --iso-8601=seconds)"
trap 'status=$?; printf "STAGE_END\tstage=%s\tstatus=%s\telapsed_seconds=%s\thost=%s\ttime=%s\n" "${stage_name}" "${status}" "$((SECONDS - start_seconds))" "$(hostname)" "$(date --iso-8601=seconds)"' EXIT

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

"${python_bin}" \
    "${project}/scripts/173_genomewide_mixing_models.py" \
    jackknife \
    --task "${SLURM_ARRAY_TASK_ID}"
