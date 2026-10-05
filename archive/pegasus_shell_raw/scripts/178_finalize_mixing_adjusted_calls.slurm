#!/usr/bin/env bash
#SBATCH --job-name=gyn178ad
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=logs/178_finalize_mixing_adjusted_calls.%A_%a.out
#SBATCH --error=logs/178_finalize_mixing_adjusted_calls.%A_%a.err

set -euo pipefail

stage_name=finalize_mixing_adjusted_calls
start_seconds=${SECONDS}
printf 'STAGE_START\tstage=%s\tjob=%s\ttask=%s\thost=%s\ttime=%s\n' \
    "${stage_name}" "${SLURM_JOB_ID:-NA}" "${SLURM_ARRAY_TASK_ID:-NA}" \
    "$(hostname)" "$(date --iso-8601=seconds)"
trap 'status=$?; printf "STAGE_END\tstage=%s\tstatus=%s\telapsed_seconds=%s\thost=%s\ttime=%s\n" "${stage_name}" "${status}" "$((SECONDS - start_seconds))" "$(hostname)" "$(date --iso-8601=seconds)"' EXIT

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

"${python_bin}" \
    "${project}/scripts/173_genomewide_mixing_models.py" \
    adjust \
    --task "${SLURM_ARRAY_TASK_ID}"
