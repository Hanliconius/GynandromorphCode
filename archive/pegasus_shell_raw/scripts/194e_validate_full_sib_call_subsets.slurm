#!/usr/bin/env bash
#SBATCH --job-name=gyn194e
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --array=0-20%21
#SBATCH --chdir=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
#SBATCH --output=logs/194e_validate_full_sib_call_subsets.%A_%a.out
#SBATCH --error=logs/194e_validate_full_sib_call_subsets.%A_%a.err

set -euo pipefail

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

printf 'STAGE_START\tstage=validate_full_sib_call_subset\tjob=%s\ttask=%s\tchromosome=%s\thost=%s\ttime=%s\n' \
    "${SLURM_ARRAY_JOB_ID}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "$((SLURM_ARRAY_TASK_ID + 1))" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/194e_validate_full_sib_call_subsets.py \
    --task-id "${SLURM_ARRAY_TASK_ID}"

printf 'STAGE_END\tstage=validate_full_sib_call_subset\tstatus=0\ttask=%s\ttime=%s\n' \
    "${SLURM_ARRAY_TASK_ID}" \
    "$(date --iso-8601=seconds)"
