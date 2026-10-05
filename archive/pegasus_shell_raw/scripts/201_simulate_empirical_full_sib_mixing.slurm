#!/usr/bin/env bash
#SBATCH --job-name=gyn201sib
#SBATCH --partition=cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=0-659%100
#SBATCH --chdir=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
#SBATCH --output=logs/201_simulate_empirical_full_sib_mixing.%A_%a.out
#SBATCH --error=logs/201_simulate_empirical_full_sib_mixing.%A_%a.err

set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

printf 'STAGE_START\tstage=simulate_empirical_full_sib_mixing\tjob=%s\ttask=%s\thost=%s\ttime=%s\n' \
    "${SLURM_ARRAY_JOB_ID}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/201_simulate_empirical_full_sib_mixing.py \
    --task "${SLURM_ARRAY_TASK_ID}"

printf 'STAGE_END\tstage=simulate_empirical_full_sib_mixing\tstatus=0\ttask=%s\ttime=%s\n' \
    "${SLURM_ARRAY_TASK_ID}" \
    "$(date --iso-8601=seconds)"
