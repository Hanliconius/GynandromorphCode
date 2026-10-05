#!/usr/bin/env bash
#SBATCH --job-name=gyn200sib
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --chdir=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
#SBATCH --output=logs/200_prepare_empirical_sibling_model_inputs.%j.out
#SBATCH --error=logs/200_prepare_empirical_sibling_model_inputs.%j.err

set -euo pipefail

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

printf 'STAGE_START\tstage=prepare_empirical_sibling_model_inputs\tjob=%s\thost=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/200_prepare_empirical_sibling_model_inputs.py

printf 'STAGE_END\tstage=prepare_empirical_sibling_model_inputs\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
