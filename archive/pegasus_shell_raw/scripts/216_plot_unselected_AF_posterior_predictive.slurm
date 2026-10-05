#!/usr/bin/env bash
#SBATCH --job-name=gyn216rawAF
#SBATCH --partition=cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/216_plot_unselected_AF_posterior_predictive.%j.out
#SBATCH --error=logs/216_plot_unselected_AF_posterior_predictive.%j.err

set -euo pipefail

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

printf 'STAGE_START\tstage=unselected_AF_posterior_predictive\tjob=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/216_plot_unselected_AF_posterior_predictive.py

printf 'STAGE_END\tstage=unselected_AF_posterior_predictive\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
