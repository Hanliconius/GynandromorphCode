#!/usr/bin/env bash
#SBATCH --job-name=gyn215matrix
#SBATCH --partition=nano
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/215_plot_founder_distribution_matrix.%j.out
#SBATCH --error=logs/215_plot_founder_distribution_matrix.%j.err

set -euo pipefail

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

printf 'STAGE_START\tstage=founder_distribution_matrix\tjob=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/215_plot_founder_distribution_matrix.py

printf 'STAGE_END\tstage=founder_distribution_matrix\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
