#!/usr/bin/env bash
#SBATCH --job-name=gyn214dist
#SBATCH --partition=nano
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/214_plot_founder_distribution_ladder.%j.out
#SBATCH --error=logs/214_plot_founder_distribution_ladder.%j.err

set -euo pipefail

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

printf 'STAGE_START\tstage=founder_distribution_ladder\tjob=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/214_plot_founder_distribution_ladder.py

printf 'STAGE_END\tstage=founder_distribution_ladder\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
