#!/usr/bin/env bash
#SBATCH --job-name=gyn211D4
#SBATCH --partition=nano
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/211_plot_collapsed_Dennis_genotypes.%j.out
#SBATCH --error=logs/211_plot_collapsed_Dennis_genotypes.%j.err

set -euo pipefail

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

printf 'STAGE_START\tstage=collapsed_Dennis_genotypes\tjob=%s\thost=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/211_plot_collapsed_Dennis_genotypes.py

printf 'STAGE_END\tstage=collapsed_Dennis_genotypes\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
