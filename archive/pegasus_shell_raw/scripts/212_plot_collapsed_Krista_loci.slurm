#!/usr/bin/env bash
#SBATCH --job-name=gyn212K5
#SBATCH --partition=nano
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/212_plot_collapsed_Krista_loci.%j.out
#SBATCH --error=logs/212_plot_collapsed_Krista_loci.%j.err

set -euo pipefail

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

printf 'STAGE_START\tstage=collapsed_Krista_loci\tjob=%s\thost=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/212_plot_collapsed_Krista_loci.py

printf 'STAGE_END\tstage=collapsed_Krista_loci\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
