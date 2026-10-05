#!/usr/bin/env bash
#SBATCH --job-name=gyn213cmp
#SBATCH --partition=nano
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/213_plot_sequential_model_compatibility.%j.out
#SBATCH --error=logs/213_plot_sequential_model_compatibility.%j.err

set -euo pipefail

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

printf 'STAGE_START\tstage=sequential_model_compatibility\tjob=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/213_plot_sequential_model_compatibility.py

printf 'STAGE_END\tstage=sequential_model_compatibility\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
