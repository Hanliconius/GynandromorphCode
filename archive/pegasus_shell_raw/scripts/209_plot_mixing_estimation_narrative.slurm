#!/usr/bin/env bash
#SBATCH --job-name=gyn209mix
#SBATCH --partition=cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --output=logs/209_plot_mixing_estimation_narrative.%j.out
#SBATCH --error=logs/209_plot_mixing_estimation_narrative.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

cd "${project}"

printf 'STAGE_START\tstage=plot_mixing_estimation_narrative\tjob=%s\thost=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/209_plot_mixing_estimation_narrative.py

printf 'STAGE_END\tstage=plot_mixing_estimation_narrative\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
