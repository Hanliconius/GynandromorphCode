#!/usr/bin/env bash
#SBATCH --job-name=gyn208mix
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/208_plot_mixing_number_mechanism.%j.out
#SBATCH --error=logs/208_plot_mixing_number_mechanism.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

cd "${project}"

printf 'STAGE_START\tstage=plot_mixing_number_mechanism\tjob=%s\thost=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/208_plot_mixing_number_mechanism.py

printf 'STAGE_END\tstage=plot_mixing_number_mechanism\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
