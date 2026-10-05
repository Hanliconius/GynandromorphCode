#!/usr/bin/env bash
#SBATCH --job-name=gyn199sibpl
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=logs/199_collect_and_plot_full_sib_ibd.%j.out
#SBATCH --error=logs/199_collect_and_plot_full_sib_ibd.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

cd "${project}"

started_epoch=$(date +%s)

printf 'STAGE_START\tstage=collect_and_plot_full_sib_ibd\tjob=%s\thost=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/199_collect_and_plot_full_sib_ibd.py

ended_epoch=$(date +%s)

printf 'STAGE_END\tstage=collect_and_plot_full_sib_ibd\tstatus=0\telapsed_seconds=%s\thost=%s\ttime=%s\n' \
    "$((ended_epoch - started_epoch))" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"
