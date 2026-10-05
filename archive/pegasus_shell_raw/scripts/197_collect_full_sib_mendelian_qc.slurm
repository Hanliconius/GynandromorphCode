#!/usr/bin/env bash
#SBATCH --job-name=gyn197sibag
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/197_collect_full_sib_mendelian_qc.%j.out
#SBATCH --error=logs/197_collect_full_sib_mendelian_qc.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

cd "${project}"

started_epoch=$(date +%s)

printf 'STAGE_START\tstage=collect_full_sib_mendelian_qc\tjob=%s\thost=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/197_collect_full_sib_mendelian_qc.py

ended_epoch=$(date +%s)

printf 'STAGE_END\tstage=collect_full_sib_mendelian_qc\tstatus=0\telapsed_seconds=%s\thost=%s\ttime=%s\n' \
    "$((ended_epoch - started_epoch))" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"
