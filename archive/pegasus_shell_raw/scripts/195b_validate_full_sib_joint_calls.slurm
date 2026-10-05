#!/usr/bin/env bash
#SBATCH --job-name=gyn195sibv2
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:30:00
#SBATCH --output=logs/195b_validate_full_sib_joint_calls.%j.out
#SBATCH --error=logs/195b_validate_full_sib_joint_calls.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

cd "${project}"

started_epoch=$(date +%s)

printf 'STAGE_START\tstage=validate_full_sib_joint_calls_v2\tjob=%s\thost=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/195b_validate_full_sib_joint_calls.py

ended_epoch=$(date +%s)

printf 'STAGE_END\tstage=validate_full_sib_joint_calls_v2\tstatus=0\telapsed_seconds=%s\thost=%s\ttime=%s\n' \
    "$((ended_epoch - started_epoch))" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"
