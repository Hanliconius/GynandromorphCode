#!/usr/bin/env bash
#SBATCH --job-name=gyn195c
#SBATCH --partition=nano
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --chdir=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
#SBATCH --output=logs/195c_validate_full_sib_joint_calls.%j.out
#SBATCH --error=logs/195c_validate_full_sib_joint_calls.%j.err

set -euo pipefail

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

printf 'STAGE_START\tstage=aggregate_full_sib_joint_call_validation\tjob=%s\thost=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/194e_validate_full_sib_call_subsets.py \
    --aggregate

printf 'ENDPOINT=FULL_SIB_JOINT_CALLS_READY\n'

printf 'STAGE_END\tstage=aggregate_full_sib_joint_call_validation\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
