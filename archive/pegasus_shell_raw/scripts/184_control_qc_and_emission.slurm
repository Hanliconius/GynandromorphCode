#!/usr/bin/env bash
#SBATCH --job-name=gyn184qc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=08:00:00
#SBATCH --output=logs/184_control_qc.%A_%a.out
#SBATCH --error=logs/184_control_qc.%A_%a.err

set -euo pipefail

start=$SECONDS
trap 'rc=$?; printf "STAGE_END status=%s elapsed_seconds=%s\n" "$rc" "$((SECONDS-start))"' EXIT

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

printf 'STAGE_START mode=%s job=%s host=%s time=%s\n' \
    "${MODE:-scan}" "$SLURM_JOB_ID" "$(hostname)" "$(date -Is)"

if [[ "${MODE:-scan}" == fit ]]
then
    "$python_bin" -u scripts/184_control_qc_and_emission.py fit
else
    "$python_bin" -u scripts/184_control_qc_and_emission.py scan \
        --chromosome "$((SLURM_ARRAY_TASK_ID + 1))"
fi
