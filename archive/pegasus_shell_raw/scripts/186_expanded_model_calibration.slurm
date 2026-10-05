#!/usr/bin/env bash
#SBATCH --job-name=gyn186
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/186_expanded_calibration.%A_%a.out
#SBATCH --error=logs/186_expanded_calibration.%A_%a.err

set -euo pipefail

start=$SECONDS
trap 'rc=$?; printf "STAGE_END status=%s elapsed_seconds=%s\n" "$rc" "$((SECONDS-start))"' EXIT

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPLBACKEND=Agg

mode=${MODE:?MODE not supplied}

printf 'STAGE_START mode=%s job=%s task=%s host=%s time=%s\n' \
    "$mode" "$SLURM_JOB_ID" "${SLURM_ARRAY_TASK_ID:-NA}" \
    "$(hostname)" "$(date -Is)"

if [[ "$mode" == prepare || "$mode" == simulate ]]
then
    "$python_bin" -u \
        scripts/186_expanded_model_calibration.py \
        "$mode" \
        --task "$SLURM_ARRAY_TASK_ID"
elif [[ "$mode" == observed || "$mode" == collect || "$mode" == plot ]]
then
    "$python_bin" -u \
        scripts/186_expanded_model_calibration.py \
        "$mode"
else
    printf 'ERROR: unsupported MODE=%s\n' "$mode" >&2
    false
fi
