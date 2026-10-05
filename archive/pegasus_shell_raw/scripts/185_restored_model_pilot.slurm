#!/usr/bin/env bash
#SBATCH --job-name=gyn185
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/185_restored_model_pilot.%A_%a.out
#SBATCH --error=logs/185_restored_model_pilot.%A_%a.err

set -euo pipefail

started=$SECONDS

trap 'status=$?; printf "STAGE_END status=%s elapsed_seconds=%s\n" "$status" "$((SECONDS-started))"' EXIT

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPLBACKEND=Agg

mode=${MODE:?MODE was not supplied}

printf 'STAGE_START mode=%s job=%s task=%s host=%s time=%s\n' \
    "$mode" \
    "$SLURM_JOB_ID" \
    "${SLURM_ARRAY_TASK_ID:-NA}" \
    "$(hostname)" \
    "$(date -Is)"

if [[ "$mode" == prepare ]]
then
    "$python_bin" -u \
        scripts/185_restored_model_pilot.py \
        prepare \
        --task "$SLURM_ARRAY_TASK_ID"
elif [[ "$mode" == observed ]]
then
    "$python_bin" -u \
        scripts/185_restored_model_pilot.py \
        observed
elif [[ "$mode" == pilot ]]
then
    "$python_bin" -u \
        scripts/185_restored_model_pilot.py \
        pilot \
        --task "$SLURM_ARRAY_TASK_ID"
elif [[ "$mode" == plot ]]
then
    "$python_bin" -u \
        scripts/185_restored_model_pilot.py \
        plot
else
    printf 'ERROR: unsupported MODE=%s\n' "$mode" >&2
    false
fi
