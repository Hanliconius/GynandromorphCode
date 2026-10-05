#!/usr/bin/env bash
#SBATCH --job-name=gyn187
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/187_188_expanded_calibration.%A_%a.out
#SBATCH --error=logs/187_188_expanded_calibration.%A_%a.err

set -euo pipefail

started=$SECONDS
trap 'rc=$?; printf "STAGE_END status=%s elapsed_seconds=%s\n" "$rc" "$((SECONDS-started))"' EXIT

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

if [[ "$mode" == simulate ]]
then
    "$python_bin" -u \
        scripts/187_expand_simulations_and_window_evidence.py \
        simulate \
        --task "$SLURM_ARRAY_TASK_ID"
elif [[ "$mode" == collect ||
        "$mode" == evidence ||
        "$mode" == zdiag ]]
then
    "$python_bin" -u \
        scripts/187_expand_simulations_and_window_evidence.py \
        "$mode"
elif [[ "$mode" == narrative ]]
then
    "$python_bin" -u \
        scripts/188_plot_corrected_genomewide_narrative.py
else
    printf 'ERROR: unsupported MODE=%s\n' "$mode" >&2
    false
fi
