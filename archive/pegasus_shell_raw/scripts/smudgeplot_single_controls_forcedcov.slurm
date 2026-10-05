#!/bin/bash
#SBATCH --job-name=smudgeCov
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/smudgeCov_%A_%a.out
#SBATCH --error=logs/smudgeCov_%A_%a.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
BASE="${ROOT}/analysis/smudgeplot_arrangements"

export PATH="$HOME/.local/share/mamba/envs/smudgeplot/bin:${PATH}"

DIRS=(
    "${BASE}/ERR1307111_single"
    "${BASE}/ERR1307112_single"
)

LOGS=(
    "${ROOT}/logs/smudgeHiR_138358636_0.err"
    "${ROOT}/logs/smudgeHiR_138358636_1.err"
)

dir="${DIRS[$SLURM_ARRAY_TASK_ID]}"
log="${LOGS[$SLURM_ARRAY_TASK_ID]}"

cov=$(
    grep 'Best coverage to precision of 1/100 (just to be sure):' "$log" \
    | tail -1 \
    | awk '{print $NF}'
)

[[ -n "$cov" ]] || {
    echo "ERROR: could not recover coverage from $log"
    exit 1
}

cd "$dir"

echo "SAMPLE = $(basename "$dir")"
echo "1n coverage = $cov"

test -s kmerpairs.smu

smudgeplot all \
    -cov "$cov" \
    --format pdf \
    -o "output_cov${cov}" \
    kmerpairs.smu

echo
ls -lh output_cov"${cov}"*
