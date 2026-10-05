#!/bin/bash
#SBATCH --job-name=smudgePT
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/smudgePT_%j.out
#SBATCH --error=logs/smudgePT_%j.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
DIR="${ROOT}/analysis/smudgeplot_arrangements/pseudoTetraploid_ERR1307111_ERR1307112"

SMUDGE_BIN="$HOME/.local/share/mamba/envs/smudgeplot/bin"
export PATH="${SMUDGE_BIN}:${PATH}"

cd "$DIR"

smudgeplot all \
    -cov 24.85 \
    --format pdf \
    -o output_cov24.85 \
    kmerpairs.smu

ls -lh output_cov24.85*
