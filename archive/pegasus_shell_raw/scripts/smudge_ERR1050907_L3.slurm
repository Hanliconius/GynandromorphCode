#!/bin/bash
#SBATCH --job-name=smg0907
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/smg0907_L3_%j.out
#SBATCH --error=logs/smg0907_L3_%j.err

set -euo pipefail

export PATH="$HOME/.local/share/mamba/envs/smudgeplot/bin:${PATH}"

DIR="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/smudgeplot_arrangements/ERR1050907_krista"

cd "$DIR"

rm -f kmerpairs_L3.smu output_L3*

smudgeplot hetmers \
    -L 3 \
    -t 60 \
    -o kmerpairs_L3 \
    --verbose \
    FastK_Table

test -s kmerpairs_L3.smu

smudgeplot all \
    --format pdf \
    -o output_L3 \
    kmerpairs_L3.smu

echo
echo "===== REPORT ====="
cat output_L3.smudge_report.tsv

echo
echo "===== OUTPUTS ====="
ls -lh output_L3* kmerpairs_L3.smu
