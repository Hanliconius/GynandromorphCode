#!/bin/bash
#SBATCH --job-name=smudgeHiR
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --array=0-1
#SBATCH --output=logs/smudgeHiR_%A_%a.out
#SBATCH --error=logs/smudgeHiR_%A_%a.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
BASE="${ROOT}/analysis/smudgeplot_arrangements"

SMUDGE_BIN="$HOME/.local/share/mamba/envs/smudgeplot/bin"
export PATH="${SMUDGE_BIN}:${PATH}"

DIRS=(
    "${BASE}/ERR1307111_single"
    "${BASE}/ERR1307112_single"
)

dir="${DIRS[$SLURM_ARRAY_TASK_ID]}"
label="$(basename "$dir")"

echo "LABEL = ${label}"
echo "DIR   = ${dir}"

cd "$dir"

# Remove only partial FastK products from the failed OOM run.
rm -f FastK_Table.hist FastK_Table.ktab

echo "===== FastK ====="

FastK \
    -v \
    -t4 \
    -k22 \
    -M100 \
    -T60 \
    ./*.fastq.gz \
    -NFastK_Table

test -s FastK_Table.ktab
test -s FastK_Table.hist

Histex -A -h1:100 FastK_Table.hist > FastK_Table.hist.tsv

echo "===== hetmers ====="

smudgeplot hetmers \
    -L 12 \
    -t 60 \
    -o kmerpairs \
    --verbose \
    FastK_Table

test -s kmerpairs.smu

echo "===== smudgeplot all ====="

smudgeplot all \
    --format pdf \
    -o output \
    kmerpairs.smu

echo
echo "===== FINAL FILES ====="
ls -lh FastK_Table.* kmerpairs.smu output* || true
