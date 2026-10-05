#!/bin/bash
#SBATCH --job-name=smudgeKrista
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/smudgeKrista_%A_%a.out
#SBATCH --error=logs/smudgeKrista_%A_%a.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
FASTQROOT="${ROOT}/analysis/genomewide_mixing_input/fastq"
OUTROOT="${ROOT}/analysis/smudgeplot_arrangements"

SMUDGE_BIN="$HOME/.local/share/mamba/envs/smudgeplot/bin"
export PATH="${SMUDGE_BIN}:${PATH}"

SAMPLES=(
    ERR1050907
    ERR1050958
    ERR1050994
    ERR1051012
)

s="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
dir="${OUTROOT}/${s}_krista"

mkdir -p "$dir"
cd "$dir"

echo "SAMPLE = ${s}"
echo "DIR    = ${dir}"
echo

# Link original FASTQs
ln -sf "${FASTQROOT}/${s}/${s}_1.fastq.gz" "${s}_1.fastq.gz"
ln -sf "${FASTQROOT}/${s}/${s}_2.fastq.gz" "${s}_2.fastq.gz"

echo "INPUTS"
ls -lh ./*.fastq.gz
echo

# ------------------------------------------------------------
# 1. FastK — Krista's parameters
# ------------------------------------------------------------

if [[ ! -s FastK_Table.ktab || ! -s FastK_Table.hist ]]; then
    echo "===== FastK ====="

    FastK \
        -v \
        -t4 \
        -k22 \
        -M100 \
        -T60 \
        ./*.fastq.gz \
        -NFastK_Table
else
    echo "FastK already present; skipping."
fi

test -s FastK_Table.ktab
test -s FastK_Table.hist

# Human-readable histogram for audit
Histex -A -h1:100 FastK_Table.hist > FastK_Table.hist.tsv

echo
echo "===== HISTOGRAM 1-50x ====="
awk '$1 >= 1 && $1 <= 50 {print}' FastK_Table.hist.tsv
echo

# ------------------------------------------------------------
# 2. hetmers — Krista's parameters
# ------------------------------------------------------------

if [[ ! -s kmerpairs.smu ]]; then
    echo "===== smudgeplot hetmers ====="

    smudgeplot hetmers \
        -L 12 \
        -t 60 \
        -o kmerpairs \
        --verbose \
        FastK_Table
else
    echo "kmerpairs.smu already present; skipping."
fi

test -s kmerpairs.smu

# ------------------------------------------------------------
# 3. smudgeplot all
# ------------------------------------------------------------

echo
echo "===== smudgeplot all ====="

smudgeplot all \
    --format pdf \
    -o output \
    kmerpairs.smu

echo
echo "===== OUTPUTS ====="
ls -lh \
    FastK_Table.hist \
    FastK_Table.ktab \
    FastK_Table.hist.tsv \
    kmerpairs.smu \
    output* \
    2>/dev/null || true
