#!/bin/bash
#SBATCH --job-name=smudgeHi
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/smudgeHi_%A_%a.out
#SBATCH --error=logs/smudgeHi_%A_%a.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
BASE="${ROOT}/analysis/smudgeplot_arrangements"

SMUDGE_BIN="$HOME/.local/share/mamba/envs/smudgeplot/bin"
export PATH="${SMUDGE_BIN}:${PATH}"

DIRS=(
    "${BASE}/ERR1307111_single"
    "${BASE}/ERR1307112_single"
    "${BASE}/pseudoTetraploid_ERR1307111_ERR1307112"
)

dir="${DIRS[$SLURM_ARRAY_TASK_ID]}"
label="$(basename "$dir")"

echo "TASK       = ${SLURM_ARRAY_TASK_ID}"
echo "LABEL      = ${label}"
echo "DIRECTORY  = ${dir}"
echo "FastK      = $(command -v FastK)"
echo "smudgeplot = $(command -v smudgeplot)"
echo

cd "$dir"

echo "INPUT FASTQS"
ls -lh ./*.fastq.gz
echo

# ------------------------------------------------------------
# 1. FastK
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
    echo "FastK output already present; skipping recount."
fi

test -s FastK_Table.ktab
test -s FastK_Table.hist


# ------------------------------------------------------------
# 2. Export readable histogram for audit
# ------------------------------------------------------------

Histex -A -h1:100 FastK_Table.hist > FastK_Table.hist.tsv

echo
echo "===== TOP HISTOGRAM BINS ====="
sort -k2,2nr FastK_Table.hist.tsv | head -15
echo


# ------------------------------------------------------------
# 3. Heterozygous k-mer pairs
# ------------------------------------------------------------

if [[ ! -s kmerpairs.smu ]]; then

    echo "===== smudgeplot hetmers : L=12 ====="

    smudgeplot hetmers \
        -L 12 \
        -t 60 \
        -o kmerpairs \
        --verbose \
        FastK_Table
else
    echo "kmerpairs.smu already exists; skipping hetmers."
fi

test -s kmerpairs.smu


# ------------------------------------------------------------
# 4. Smudge inference + VECTOR PDF
# ------------------------------------------------------------

echo
echo "===== smudgeplot all ====="

rm -f \
    output.smudge_report.tsv \
    output_smudgeplot.pdf \
    output_smudgeplot_log10.pdf \
    output.sma

smudgeplot all \
    --format pdf \
    -o output \
    kmerpairs.smu

test -s output.smudge_report.tsv
test -s output_smudgeplot.pdf
test -s output_smudgeplot_log10.pdf


# ------------------------------------------------------------
# 5. Final audit
# ------------------------------------------------------------

echo
echo "===== SUCCESS: ${label} ====="

ls -lh \
    FastK_Table.hist \
    FastK_Table.ktab \
    FastK_Table.hist.tsv \
    kmerpairs.smu \
    output.smudge_report.tsv \
    output_smudgeplot.pdf \
    output_smudgeplot_log10.pdf

echo
echo "===== REPORT ====="
cat output.smudge_report.tsv
