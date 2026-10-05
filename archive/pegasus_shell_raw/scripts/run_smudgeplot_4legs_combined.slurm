#!/bin/bash
#SBATCH --job-name=smudge4legs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/smudge4legs.%j.out
#SBATCH --error=logs/smudge4legs.%j.err

set -euo pipefail

PROJECT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

WORK="$PROJECT/analysis/smudgeplot_4legs_combined"
mkdir -p "$WORK"

FASTQ_ROOT="$PROJECT/analysis/genomewide_mixing_input/fastq"

FASTQS=(
    "$FASTQ_ROOT/TN_Left_T2/TN_Left_T2_1.fastq.gz"
    "$FASTQ_ROOT/TN_Left_T2/TN_Left_T2_2.fastq.gz"

    "$FASTQ_ROOT/TN_Left_T3/TN_Left_T3_1.fastq.gz"
    "$FASTQ_ROOT/TN_Left_T3/TN_Left_T3_2.fastq.gz"

    "$FASTQ_ROOT/TN_Right_T2/TN_Right_T2_1.fastq.gz"
    "$FASTQ_ROOT/TN_Right_T2/TN_Right_T2_2.fastq.gz"

    "$FASTQ_ROOT/TN_Right_T3/TN_Right_T3_1.fastq.gz"
    "$FASTQ_ROOT/TN_Right_T3/TN_Right_T3_2.fastq.gz"
)

echo "=== INPUT FASTQS ==="

for f in "${FASTQS[@]}"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing FASTQ: $f" >&2
        exit 1
    fi

    ls -lh "$f"
done


# ------------------------------------------------------------
# Recover the same smudgeplot environment used previously.
# ------------------------------------------------------------

export PATH="$HOME/.local/share/mamba/envs/smudgeplot/bin:$PATH"

FASTK="$(command -v FastK || true)"
SMUDGE="$(command -v smudgeplot || true)"

if [[ -z "$FASTK" ]]; then
    echo "ERROR: FastK not found." >&2
    exit 1
fi

if [[ -z "$SMUDGE" ]]; then
    echo "ERROR: smudgeplot not found." >&2
    exit 1
fi

echo
echo "FASTK=$FASTK"
echo "SMUDGEPLOT=$SMUDGE"

"$SMUDGE" --version || true


# ------------------------------------------------------------
# 1. FastK
#
# Same parameters as the old workflow.
# Counting all eight FASTQs together is equivalent to the old
# merged-R1/merged-R2 setup for k-mer counting.
# ------------------------------------------------------------

echo
echo "=== FASTK ==="

"$FASTK" \
    -v \
    -t4 \
    -k22 \
    -M100 \
    -T60 \
    "${FASTQS[@]}" \
    -N"$WORK/FastK_Table"


# ------------------------------------------------------------
# 2. Heterozygous k-mer pairs
#
# Same old settings:
#   -L 12
#   -t 60
# ------------------------------------------------------------

echo
echo "=== SMUDGEPLOT HETMERS ==="

"$SMUDGE" hetmers \
    -L 12 \
    -t 60 \
    -o "$WORK/kmerpairs" \
    --verbose \
    "$WORK/FastK_Table"


if [[ ! -s "$WORK/kmerpairs.smu" ]]; then
    echo "ERROR: kmerpairs.smu was not created." >&2
    exit 1
fi


# ------------------------------------------------------------
# 3. Standard smudgeplot output as a sanity/reference copy
# ------------------------------------------------------------

echo
echo "=== SMUDGEPLOT ALL ==="

"$SMUDGE" all \
    -o "$WORK/output" \
    "$WORK/kmerpairs.smu"


# ------------------------------------------------------------
# Audit
# ------------------------------------------------------------

echo
echo "=== OUTPUTS ==="

ls -lh \
    "$WORK"/FastK_Table* \
    "$WORK"/kmerpairs* \
    "$WORK"/output* \
    2>/dev/null || true

echo
echo "=== FIRST SMU ROWS ==="
head -10 "$WORK/kmerpairs.smu"

echo
echo "COMPLETE"
