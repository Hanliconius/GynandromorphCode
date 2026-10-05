#!/bin/bash
#SBATCH --job-name=smudge4leg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/smudge4leg_%j.out
#SBATCH --error=logs/smudge4leg_%j.err

set -u

ROOT=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
export PATH="$HOME/.local/share/mamba/envs/smudgeplot/bin:$PATH"

FOUR="$ROOT/analysis/smudgeplot_4legs_combined"
OUT="$ROOT/analysis/smudgeplot_standard_rerun_20261002/TN_all4_legs/L12"

# Final comparison/download directory we were already using
DL="$ROOT/analysis/smudgeplot_lowerL_20261002/L6_download"

mkdir -p "$OUT" "$DL"

echo "============================================================"
echo "FOUR LEGS COMBINED — STANDARD SMUDGEPLOT, L=12"
echo "============================================================"

cd "$OUT"

rm -f kmerpairs.smu

smudgeplot hetmers \
    -L 12 \
    -t 16 \
    -o kmerpairs \
    --verbose \
    "$FOUR/FastK_Table" \
    > hetmers.log 2>&1

if [[ ! -s kmerpairs.smu ]]; then
    echo "ERROR: four-leg L12 extraction produced no usable SMU"
    cat hetmers.log
    exit 1
fi

ls -lh kmerpairs.smu

smudgeplot all \
    --format pdf \
    --json_report \
    -ylim 240 \
    -t "TN all four legs; L=12; automatic coverage" \
    -o output_auto \
    kmerpairs.smu \
    > smudgeplot_auto.log 2>&1

if [[ -s output_auto_smudgeplot.pdf ]]; then
    cp -f output_auto_smudgeplot.pdf \
        "$DL/TN_all4_legs_L12_smudgeplot.pdf"
fi

# The already-successful final diploid control:
cp -f \
"$ROOT/analysis/smudgeplot_standard_rerun_20261002/ERR1051012_diploid/L6/output_auto_smudgeplot.pdf" \
"$DL/ERR1051012_diploid_control_L6_smudgeplot.pdf"

echo
echo "============================================================"
echo "FOUR-LEG LOG"
echo "============================================================"
cat smudgeplot_auto.log

echo
echo "============================================================"
echo "FINAL FILES"
echo "============================================================"
ls -lh \
"$DL/ERR1051012_diploid_control_L6_smudgeplot.pdf" \
"$DL/TN_all4_legs_L12_smudgeplot.pdf"
