#!/bin/bash
#SBATCH --job-name=smudgeFix
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/smudgeFix_%j.out
#SBATCH --error=logs/smudgeFix_%j.err

set -euo pipefail

ROOT=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
export PATH="$HOME/.local/share/mamba/envs/smudgeplot/bin:$PATH"

OUT="$ROOT/analysis/smudgeplot_standard_rerun_20261002"
CTRL="$ROOT/analysis/smudgeplot_arrangements/ERR1051012_krista"
FOUR="$ROOT/analysis/smudgeplot_4legs_combined"

mkdir -p "$OUT"/{ERR1051012_diploid,TN_all4_legs}
mkdir -p "$ROOT/logs"

echo "smudgeplot:"
smudgeplot -v

########################################################################
# Preconditions
########################################################################

for x in \
    "$CTRL/FastK_Table" \
    "$FOUR/FastK_Table"
do
    if [[ ! -e "$x" && ! -e "${x}.ktab" ]]; then
        echo "ERROR: cannot find FastK database: $x" >&2
        exit 1
    fi
done


########################################################################
# 1. DIPLOID CONTROL: ERR1051012
#
# Re-extract at progressively cleaner lower cutoffs.
# For each extraction make:
#   - automatic coverage plot
#   - manual 1n = 8.01 plot
########################################################################

for L in 6 8 12
do
    D="$OUT/ERR1051012_diploid/L${L}"
    mkdir -p "$D"

    echo
    echo "============================================================"
    echo "ERR1051012 diploid: L=$L"
    echo "============================================================"

    if [[ ! -s "$D/kmerpairs.smu" ]]; then
        (
            cd "$D"
            smudgeplot hetmers \
                -L "$L" \
                -t 16 \
                -o kmerpairs \
                --verbose \
                "$CTRL/FastK_Table" \
                > hetmers.log 2>&1
        )
    fi

    test -s "$D/kmerpairs.smu"

    (
        cd "$D"

        smudgeplot all \
            --format pdf \
            --json_report \
            -ylim 120 \
            -t "ERR1051012 diploid; L=${L}; automatic coverage" \
            -o output_auto \
            kmerpairs.smu \
            > smudgeplot_auto.log 2>&1

        smudgeplot all \
            --format pdf \
            --json_report \
            -cov 8.01 \
            -ylim 120 \
            -t "ERR1051012 diploid; L=${L}; 1n=8.01" \
            -o output_cov8.01 \
            kmerpairs.smu \
            > smudgeplot_cov8.01.log 2>&1
    )

done


########################################################################
# 2. ALL FOUR LEGS TOGETHER
#
# Regenerate the SMU from the existing four-leg FastK table at L=12,
# then create stock Smudgeplot output rather than the old custom figure.
########################################################################

D="$OUT/TN_all4_legs/L12"
mkdir -p "$D"

echo
echo "============================================================"
echo "TN all four legs: L=12"
echo "============================================================"

if [[ ! -s "$D/kmerpairs.smu" ]]; then
    (
        cd "$D"
        smudgeplot hetmers \
            -L 12 \
            -t 16 \
            -o kmerpairs \
            --verbose \
            "$FOUR/FastK_Table" \
            > hetmers.log 2>&1
    )
fi

test -s "$D/kmerpairs.smu"

(
    cd "$D"

    smudgeplot all \
        --format pdf \
        --json_report \
        -ylim 240 \
        -t "TN all four legs; L=12; automatic coverage" \
        -o output_auto \
        kmerpairs.smu \
        > smudgeplot_auto.log 2>&1

    smudgeplot all \
        --format pdf \
        --json_report \
        -cov 20 \
        -ylim 240 \
        -t "TN all four legs; L=12; 1n=20" \
        -o output_cov20 \
        kmerpairs.smu \
        > smudgeplot_cov20.log 2>&1
)

########################################################################
# Collect just the useful standard PDFs in one download directory.
########################################################################

DL="$OUT/download"
mkdir -p "$DL"
rm -f "$DL"/*.pdf

for L in 6 8 12
do
    cp "$OUT/ERR1051012_diploid/L${L}/output_auto_smudgeplot.pdf" \
       "$DL/ERR1051012_diploid_L${L}_AUTO_smudgeplot.pdf"

    cp "$OUT/ERR1051012_diploid/L${L}/output_cov8.01_smudgeplot.pdf" \
       "$DL/ERR1051012_diploid_L${L}_cov8.01_smudgeplot.pdf"
done

cp "$OUT/TN_all4_legs/L12/output_auto_smudgeplot.pdf" \
   "$DL/TN_all4_legs_L12_AUTO_smudgeplot.pdf"

cp "$OUT/TN_all4_legs/L12/output_cov20_smudgeplot.pdf" \
   "$DL/TN_all4_legs_L12_cov20_smudgeplot.pdf"

echo
echo "============================================================"
echo "FINISHED"
echo "============================================================"
ls -lh "$DL"
