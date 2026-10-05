#!/bin/bash
#SBATCH --job-name=smudgeFull
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/smudgeFull_%j.out
#SBATCH --error=logs/smudgeFull_%j.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
ARR="${ROOT}/analysis/smudgeplot_arrangements"
OUT="${ROOT}/analysis/smudgeplot_full_detail_bundle"

export PATH="$HOME/.local/share/mamba/envs/smudgeplot/bin:${PATH}"

mkdir -p "$OUT"
cd "$ROOT"

echo "START = $(date)"
echo "HOST  = $(hostname)"
echo


# ============================================================
# Generate native full-detail PDF outputs for the four TN sets
#
# We use the existing kmerpairs.smu files.
# These retain the original L=12 hetmer filtering.
# smudgeplot all is rerun only to request PDF output.
# ============================================================

run_tn () {
    name="$1"
    dir="$2"

    echo
    echo "============================================================"
    echo "TN DATASET: ${name}"
    echo "DIR: ${dir}"
    echo "============================================================"

    cd "$dir"

    test -s kmerpairs.smu

    if [[ ! -s full_detail_smudgeplot.pdf ]]; then
        smudgeplot all \
            --format pdf \
            -o full_detail \
            kmerpairs.smu \
            2>&1 | tee "${OUT}/${name}.smudgeplot_all.log"
    else
        echo "Existing full_detail_smudgeplot.pdf found; not rerunning."
    fi

    test -s full_detail_smudgeplot.pdf
    test -s full_detail_smudgeplot_log10.pdf

    cd "$ROOT"
}


run_tn \
    TN_Left \
    "${ARR}/TN_Left"

run_tn \
    TN_Right \
    "${ARR}/TN_Right"

run_tn \
    TN_T2_merge \
    "${ARR}/TN_T2_merge"

run_tn \
    TN_T3_merge \
    "${ARR}/TN_T3_merge"


# ============================================================
# Calibrated controls already generated
# ============================================================

DIPDIR="${ARR}/ERR1051012_krista"
PSEUDODIR="${ARR}/pseudoTetraploid_ERR1307111_ERR1307112"

DIP_MAIN="${DIPDIR}/output_L4_cov8.01_smudgeplot.pdf"
DIP_LOG="${DIPDIR}/output_L4_cov8.01_smudgeplot_log10.pdf"

PSEUDO_MAIN="${PSEUDODIR}/output_cov24.85_smudgeplot.pdf"
PSEUDO_LOG="${PSEUDODIR}/output_cov24.85_smudgeplot_log10.pdf"

test -s "$DIP_MAIN"
test -s "$DIP_LOG"
test -s "$PSEUDO_MAIN"
test -s "$PSEUDO_LOG"


# ============================================================
# Copy into one clearly named bundle directory
# ============================================================

cp -f \
    "$DIP_MAIN" \
    "${OUT}/01_diploid_ERR1051012_full.pdf"

cp -f \
    "$DIP_LOG" \
    "${OUT}/01_diploid_ERR1051012_log10.pdf"

cp -f \
    "$PSEUDO_MAIN" \
    "${OUT}/02_pseudoTet_ERR1307111_ERR1307112_full.pdf"

cp -f \
    "$PSEUDO_LOG" \
    "${OUT}/02_pseudoTet_ERR1307111_ERR1307112_log10.pdf"

cp -f \
    "${ARR}/TN_Left/full_detail_smudgeplot.pdf" \
    "${OUT}/03_TN_Left_full.pdf"

cp -f \
    "${ARR}/TN_Left/full_detail_smudgeplot_log10.pdf" \
    "${OUT}/03_TN_Left_log10.pdf"

cp -f \
    "${ARR}/TN_Right/full_detail_smudgeplot.pdf" \
    "${OUT}/04_TN_Right_full.pdf"

cp -f \
    "${ARR}/TN_Right/full_detail_smudgeplot_log10.pdf" \
    "${OUT}/04_TN_Right_log10.pdf"

cp -f \
    "${ARR}/TN_T2_merge/full_detail_smudgeplot.pdf" \
    "${OUT}/05_TN_T2_merge_full.pdf"

cp -f \
    "${ARR}/TN_T2_merge/full_detail_smudgeplot_log10.pdf" \
    "${OUT}/05_TN_T2_merge_log10.pdf"

cp -f \
    "${ARR}/TN_T3_merge/full_detail_smudgeplot.pdf" \
    "${OUT}/06_TN_T3_merge_full.pdf"

cp -f \
    "${ARR}/TN_T3_merge/full_detail_smudgeplot_log10.pdf" \
    "${OUT}/06_TN_T3_merge_log10.pdf"


# ============================================================
# Manifest
# ============================================================

cat > "${OUT}/manifest.tsv" <<EOF
orderdatasetrolehetmer_Lcoverage_handlinginput
1ERR1051012diploid_control4manual_1n=8.01${DIPDIR}/kmerpairs_L4.smu
2ERR1307111+ERR1307112pseudo_tetraploid12manual_1n=24.85${PSEUDODIR}/kmerpairs.smu
3TN_Leftgynandromorph_left12auto${ARR}/TN_Left/kmerpairs.smu
4TN_Rightgynandromorph_right12auto${ARR}/TN_Right/kmerpairs.smu
5TN_T2_mergegynandromorph_T2_L+R12auto${ARR}/TN_T2_merge/kmerpairs.smu
6TN_T3_mergegynandromorph_T3_L+R12auto${ARR}/TN_T3_merge/kmerpairs.smu
EOF


# ============================================================
# Concatenate available reports into one readable text file
# ============================================================

REPORT="${OUT}/all_smudge_reports.txt"
: > "$REPORT"

add_report () {
    label="$1"
    file="$2"

    {
        echo
        echo "================================================================"
        echo "$label"
        echo "================================================================"
        echo "FILE: $file"
        echo

        if [[ -s "$file" ]]; then
            cat "$file"
        else
            echo "REPORT FILE NOT FOUND"
        fi
    } >> "$REPORT"
}


add_report \
    "01 DIPLOID CONTROL - ERR1051012, L4, cov 8.01" \
    "${DIPDIR}/output_L4_cov8.01.smudge_report.tsv"

add_report \
    "02 PSEUDO-TETRAPLOID - ERR1307111 + ERR1307112, cov 24.85" \
    "${PSEUDODIR}/output_cov24.85.smudge_report.tsv"

add_report \
    "03 TN LEFT" \
    "${ARR}/TN_Left/output.smudge_report.tsv"

add_report \
    "04 TN RIGHT" \
    "${ARR}/TN_Right/output.smudge_report.tsv"

add_report \
    "05 TN T2 MERGE" \
    "${ARR}/TN_T2_merge/output.smudge_report.tsv"

add_report \
    "06 TN T3 MERGE" \
    "${ARR}/TN_T3_merge/output.smudge_report.tsv"


# ============================================================
# PDF merge helper
# ============================================================

merge_pdf () {
    outfile="$1"
    shift

    rm -f "$outfile"

    if command -v pdfunite >/dev/null 2>&1; then
        pdfunite "$@" "$outfile"

    elif command -v qpdf >/dev/null 2>&1; then
        qpdf --empty --pages "$@" -- "$outfile"

    elif command -v gs >/dev/null 2>&1; then
        gs \
            -dBATCH \
            -dNOPAUSE \
            -q \
            -sDEVICE=pdfwrite \
            -sOutputFile="$outfile" \
            "$@"
    else
        echo "ERROR: no PDF merging program found." >&2
        return 1
    fi

    test -s "$outfile"
}


# ============================================================
# Six-page scrollable PDF:
# one COMPLETE native Smudgeplot page per dataset
# ============================================================

merge_pdf \
    "${OUT}/smudge_full_detail_6page.pdf" \
    "${OUT}/01_diploid_ERR1051012_full.pdf" \
    "${OUT}/02_pseudoTet_ERR1307111_ERR1307112_full.pdf" \
    "${OUT}/03_TN_Left_full.pdf" \
    "${OUT}/04_TN_Right_full.pdf" \
    "${OUT}/05_TN_T2_merge_full.pdf" \
    "${OUT}/06_TN_T3_merge_full.pdf"


# ============================================================
# Twelve-page diagnostic version:
# normal + log10 for every dataset
# ============================================================

merge_pdf \
    "${OUT}/smudge_full_detail_with_log10_12page.pdf" \
    "${OUT}/01_diploid_ERR1051012_full.pdf" \
    "${OUT}/01_diploid_ERR1051012_log10.pdf" \
    "${OUT}/02_pseudoTet_ERR1307111_ERR1307112_full.pdf" \
    "${OUT}/02_pseudoTet_ERR1307111_ERR1307112_log10.pdf" \
    "${OUT}/03_TN_Left_full.pdf" \
    "${OUT}/03_TN_Left_log10.pdf" \
    "${OUT}/04_TN_Right_full.pdf" \
    "${OUT}/04_TN_Right_log10.pdf" \
    "${OUT}/05_TN_T2_merge_full.pdf" \
    "${OUT}/05_TN_T2_merge_log10.pdf" \
    "${OUT}/06_TN_T3_merge_full.pdf" \
    "${OUT}/06_TN_T3_merge_log10.pdf"


# ============================================================
# Audit
# ============================================================

echo
echo "================ BUNDLE CONTENTS ================"
ls -lh "$OUT"

echo
echo "================ MANIFEST ================"
column -t -s $'\t' "${OUT}/manifest.tsv" || cat "${OUT}/manifest.tsv"

echo
echo "FINISH = $(date)"
