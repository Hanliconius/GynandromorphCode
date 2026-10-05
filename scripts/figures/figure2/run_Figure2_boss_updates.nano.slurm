#!/usr/bin/env bash
#SBATCH --job-name=Fig2boss
#SBATCH --partition=nano
#SBATCH --array=0-7%4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${FIGURE2_CODE_DIR:?FIGURE2_CODE_DIR must be exported by the submission helper}
PYTHON=${PYTHON:-$HOME/.local/share/mamba/envs/gynandrocalling/bin/python}
SMUDGE_BIN="$HOME/.local/share/mamba/envs/smudgeplot/bin"

OUTDIR="$PROJECT_ROOT/analysis/figure2/panels"
WORKROOT="$PROJECT_ROOT/analysis/figure2/boss_updates"
mkdir -p "$OUTDIR" "$WORKROOT"

task=${SLURM_ARRAY_TASK_ID}

echo "FIGURE2_BOSS_UPDATE_START"
echo "task=$task"
echo "host=$(hostname)"
echo "PROJECT_ROOT=$PROJECT_ROOT"
echo "CODE_DIR=$CODE_DIR"

validate_pdf_text () {
    local pdf="$1"

    if [[ ! -s "$pdf" ]]; then
        echo "ERROR: missing/empty PDF: $pdf" >&2
        exit 1
    fi

    if command -v pdftotext >/dev/null 2>&1; then
        if ! pdftotext "$pdf" - 2>/dev/null | grep -q '[[:alnum:]]'; then
            echo "ERROR: PDF contains no extractable text: $pdf" >&2
            exit 1
        fi
        echo "TEXT_CHECK=PASS $pdf"
    else
        echo "TEXT_CHECK=SKIPPED pdftotext_not_available $pdf"
    fi

    if command -v pdffonts >/dev/null 2>&1; then
        echo "PDF_FONTS $pdf"
        pdffonts "$pdf" || true
    fi
}

if (( task <= 5 )); then
    export PATH="$SMUDGE_BIN:$PATH"
    export MPLBACKEND=Agg

    case "$task" in
        0)
            key="C1_diploid_control"
            title="Diploid control"
            smu="$PROJECT_ROOT/analysis/smudgeplot_standard_rerun_20261002/ERR1051012_diploid/L6/kmerpairs.smu"
            dest="$OUTDIR/Fig2C1_diploid_control.pdf"
            ;;
        1)
            key="C2_pseudo_tetraploid"
            title="Pseudo-tetraploid"
            smu="$PROJECT_ROOT/analysis/smudgeplot_lowerL_20261002/pseudoTetraploid_ERR1307111_ERR1307112/L6/kmerpairs.smu"
            dest="$OUTDIR/Fig2C2_pseudo_tetraploid.pdf"
            ;;
        2)
            key="D1_left"
            title="Left"
            smu="$PROJECT_ROOT/analysis/smudgeplot_arrangements/TN_Left/kmerpairs.smu"
            dest="$OUTDIR/Fig2D1_left.pdf"
            ;;
        3)
            key="D2_right"
            title="Right"
            smu="$PROJECT_ROOT/analysis/smudgeplot_arrangements/TN_Right/kmerpairs.smu"
            dest="$OUTDIR/Fig2D2_right.pdf"
            ;;
        4)
            key="D3_T2_merge"
            title="T2 merge"
            smu="$PROJECT_ROOT/analysis/smudgeplot_lowerL_20261002/TN_T2_merge/L6/kmerpairs.smu"
            dest="$OUTDIR/Fig2D3_T2_merge.pdf"
            ;;
        5)
            key="D4_T3_merge"
            title="T3 merge"
            smu="$PROJECT_ROOT/analysis/smudgeplot_lowerL_20261002/TN_T3_merge/L6/kmerpairs.smu"
            dest="$OUTDIR/Fig2D4_T3_merge.pdf"
            ;;
    esac

    if [[ ! -s "$smu" ]]; then
        echo "ERROR: missing smudge input: $smu" >&2
        exit 1
    fi

    work="$WORKROOT/$key"
    rm -rf "$work"
    mkdir -p "$work/mplconfig"

    # Preserve native Smudgeplot aesthetics while forcing font embedding
    # as text in PDF/SVG outputs.
    cat > "$work/mplconfig/matplotlibrc" <<'MPL'
font.family: DejaVu Sans
pdf.fonttype: 42
ps.fonttype: 42
svg.fonttype: none
text.usetex: False
MPL

    export MPLCONFIGDIR="$work/mplconfig"

    prefix="$work/output"

    echo "SMUDGE_INPUT=$smu"
    echo "SMUDGE_AXES=x:0..0.5 y:0..120"
    echo "SMUDGE_STYLE=native smudgeplot all"
    echo "SMUDGE_TEXT=Type42/extractable"

    smudgeplot all         --format pdf         --json_report         -ylim 120         -col_ramp viridis         -t "$title"         -o "$prefix"         "$smu"

    source_pdf="${prefix}_smudgeplot.pdf"

    if [[ ! -s "$source_pdf" ]]; then
        echo "ERROR: Smudgeplot did not create $source_pdf" >&2
        find "$work" -maxdepth 1 -type f -printf '%f\t%s bytes\n' | sort >&2
        exit 1
    fi

    tmp_dest="${dest}.partial"
    cp "$source_pdf" "$tmp_dest"
    mv "$tmp_dest" "$dest"

    validate_pdf_text "$dest"

    echo "WROTE=$dest"

else
    if [[ ! -x "$PYTHON" ]]; then
        echo "ERROR: Python not executable: $PYTHON" >&2
        exit 1
    fi

    if [[ "$task" -eq 6 ]]; then
        script="$CODE_DIR/Fig2G_simulated_AF.py"
        dest="$OUTDIR/Fig2G.pdf"
    else
        script="$CODE_DIR/Fig2H_observed_AF.py"
        dest="$OUTDIR/Fig2H.pdf"
    fi

    if [[ ! -s "$script" ]]; then
        echo "ERROR: missing AF plotting script: $script" >&2
        exit 1
    fi

    rm -f "$dest"

    "$PYTHON" "$script"         --project-root "$PROJECT_ROOT"         --density-quantile 0.995         --origin-exclude 0.05

    validate_pdf_text "$dest"

    echo "AF_DENSITY_SCALE=LogNorm with upper saturation chosen from q=0.995 outside 0,0 corner"
    echo "WROTE=$dest"
fi

echo "FIGURE2_BOSS_UPDATE_COMPLETE"
