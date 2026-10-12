#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

[[ -x "$PYTHON" ]] || { echo "ERROR: Python not executable: $PYTHON" >&2; exit 2; }
mkdir -p "$PROJECT_ROOT/logs" "$PROJECT_ROOT/analysis/yb_raw_from_scratch"
bash -n "$HERE/run_Yb_raw_download.nano.slurm"
bash -n "$HERE/run_Yb_hmros_reference.nano.slurm"
bash -n "$HERE/run_Yb_raw_align_coverage.nano.slurm"
bash -n "$HERE/run_Yb_raw_plot.nano.slurm"
"$PYTHON" -m py_compile "$HERE/download_Yb_raw_fastq.py" "$HERE/plot_Yb_raw_Krista_style.py"

REF_JOB=$(sbatch --parsable \
    --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT" \
    "$HERE/run_Yb_hmros_reference.nano.slurm")
DL_JOB=$(sbatch --parsable \
    --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,YB_CODE_DIR=$HERE" \
    "$HERE/run_Yb_raw_download.nano.slurm")
ALIGN_JOB=$(sbatch --parsable \
    --dependency="afterok:${REF_JOB}:${DL_JOB}" \
    --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,YB_CODE_DIR=$HERE" \
    "$HERE/run_Yb_raw_align_coverage.nano.slurm")
PLOT_JOB=$(sbatch --parsable \
    --dependency="afterok:${ALIGN_JOB}" \
    --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,YB_CODE_DIR=$HERE" \
    "$HERE/run_Yb_raw_plot.nano.slurm")

printf 'REFERENCE_JOB\t%s\nDOWNLOAD_ARRAY\t%s\nALIGN_ARRAY\t%s\nPLOT_JOB\t%s\n' \
    "$REF_JOB" "$DL_JOB" "$ALIGN_JOB" "$PLOT_JOB"
echo "OUTPUT_DIR\t$PROJECT_ROOT/analysis/yb_raw_from_scratch/HmRos/plots"
