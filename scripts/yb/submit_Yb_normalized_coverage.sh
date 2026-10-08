#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

mkdir -p "$PROJECT_ROOT/logs"

echo "===== Yb NORMALIZED COVERAGE PREFLIGHT ====="

"$PYTHON" -m py_compile     "$HERE/plot_Yb_normalized_coverage.py"

bash -n     "$HERE/run_Yb_normalized_coverage.nano.slurm"

echo "PREFLIGHT=PASS"
echo "Historical inputs: Krista-Pipho/gynandromorph_2024 @ 99376deebc25f4e0bb7aa9927c7969c693591ef4"
echo "Coordinates: Hmel2.5 Hmel215003o"
echo "Broad view: 1.380-1.600 Mb"
echo "Normalization: each library / median coverage of local flanking windows"
echo

job=$(
    sbatch         --parsable         --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,YB_CODE_DIR=$HERE"         "$HERE/run_Yb_normalized_coverage.nano.slurm"
)
job=${job%%;*}

echo "JOBID=$job"
echo "OUT=$PROJECT_ROOT/logs/Yb_normalized_coverage.$job.out"
echo "ERR=$PROJECT_ROOT/logs/Yb_normalized_coverage.$job.err"
echo "PDF=$PROJECT_ROOT/analysis/yb_normalized_coverage/results/Yb_normalized_coverage.pdf"
echo "SUMMARY=$PROJECT_ROOT/analysis/yb_normalized_coverage/results/Yb_normalized_coverage.summary.tsv"
echo
squeue -j "$job"
