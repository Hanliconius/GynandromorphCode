#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

echo "===== Hmel2.5 Yb FIGURE-8 NORMALIZATION PREFLIGHT ====="

"$PYTHON" -m py_compile     "$HERE/plot_Yb_Hmel2_5_Figure8_normalized.py"

bash -n     "$HERE/run_Yb_Hmel2_5_Figure8_normalized.nano.slurm"

echo "REFERENCE=Hmel2.5 / H. m. melpomene"
echo "INPUT=Krista historical bigWigs"
echo "BIN_SIZE=50bp"
echo "NORMALIZATION=each sample / its mean depth across full Hmel215003o"
echo "IMPORTANT=historical bigWigs cannot retroactively enforce Figure-8 BAM MAPQ>=30 filtering"
echo "PREFLIGHT=PASS"
echo

job=$(
    sbatch         --parsable         --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,YB_CODE_DIR=$HERE"         "$HERE/run_Yb_Hmel2_5_Figure8_normalized.nano.slurm"
)
job=${job%%;*}

echo "JOBID=$job"
echo "OUT=$PROJECT_ROOT/logs/Yb_Hmel2_5_normalized.$job.out"
echo "ERR=$PROJECT_ROOT/logs/Yb_Hmel2_5_normalized.$job.err"
echo "PDF=$PROJECT_ROOT/analysis/yb_Hmel2_5_normalized_coverage/results/Yb_Hmel2.5_Figure8_style_normalized.pdf"
echo "SUMMARY=$PROJECT_ROOT/analysis/yb_Hmel2_5_normalized_coverage/results/Yb_Hmel2.5_Figure8_normalized.summary.tsv"
echo
squeue -j "$job"
