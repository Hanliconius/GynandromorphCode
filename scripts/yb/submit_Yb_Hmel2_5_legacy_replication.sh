#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"
INPUT="$PROJECT_ROOT/analysis/yb_Hmel2_5_normalized_coverage/results/Yb_Hmel2.5_Figure8_normalized_50bp.long.tsv.gz"

[[ -s "$INPUT" ]] || {
    echo "ERROR: required normalized Hmel2.5 table is missing: $INPUT" >&2
    exit 1
}

"$PYTHON" -m py_compile "$HERE/plot_Yb_Hmel2_5_legacy_replication.py"
bash -n "$HERE/run_Yb_Hmel2_5_legacy_replication.nano.slurm"

job=$(
    sbatch --parsable \
        --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,YB_CODE_DIR=$HERE" \
        "$HERE/run_Yb_Hmel2_5_legacy_replication.nano.slurm"
)
job=${job%%;*}

echo "JOBID=$job"
echo "OUT=$PROJECT_ROOT/logs/Yb_legacy_plot.$job.out"
echo "ERR=$PROJECT_ROOT/logs/Yb_legacy_plot.$job.err"
echo "PDF_DIR=$PROJECT_ROOT/analysis/yb_legacy_replication/Hmel2.5"
squeue -j "$job"
