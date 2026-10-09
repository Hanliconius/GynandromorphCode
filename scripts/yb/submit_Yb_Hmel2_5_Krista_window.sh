#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"
BEDGRAPH_DIR="$PROJECT_ROOT/analysis/yb_Hmel2_5_normalized_coverage/Hmel215003o_bedgraph"

"$PYTHON" -m py_compile "$HERE/plot_Yb_Hmel2_5_Krista_window.py"
bash -n "$HERE/run_Yb_Hmel2_5_Krista_window.nano.slurm"

for sample in TN2L_S1 TN2R_S2 TN3L_S3 TN3R_S4
do
    file="$BEDGRAPH_DIR/${sample}.Hmel215003o.bedGraph"
    [[ -s "$file" ]] || {
        echo "ERROR: required historical bedGraph is missing: $file"
        exit 1
    }
done

job=$(
    sbatch --parsable \
        --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,YB_CODE_DIR=$HERE" \
        "$HERE/run_Yb_Hmel2_5_Krista_window.nano.slurm"
)
job=${job%%;*}

echo "JOBID=$job"
echo "OUT=$PROJECT_ROOT/logs/Yb_Krista_window.$job.out"
echo "ERR=$PROJECT_ROOT/logs/Yb_Krista_window.$job.err"
echo "PDF_DIR=$PROJECT_ROOT/analysis/yb_legacy_replication/Hmel2.5_Krista_window"
squeue -j "$job"
