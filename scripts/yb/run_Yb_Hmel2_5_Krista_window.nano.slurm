#!/usr/bin/env bash
#SBATCH --job-name=YbKristaWindow
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_Krista_window.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_Krista_window.%j.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${YB_CODE_DIR:?YB_CODE_DIR must be exported}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"
BEDGRAPH_DIR="$PROJECT_ROOT/analysis/yb_Hmel2_5_normalized_coverage/Hmel215003o_bedgraph"
OUTPUT="$PROJECT_ROOT/analysis/yb_legacy_replication/Hmel2.5_Krista_window"

"$PYTHON" -m py_compile "$CODE_DIR/plot_Yb_Hmel2_5_Krista_window.py"
mkdir -p "$OUTPUT"

echo "YB_HMEL2_5_KRISTA_WINDOW_START"
echo "JOB_ID=$SLURM_JOB_ID"
echo "BEDGRAPH_DIR=$BEDGRAPH_DIR"
echo "OUTPUT=$OUTPUT"

"$PYTHON" "$CODE_DIR/plot_Yb_Hmel2_5_Krista_window.py" \
    --bedgraph-dir "$BEDGRAPH_DIR" \
    --output-dir "$OUTPUT"

echo "YB_HMEL2_5_KRISTA_WINDOW_FINISHED"
