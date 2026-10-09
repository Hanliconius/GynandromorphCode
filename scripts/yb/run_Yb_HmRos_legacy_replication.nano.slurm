#!/usr/bin/env bash
#SBATCH --job-name=YbHmRosPlot
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_HmRos_legacy_plot.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_HmRos_legacy_plot.%j.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${YB_CODE_DIR:?YB_CODE_DIR must be exported}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

INPUT="$PROJECT_ROOT/analysis/yb_hmros_normalized_coverage/results/Yb_hmros_panel_normalized_250bp.long.tsv.gz"
OUTPUT="$PROJECT_ROOT/analysis/yb_legacy_replication/HmRos"

[[ -s "$INPUT" ]] || {
    echo "ERROR: normalized HmRos table is missing: $INPUT" >&2
    exit 1
}

"$PYTHON" -m py_compile "$CODE_DIR/plot_Yb_HmRos_legacy_replication.py"
mkdir -p "$OUTPUT"

echo "YB_HMROS_LEGACY_PLOT_START"
echo "JOB_ID=$SLURM_JOB_ID"
echo "INPUT=$INPUT"
echo "OUTPUT=$OUTPUT"

"$PYTHON" "$CODE_DIR/plot_Yb_HmRos_legacy_replication.py" \
    --input "$INPUT" \
    --output-dir "$OUTPUT"

echo "YB_HMROS_LEGACY_PLOT_FINISHED"
