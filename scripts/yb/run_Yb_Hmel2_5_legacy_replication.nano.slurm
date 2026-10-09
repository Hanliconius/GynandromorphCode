#!/usr/bin/env bash
#SBATCH --job-name=YbLegacyPlot
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_legacy_plot.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_legacy_plot.%j.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${YB_CODE_DIR:?YB_CODE_DIR must be exported}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

INPUT="$PROJECT_ROOT/analysis/yb_Hmel2_5_normalized_coverage/results/Yb_Hmel2.5_Figure8_normalized_50bp.long.tsv.gz"
OUTPUT="$PROJECT_ROOT/analysis/yb_legacy_replication/Hmel2.5"

[[ -s "$INPUT" ]] || {
    echo "ERROR: normalized Hmel2.5 table is missing: $INPUT" >&2
    exit 1
}

"$PYTHON" -m py_compile "$CODE_DIR/plot_Yb_Hmel2_5_legacy_replication.py"
mkdir -p "$OUTPUT"

echo "YB_HMEL2_5_LEGACY_PLOT_START"
echo "JOB_ID=$SLURM_JOB_ID"
echo "INPUT=$INPUT"
echo "OUTPUT=$OUTPUT"

"$PYTHON" "$CODE_DIR/plot_Yb_Hmel2_5_legacy_replication.py" \
    --input "$INPUT" \
    --output-dir "$OUTPUT"

echo "YB_HMEL2_5_LEGACY_PLOT_FINISHED"
