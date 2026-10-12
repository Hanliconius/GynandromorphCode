#!/usr/bin/env bash
#SBATCH --job-name=YbRawPlot
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_raw_plot.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_raw_plot.%j.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${YB_CODE_DIR:?YB_CODE_DIR must be exported}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"
INPUT="$PROJECT_ROOT/analysis/yb_raw_from_scratch/HmRos"
OUTPUT="$INPUT/plots"

"$PYTHON" "$CODE_DIR/plot_Yb_raw_Krista_style.py" \
    --input-dir "$INPUT" --output-dir "$OUTPUT"
