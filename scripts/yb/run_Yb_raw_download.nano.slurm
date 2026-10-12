#!/usr/bin/env bash
#SBATCH --job-name=YbRawDl
#SBATCH --partition=nano
#SBATCH --array=0-3%4
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_raw_download.%A_%a.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_raw_download.%A_%a.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${YB_CODE_DIR:?YB_CODE_DIR must be exported}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"
OUT="$PROJECT_ROOT/analysis/yb_raw_from_scratch/fastq"

[[ -x "$PYTHON" ]] || { echo "ERROR: Python not executable: $PYTHON" >&2; exit 2; }
mkdir -p "$OUT"

"$PYTHON" "$CODE_DIR/download_Yb_raw_fastq.py" \
    --manifest "$CODE_DIR/yb_raw_samples.tsv" \
    --task-id "$SLURM_ARRAY_TASK_ID" \
    --output-dir "$OUT"
