#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"
INPUT="$PROJECT_ROOT/analysis/yb_hmros_normalized_coverage/results/Yb_hmros_panel_normalized_250bp.long.tsv.gz"
SAMPLES="$PROJECT_ROOT/analysis/yb_hmros_normalized_coverage/metadata/samples.tsv"
COORDINATES="$PROJECT_ROOT/analysis/yb_hmros_normalized_coverage/metadata/Yb_hmros_coordinates.tsv"
SAMTOOLS="$HOME/.local/share/mamba/envs/colias_w_popgen/bin/samtools"

[[ -s "$INPUT" ]] || {
    echo "ERROR: required normalized HmRos table is missing: $INPUT"
    exit 1
}
for f in "$SAMPLES" "$COORDINATES" "$SAMTOOLS"
do
    [[ -s "$f" ]] || {
        echo "ERROR: required HmRos input is missing: $f"
        exit 1
    }
done

"$PYTHON" -m py_compile "$HERE/plot_Yb_HmRos_legacy_replication.py"
bash -n "$HERE/run_Yb_HmRos_legacy_replication.nano.slurm"

job=$(
    sbatch --parsable \
        --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,YB_CODE_DIR=$HERE" \
        "$HERE/run_Yb_HmRos_legacy_replication.nano.slurm"
)
job=${job%%;*}

echo "JOBID=$job"
echo "OUT=$PROJECT_ROOT/logs/Yb_HmRos_legacy_plot.$job.out"
echo "ERR=$PROJECT_ROOT/logs/Yb_HmRos_legacy_plot.$job.err"
echo "PDF_DIR=$PROJECT_ROOT/analysis/yb_legacy_replication/HmRos"
squeue -j "$job"
