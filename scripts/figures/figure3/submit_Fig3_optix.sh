#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

mkdir -p "$PROJECT_ROOT/logs"

echo "===== FIGURE 3 OPTIX PREFLIGHT ====="
"$PYTHON" -m py_compile "$HERE/Fig3_optix_collapsed_pattern_dosage.py"
bash -n "$HERE/run_Fig3_optix.nano.slurm"
echo "PREFLIGHT=PASS"

job=$(
    sbatch         --parsable         --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,FIGURE3_CODE_DIR=$HERE"         "$HERE/run_Fig3_optix.nano.slurm"
)
job=${job%%;*}

echo "FIGURE3_OPTIX_JOB=$job"
echo "LOG_OUT=$PROJECT_ROOT/logs/Fig3_optix.$job.out"
echo "LOG_ERR=$PROJECT_ROOT/logs/Fig3_optix.$job.err"
echo "OUTPUT_DIR=$PROJECT_ROOT/analysis/figure3/optix_pattern_dosage"
squeue -j "$job"
