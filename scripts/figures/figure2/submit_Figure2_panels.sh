#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}

mkdir -p "${PROJECT_ROOT}/logs" "${PROJECT_ROOT}/analysis/figure2/panels"

job=$(sbatch --parsable     --output="${PROJECT_ROOT}/logs/Fig2_panel.%A_%a.out"     --error="${PROJECT_ROOT}/logs/Fig2_panel.%A_%a.err"     --export="ALL,GYNANDRO_ROOT=${PROJECT_ROOT}"     "${HERE}/run_Figure2_panels.nano.slurm")

job=${job%%;*}
echo "FIGURE2_ARRAY_JOB=${job}"
echo "squeue -j ${job}"
