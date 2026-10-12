#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_ROOT=${GYNANDRO_CODE_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo}
mkdir -p "$PROJECT_ROOT/logs" "$PROJECT_ROOT/analysis/yb_legacy_replication/rebuild_inventory"

job=$(sbatch --parsable \
  --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT" \
  "$CODE_ROOT/scripts/yb/run_Yb_existing_bams_audit.nano.slurm")
job=${job%%;*}
echo "JOBID=$job"
echo "OUT=$PROJECT_ROOT/logs/Yb_existing_BAM_audit.$job.out"
echo "ERR=$PROJECT_ROOT/logs/Yb_existing_BAM_audit.$job.err"
squeue -j "$job"
