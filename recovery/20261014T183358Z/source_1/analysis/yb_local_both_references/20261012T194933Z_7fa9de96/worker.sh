#!/usr/bin/env bash
set -euo pipefail
exec /gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_local_both_references/20261012T194933Z_7fa9de96/rebuild_Yb_local_both_references.py --config /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_local_both_references/20261012T194933Z_7fa9de96/run.json --stage "$1" --task "$SLURM_ARRAY_TASK_ID"
