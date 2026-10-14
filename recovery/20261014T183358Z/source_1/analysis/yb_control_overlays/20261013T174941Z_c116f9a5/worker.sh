#!/usr/bin/env bash
set -euo pipefail
exec /gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_control_overlays/20261013T174941Z_c116f9a5/Yb_control_overlays.py --config /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_control_overlays/20261013T174941Z_c116f9a5/run.json --stage "$1" --task "$SLURM_ARRAY_TASK_ID"
