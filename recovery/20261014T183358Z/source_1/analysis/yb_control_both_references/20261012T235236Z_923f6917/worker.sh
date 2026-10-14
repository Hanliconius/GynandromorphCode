#!/usr/bin/env bash
set -euo pipefail
exec /gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_control_both_references/20261012T235236Z_923f6917/rebuild_Yb_local_both_references.py --config /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_control_both_references/20261012T235236Z_923f6917/run.json --stage align --task "$SLURM_ARRAY_TASK_ID"
