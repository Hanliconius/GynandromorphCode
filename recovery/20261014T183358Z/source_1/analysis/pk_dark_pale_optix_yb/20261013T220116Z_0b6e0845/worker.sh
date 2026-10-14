#!/usr/bin/env bash
set -euo pipefail
exec /gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/pk_dark_pale_optix_yb/20261013T220116Z_0b6e0845/extend_PK_optix_Yb.py --config /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/pk_dark_pale_optix_yb/20261013T220116Z_0b6e0845/run.json --stage "$1" --task "${SLURM_ARRAY_TASK_ID:-0}"
