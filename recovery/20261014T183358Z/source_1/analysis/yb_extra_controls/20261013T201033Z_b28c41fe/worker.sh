#!/usr/bin/env bash
set -euo pipefail
if [[ "$1" == download ]]; then
 exec /gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_extra_controls/20261013T201033Z_b28c41fe/submit_Yb_extra_controls.py --config /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_extra_controls/20261013T201033Z_b28c41fe/run.json --task "$SLURM_ARRAY_TASK_ID"
else
 exec /gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_extra_controls/20261013T201033Z_b28c41fe/rebuild_Yb_local_both_references.py --config /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_extra_controls/20261013T201033Z_b28c41fe/run.json --stage align --task "$SLURM_ARRAY_TASK_ID"
fi
