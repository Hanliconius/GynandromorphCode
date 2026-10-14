#!/usr/bin/env bash
set -euo pipefail
RUN=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_local_both_references/20261012T194933Z_7fa9de96
PYTHON=/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python
PLOT=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_local_both_references/20261012T194933Z_7fa9de96/plot_retry_20261012T232027Z_2543983/plot_Yb_local_both_references.py
refs=(Hmel2.5 HmRos)
ref=${refs[$SLURM_ARRAY_TASK_ID]}
exec "$PYTHON" "$PLOT" --reference "$ref" --input-dir "$RUN/$ref/coverage" --output-dir "$RUN/$ref/plots"
