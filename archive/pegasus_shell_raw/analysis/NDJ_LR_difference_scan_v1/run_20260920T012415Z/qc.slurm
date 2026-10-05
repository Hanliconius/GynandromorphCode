#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/NDJ_LR_difference_scan_v1/run_20260920T012415Z/mplcache
/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/pairwise-af/.venv-af/bin/python /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/NDJ_LR_difference_scan_v1/run_20260920T012415Z/code/run_ndj_scan.py worker --run /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/NDJ_LR_difference_scan_v1/run_20260920T012415Z --stage qc --index "${SLURM_ARRAY_TASK_ID}"
