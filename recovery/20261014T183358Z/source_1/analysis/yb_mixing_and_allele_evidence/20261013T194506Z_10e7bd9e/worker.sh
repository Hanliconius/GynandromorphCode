#!/usr/bin/env bash
set -euo pipefail
exec /gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_mixing_and_allele_evidence/20261013T194506Z_10e7bd9e/run_Yb_mixing_and_allele_evidence.py --config /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_mixing_and_allele_evidence/20261013T194506Z_10e7bd9e/run.json --stage "$1" --task "$SLURM_ARRAY_TASK_ID"
