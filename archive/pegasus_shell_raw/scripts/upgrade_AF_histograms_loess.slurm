#!/bin/bash
#SBATCH --job-name=AFhistLOESS
#SBATCH --partition=cpu
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/AFhistLOESS.%A_%a.out
#SBATCH --error=logs/AFhistLOESS.%A_%a.err

set -euo pipefail
cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

pairwise-af/.venv-af/bin/python \
    scripts/upgrade_AF_histograms_loess.py \
    extract --task "${SLURM_ARRAY_TASK_ID}"
