#!/bin/bash
#SBATCH --job-name=TN_AF40
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/TN_AF40.%j.out
#SBATCH --error=logs/TN_AF40.%j.err

set -euo pipefail
cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

pairwise-af/.venv-af/bin/python scripts/regenerate_TN_AF_hex40.py
