#!/usr/bin/env bash
#SBATCH --job-name=gyn203sib
#SBATCH --partition=cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --chdir=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
#SBATCH --output=logs/203_calculate_empirical_sibling_window_evidence.%j.out
#SBATCH --error=logs/203_calculate_empirical_sibling_window_evidence.%j.err

set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

"${python_bin}" \
    scripts/202_204_finish_empirical_sibling_calibration.py \
    evidence
