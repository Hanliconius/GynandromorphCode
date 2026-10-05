#!/usr/bin/env bash
#SBATCH --job-name=gyn204sib
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --chdir=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
#SBATCH --output=logs/204_plot_empirical_sibling_narrative.%j.out
#SBATCH --error=logs/204_plot_empirical_sibling_narrative.%j.err

set -euo pipefail

export MPLBACKEND=Agg
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

"${python_bin}" \
    scripts/202_204_finish_empirical_sibling_calibration.py \
    plot
