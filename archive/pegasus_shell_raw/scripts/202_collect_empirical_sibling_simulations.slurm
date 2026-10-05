#!/usr/bin/env bash
#SBATCH --job-name=gyn202sib
#SBATCH --partition=nano
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --chdir=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
#SBATCH --output=logs/202_collect_empirical_sibling_simulations.%j.out
#SBATCH --error=logs/202_collect_empirical_sibling_simulations.%j.err

set -euo pipefail

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

"${python_bin}" \
    scripts/202_204_finish_empirical_sibling_calibration.py \
    collect
