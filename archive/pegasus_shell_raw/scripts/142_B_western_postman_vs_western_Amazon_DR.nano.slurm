#!/usr/bin/env bash
#SBATCH --job-name=BwestDR
#SBATCH --partition=nano
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=logs/142_B_western_postman_vs_western_Amazon_DR.nano.%j.out
#SBATCH --error=logs/142_B_western_postman_vs_western_Amazon_DR.nano.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
gyn_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python
gyn_rscript=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/Rscript

cd "${project}"

hostname
date --iso-8601=seconds

"${gyn_python}" \
    scripts/142_B_western_postman_vs_western_Amazon_DR_adjusted_LR.py

"${gyn_rscript}" \
    -e 'invisible(parse(file="scripts/142_plot_B_western_postman_vs_western_Amazon_DR_adjusted_LR.R")); cat("R syntax check passed.\n")'

"${gyn_rscript}" \
    scripts/142_plot_B_western_postman_vs_western_Amazon_DR_adjusted_LR.R

date --iso-8601=seconds
