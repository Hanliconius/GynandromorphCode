#!/usr/bin/env bash
#SBATCH --job-name=BwestHeur
#SBATCH --partition=nano
#SBATCH --time=00:29:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G

set -euo pipefail

project=/scratch/martinlab/joe/GynandroCalling
gyn_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python
gyn_rscript=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/Rscript

cd "${project}"

hostname
date --iso-8601=seconds

"${gyn_python}" \
    scripts/141_B_western_plus_heurippa_adjusted_LR.py

"${gyn_rscript}" \
    -e 'invisible(parse(file="scripts/141_plot_B_western_plus_heurippa_adjusted_LR.R")); cat("R syntax check passed.\n")'

"${gyn_rscript}" \
    scripts/141_plot_B_western_plus_heurippa_adjusted_LR.R

date --iso-8601=seconds
