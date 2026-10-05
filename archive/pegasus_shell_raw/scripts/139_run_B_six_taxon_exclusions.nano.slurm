#!/bin/bash
#SBATCH --job-name=B6taxa
#SBATCH --partition=nano
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --chdir=/scratch/martinlab/joe/GynandroCalling
#SBATCH --output=logs/139_B_six_taxon_exclusions.nano.%j.out
#SBATCH --error=logs/139_B_six_taxon_exclusions.nano.%j.err
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

gyn_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python
gyn_rscript=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/Rscript

hostname
date -Is

"${gyn_python}" -u scripts/139_recalculate_B_six_taxon_exclusions.py

"${gyn_rscript}" \
    -e 'invisible(parse(file="scripts/139_plot_B_six_taxon_exclusions.R"))'

"${gyn_rscript}" scripts/139_plot_B_six_taxon_exclusions.R

date -Is
