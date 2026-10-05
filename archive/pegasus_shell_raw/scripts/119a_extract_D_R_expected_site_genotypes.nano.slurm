#!/bin/bash
#SBATCH --job-name=DRgeno
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=logs/119a_extract_D_R_expected_site_genotypes.nano.%j.out
#SBATCH --error=logs/119a_extract_D_R_expected_site_genotypes.nano.%j.err

cd /scratch/martinlab/joe/GynandroCalling

gyn_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python

"${gyn_python}" \
    scripts/119a_extract_D_R_expected_site_genotypes.py
