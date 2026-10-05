#!/bin/bash
#SBATCH --job-name=R_BA90
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=logs/121a_extract_all_R_BA90_genotypes.nano.%j.out
#SBATCH --error=logs/121a_extract_all_R_BA90_genotypes.nano.%j.err

cd /scratch/martinlab/joe/GynandroCalling

gyn_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python

"${gyn_python}" \
    scripts/121a_extract_all_R_BA90_genotypes.py
