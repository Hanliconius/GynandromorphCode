#!/bin/bash
#SBATCH --job-name=optixAF
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=logs/132_extract_full_optix_AF.nano.%j.out
#SBATCH --error=logs/132_extract_full_optix_AF.nano.%j.err

cd /scratch/martinlab/joe/GynandroCalling

gyn_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python

"${gyn_python}" \
    scripts/132_extract_full_optix_AF_focal_and_WGS_controls.py
