#!/usr/bin/env bash
#SBATCH --job-name=focalDR
#SBATCH --partition=nano
#SBATCH --time=00:25:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --output=logs/106a_extract_focal_rows_for_perfect_D_R_genoplots.nano.%j.out
#SBATCH --error=logs/106a_extract_focal_rows_for_perfect_D_R_genoplots.nano.%j.err

project_dir="/scratch/martinlab/joe/GynandroCalling"
python_bin="/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python"

cd "${project_dir}"

"${python_bin}" \
    scripts/106a_extract_focal_rows_for_perfect_D_R_genoplots.py
