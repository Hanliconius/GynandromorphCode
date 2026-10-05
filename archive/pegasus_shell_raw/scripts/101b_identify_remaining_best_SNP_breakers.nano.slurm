#!/usr/bin/env bash
#SBATCH --job-name=bestbreak
#SBATCH --partition=nano
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/101b_identify_remaining_best_SNP_breakers.nano.%j.out
#SBATCH --error=logs/101b_identify_remaining_best_SNP_breakers.nano.%j.err

project_dir="/scratch/martinlab/joe/GynandroCalling"
python_executable="/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python"

cd "${project_dir}"

"${python_executable}" \
    scripts/101b_identify_remaining_best_SNP_breakers.py
