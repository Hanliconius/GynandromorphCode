#!/usr/bin/env bash
#SBATCH --job-name=DRgenoplot
#SBATCH --partition=nano
#SBATCH --time=00:25:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=logs/104_105_make_perfect_D_R_all_comparative_genoplots.nano.%j.out
#SBATCH --error=logs/104_105_make_perfect_D_R_all_comparative_genoplots.nano.%j.err

project_dir="/scratch/martinlab/joe/GynandroCalling"
environment_bin="/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin"

cd "${project_dir}"

"${environment_bin}/python" \
    scripts/104_extract_perfect_D_R_all_comparative_genotypes.py

"${environment_bin}/Rscript" \
    scripts/105_plot_perfect_D_R_all_comparative_genotypes.R
