#!/usr/bin/env bash
#SBATCH --job-name=BnoDRpos
#SBATCH --partition=nano
#SBATCH --time=00:25:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/110c_scan_B_no_DR_positive.nano.%j.out
#SBATCH --error=logs/110c_scan_B_no_DR_positive.nano.%j.err

project_dir="/scratch/martinlab/joe/GynandroCalling"
python_bin="/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python"

cd "${project_dir}"

"${python_bin}" \
    scripts/110c_scan_B_no_DR_positive.py
