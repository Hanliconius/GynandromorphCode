#!/usr/bin/env bash
#SBATCH --job-name=gynBdrop6
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=/scratch/martinlab/joe/GynandroCalling/logs/100c_scan_band_drop_problematic.nano.%j.out
#SBATCH --error=/scratch/martinlab/joe/GynandroCalling/logs/100c_scan_band_drop_problematic.nano.%j.err

project_dir="/scratch/martinlab/joe/GynandroCalling"
env_python="/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python"

cd "${project_dir}"

if [[ -x "${env_python}" ]]
then
    "${env_python}" \
        scripts/100c_scan_band_drop_problematic.py
else
    printf 'ERROR: environment Python is missing: %s\n' \
        "${env_python}" >&2
    false
fi
