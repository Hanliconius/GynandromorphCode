#!/usr/bin/env bash
#SBATCH --job-name=gynBestErr
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=/scratch/martinlab/joe/GynandroCalling/logs/101a_identify_best_B_D_R_site_breakers.nano.%j.out
#SBATCH --error=/scratch/martinlab/joe/GynandroCalling/logs/101a_identify_best_B_D_R_site_breakers.nano.%j.err

project_dir="/scratch/martinlab/joe/GynandroCalling"
env_python="/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python"

cd "${project_dir}"

if [[ -x "${env_python}" ]]
then
    "${env_python}" \
        scripts/101a_identify_best_B_D_R_site_breakers.py
else
    printf 'ERROR: environment Python is missing: %s\n' \
        "${env_python}" >&2
    false
fi
