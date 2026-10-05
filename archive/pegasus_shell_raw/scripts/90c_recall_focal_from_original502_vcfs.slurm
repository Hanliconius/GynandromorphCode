#!/bin/bash -l
#SBATCH --job-name=gyn90c
#SBATCH --partition=nano
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/90c_recall_focal_from_original502_vcfs.%j.out
#SBATCH --error=logs/90c_recall_focal_from_original502_vcfs.%j.err

project_root="/scratch/martinlab/joe/GynandroCalling"

if [[ -d "${project_root}" ]]
then
    cd "${project_root}"

    micromamba run \
        -n gynandrocalling \
        python scripts/90c_recall_focal_from_original502_vcfs.py
else
    printf 'ERROR: project directory not found: %s\n' \
        "${project_root}" >&2
fi
