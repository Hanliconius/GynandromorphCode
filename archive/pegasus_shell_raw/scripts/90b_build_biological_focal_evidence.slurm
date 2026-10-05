#!/bin/bash -l
#SBATCH --job-name=gyn90b
#SBATCH --partition=nano
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/90b_build_biological_focal_evidence.%j.out
#SBATCH --error=logs/90b_build_biological_focal_evidence.%j.err

project_root="/scratch/martinlab/joe/GynandroCalling"

if [[ -d "${project_root}" ]]
then
    cd "${project_root}"

    micromamba run \
        -n gynandrocalling \
        python scripts/90b_build_biological_focal_evidence.py
else
    printf 'ERROR: project directory not found: %s\n' "${project_root}" >&2
fi
