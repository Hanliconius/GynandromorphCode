#!/bin/bash -l
#SBATCH --job-name=gyn90d
#SBATCH --partition=nano
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/90d_audit_original502_focal_vcf_sources.%j.out
#SBATCH --error=logs/90d_audit_original502_focal_vcf_sources.%j.err

project_root="/scratch/martinlab/joe/GynandroCalling"

if [[ -d "${project_root}" ]]
then
    cd "${project_root}"

    micromamba run \
        -n gynandrocalling \
        python scripts/90d_audit_original502_focal_vcf_sources.py
else
    printf 'ERROR: project directory not found: %s\n' \
        "${project_root}" >&2
fi
