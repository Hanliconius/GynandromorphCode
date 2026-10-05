#!/usr/bin/env bash
#SBATCH --job-name=optixDAF
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=logs/117a_score_optix_deltaAF.nano.%j.out
#SBATCH --error=logs/117a_score_optix_deltaAF.nano.%j.err

set -euo pipefail

project_dir=/scratch/martinlab/joe/GynandroCalling
gyn_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python

cd "${project_dir}"

printf 'Host: %s\n' "$(hostname)"
printf 'Started: %s\n' "$(date --iso-8601=seconds)"

"${gyn_python}" \
    scripts/117a_score_optix_comparative_deltaAF_focal_dosage.py

printf 'Finished: %s\n' "$(date --iso-8601=seconds)"
