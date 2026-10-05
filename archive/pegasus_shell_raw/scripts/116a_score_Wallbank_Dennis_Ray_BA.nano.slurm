#!/usr/bin/env bash
#SBATCH --job-name=wallDRBA
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=logs/116a_score_Wallbank_Dennis_Ray_BA.nano.%j.out
#SBATCH --error=logs/116a_score_Wallbank_Dennis_Ray_BA.nano.%j.err

set -euo pipefail

project_dir=/scratch/martinlab/joe/GynandroCalling
gyn_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python

cd "${project_dir}"

printf 'Host: %s\n' "$(hostname)"
printf 'Started: %s\n' "$(date --iso-8601=seconds)"

"${gyn_python}" \
    scripts/116a_score_Wallbank_Dennis_Ray_BA.py

printf 'Finished: %s\n' "$(date --iso-8601=seconds)"
