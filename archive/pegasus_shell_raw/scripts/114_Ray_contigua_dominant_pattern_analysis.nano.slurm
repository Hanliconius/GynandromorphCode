#!/usr/bin/env bash
#SBATCH --job-name=rayDom
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/114_Ray_contigua_dominant_pattern_analysis.nano.%j.out
#SBATCH --error=logs/114_Ray_contigua_dominant_pattern_analysis.nano.%j.err

set -u

project="/scratch/martinlab/joe/GynandroCalling"
python_bin="/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python"

cd "${project}"

printf 'Host: %s\n' "$(hostname)"
printf 'Started: %s\n' "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/114_Ray_contigua_dominant_pattern_analysis.py

printf 'Finished: %s\n' "$(date --iso-8601=seconds)"
