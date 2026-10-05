#!/usr/bin/env bash
#SBATCH --job-name=fixedSNP
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/112_replicate_wallbank2016_fixed_SNP_analysis.nano.%j.out
#SBATCH --error=logs/112_replicate_wallbank2016_fixed_SNP_analysis.nano.%j.err

set -u

project="/scratch/martinlab/joe/GynandroCalling"
python_bin="/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python"

cd "${project}"

printf 'Host: %s\n' "$(hostname)"
printf 'Started: %s\n' "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/112_replicate_wallbank2016_fixed_SNP_analysis.py

printf 'Finished: %s\n' "$(date --iso-8601=seconds)"
