#!/usr/bin/env bash
#SBATCH --job-name=gyn207nar
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=logs/207_build_sequential_founder_narrative.%j.out
#SBATCH --error=logs/207_build_sequential_founder_narrative.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

cd "${project}"

printf 'STAGE_START\tstage=build_sequential_founder_narrative\tjob=%s\thost=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/207_build_sequential_founder_narrative.py

printf 'STAGE_END\tstage=build_sequential_founder_narrative\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
