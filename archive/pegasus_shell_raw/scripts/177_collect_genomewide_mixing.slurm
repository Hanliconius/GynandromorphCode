#!/usr/bin/env bash
#SBATCH --job-name=gyn177cl
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=logs/177_collect_genomewide_mixing.%j.out
#SBATCH --error=logs/177_collect_genomewide_mixing.%j.err

set -euo pipefail

stage_name=collect_genomewide_mixing
start_seconds=${SECONDS}
printf 'STAGE_START\tstage=%s\tjob=%s\thost=%s\ttime=%s\n' \
    "${stage_name}" "${SLURM_JOB_ID:-NA}" "$(hostname)" "$(date --iso-8601=seconds)"
trap 'status=$?; printf "STAGE_END\tstage=%s\tstatus=%s\telapsed_seconds=%s\thost=%s\ttime=%s\n" "${stage_name}" "${status}" "$((SECONDS - start_seconds))" "$(hostname)" "$(date --iso-8601=seconds)"' EXIT

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

"${python_bin}" \
    "${project}/scripts/173_genomewide_mixing_models.py" \
    collect
