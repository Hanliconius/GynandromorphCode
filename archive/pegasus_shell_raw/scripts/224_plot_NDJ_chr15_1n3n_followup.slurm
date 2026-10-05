#!/bin/bash
#SBATCH --job-name=NDJ15fit
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/224_plot_NDJ_chr15_1n3n_followup.%j.out
#SBATCH --error=logs/224_plot_NDJ_chr15_1n3n_followup.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=${project}/pairwise-af/.venv-af/bin/python
run=${project}/analysis/NDJ_LR_difference_scan_v1/run_20260920T012415Z

cd "${project}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPLCONFIGDIR="${run}/mplcache_224"

printf 'STAGE_START\tstage=NDJ_chr15_1n3n_followup\tjob=%s\thost=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/224_plot_NDJ_chr15_1n3n_followup.py \
    --run "${run}"

printf 'STAGE_END\tstage=NDJ_chr15_1n3n_followup\tstatus=0\ttime=%s\n' \
    "$(date --iso-8601=seconds)"
