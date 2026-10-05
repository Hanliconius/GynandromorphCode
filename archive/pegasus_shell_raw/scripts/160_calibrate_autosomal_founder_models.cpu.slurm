#!/bin/bash
#SBATCH --job-name=gyn160cal
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/160_calibrate_autosomal_founder_models.%j.out
#SBATCH --error=logs/160_calibrate_autosomal_founder_models.%j.err

set -uo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
working_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

cd "${project}"
mkdir -p logs

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/gyn160_matplotlib_${SLURM_JOB_ID}"
mkdir -p "${MPLCONFIGDIR}"

run_log="logs/160_calibrate_autosomal_founder_models.${SLURM_JOB_ID}.combined.log"

"${working_python}" \
    scripts/160_simulate_calibrate_autosomal_founder_models.py \
    --replicates-per-model 196 \
    --bootstrap-replicates 200 \
    2>&1 \
    | tee "${run_log}"

run_status=${PIPESTATUS[0]}
output_status=0

for output in \
    analysis/autosomal_three_vs_four_founder/observed_nested_founder_fit.tsv \
    analysis/autosomal_three_vs_four_founder/founder_simulation_replicates.tsv.gz \
    analysis/autosomal_three_vs_four_founder/founder_simulation_summary.tsv \
    analysis/autosomal_three_vs_four_founder/mixing_calibration.tsv \
    analysis/autosomal_three_vs_four_founder/observed_10kb_block_bootstrap.tsv.gz \
    analysis/autosomal_three_vs_four_founder/calibrated_founder_and_mixing_summary.tsv \
    analysis/autosomal_three_vs_four_founder/simulation_all_47_controls.tsv \
    analysis/autosomal_three_vs_four_founder/three_vs_four_founder_simulation_calibration.pdf \
    analysis/autosomal_three_vs_four_founder/three_vs_four_founder_simulation_calibration.summary.png \
    analysis/autosomal_three_vs_four_founder/simulation_output_checksums.sha256
do
    if [[ -s "${output}" ]]
    then
        printf 'PRESENT\t%s\t%s bytes\n' \
            "${output}" \
            "$(wc -c < "${output}")"
    else
        printf 'MISSING_OR_EMPTY\t%s\n' "${output}"
        output_status=1
    fi
done

printf '\nFINAL STATUS\n'
printf 'run_status=%s\n' "${run_status}"
printf 'output_status=%s\n' "${output_status}"

if [[ "${run_status}" -ne 0 || "${output_status}" -ne 0 ]]
then
    false
fi
