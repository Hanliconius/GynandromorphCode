#!/usr/bin/env bash
#SBATCH --job-name=gyn210rec
#SBATCH --partition=nano
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/210_plot_reciprocal_mixture_recaller_narrative.%j.out
#SBATCH --error=logs/210_plot_reciprocal_mixture_recaller_narrative.%j.err

set -euo pipefail

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

printf 'STAGE_START\tstage=reciprocal_mixture_recaller_narrative\tjob=%s\thost=%s\ttime=%s\n' \
    "${SLURM_JOB_ID}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/210_plot_reciprocal_mixture_recaller_narrative.py

pdf=analysis/genomewide_mixing_corrected_v1/figures/reciprocal_mixture_snp_recaller_narrative.pdf
audit=analysis/genomewide_mixing_corrected_v1/figures/reciprocal_mixture_snp_recaller_narrative.audit.json

for output in "${pdf}" "${audit}"
do
    if [[ ! -s "${output}" ]]
    then
        printf 'ERROR: missing or empty output: %s\n' "${output}" >&2
        exit 1
    fi
done

printf 'STAGE_END\tstatus=0\tpdf=%s\ttime=%s\n' \
    "$(realpath "${pdf}")" \
    "$(date --iso-8601=seconds)"
