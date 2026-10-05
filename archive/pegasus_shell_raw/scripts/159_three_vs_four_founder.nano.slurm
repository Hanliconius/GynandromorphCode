#!/bin/bash
#SBATCH --job-name=founder_3v4
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=logs/159_three_vs_four_founder.%j.out
#SBATCH --error=logs/159_three_vs_four_founder.%j.err

set -u

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

working_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

date --iso-8601=seconds

"${working_python}" \
    scripts/159_fit_autosomal_three_vs_four_founder.py

run_status=$?

printf 'run_status=%s\n' "${run_status}"

if [[ "${run_status}" -eq 0 ]]
then
    printf '\nOBSERVED FITS\n'

    cat \
        analysis/autosomal_three_vs_four_founder/observed_model_fits.tsv

    printf '\nOUTPUT FILES\n'

    find \
        analysis/autosomal_three_vs_four_founder \
        -maxdepth 1 \
        -type f \
        -printf '%f\t%s bytes\n' \
        | sort
fi

date --iso-8601=seconds

test "${run_status}" -eq 0
