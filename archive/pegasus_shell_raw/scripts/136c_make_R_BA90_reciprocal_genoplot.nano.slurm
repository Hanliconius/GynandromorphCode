#!/bin/bash
#SBATCH --job-name=RmixGeno
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/136c_make_R_BA90_reciprocal_genoplot.nano.%j.out
#SBATCH --error=logs/136c_make_R_BA90_reciprocal_genoplot.nano.%j.err

project=/scratch/martinlab/joe/GynandroCalling

gyn_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python
gyn_rscript=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/Rscript

cd "${project}"

printf 'Host: %s\n' "$(hostname)"
printf 'Started: %s\n\n' "$(date --iso-8601=seconds)"

if "${gyn_python}" \
    scripts/136a_build_R_BA90_reciprocal_genoplot_table.py
then
    printf '\n%s\n' 'Genotype table completed.'

    if "${gyn_rscript}" \
        scripts/136b_plot_R_BA90_reciprocal_genoplot.R
    then
        printf '\nFinished: %s\n' "$(date --iso-8601=seconds)"
    else
        printf '%s\n' 'ERROR: genotype table succeeded but plotting failed.'
    fi
else
    printf '%s\n' 'ERROR: genotype-table preparation failed; plotting was not attempted.'
fi
