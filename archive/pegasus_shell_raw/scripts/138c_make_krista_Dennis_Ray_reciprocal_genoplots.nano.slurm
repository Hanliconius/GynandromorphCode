#!/usr/bin/env bash

#SBATCH --job-name=kristaDR
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=logs/138c_make_krista_Dennis_Ray_reciprocal_genoplots.nano.%j.out
#SBATCH --error=logs/138c_make_krista_Dennis_Ray_reciprocal_genoplots.nano.%j.err

project=/scratch/martinlab/joe/GynandroCalling
gyn_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python
gyn_rscript=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/Rscript

cd "${project}"

printf 'Host: %s\n' "$(hostname)"
printf 'Started: %s\n\n' "$(date --iso-8601=seconds)"

if "${gyn_python}" -m py_compile \
    scripts/138a_build_krista_Dennis_Ray_reciprocal_genoplots.py \
    &&
    "${gyn_rscript}" \
        -e 'invisible(parse(file="scripts/138b_plot_krista_Dennis_Ray_reciprocal_genoplots.R")); cat("R syntax check passed.\n")'
then
    if "${gyn_python}" \
        scripts/138a_build_krista_Dennis_Ray_reciprocal_genoplots.py
    then
        "${gyn_rscript}" \
            scripts/138b_plot_krista_Dennis_Ray_reciprocal_genoplots.R
    else
        printf '%s\n' \
            'Python table construction failed; plotting was not run.'
    fi
else
    printf '%s\n' \
        'Syntax checking failed; analysis was not run.'
fi

printf '\nFinished: %s\n' "$(date --iso-8601=seconds)"
