#!/bin/bash
#SBATCH --job-name=BvsDR161
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/161_B_vs_DR_repolarized.%j.out
#SBATCH --error=logs/161_B_vs_DR_repolarized.%j.err

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

gyn_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python
gyn_rscript=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/Rscript

printf '%s\n' \
    '===== ANALYSIS 161: RECOMPUTE B vs DR ====='

"${gyn_python}" \
    scripts/161_repolarize_optix_B_vs_DR.py

python_status=$?

printf '\npython_status=%s\n' \
    "${python_status}"

if [[ "${python_status}" -eq 0 ]]
then
    printf '\n%s\n' \
        '===== PARSE PLOTTER ====='

    "${gyn_rscript}" \
        -e 'invisible(parse(file="scripts/161_plot_B_vs_DR_repolarized.R")); cat("161 R plotter parses successfully.\n")'

    parse_status=$?

    if [[ "${parse_status}" -eq 0 ]]
    then
        printf '\n%s\n' \
            '===== PLOT ====='

        "${gyn_rscript}" \
            scripts/161_plot_B_vs_DR_repolarized.R

        plot_status=$?

        printf '\nparse_status=%s plot_status=%s\n' \
            "${parse_status}" \
            "${plot_status}"

        if [[ "${plot_status}" -eq 0 ]]
        then
            printf '\n%s\n' \
                '===== SUMMARY ====='

            column -t -s $'\t' \
                analysis/B_vs_DR_repolarization/checkpoints/B_vs_DR_repolarization_summary.tsv

            printf '\n%s\n' \
                '===== FOCAL PAIRS ====='

            column -t -s $'\t' \
                analysis/B_vs_DR_repolarization/checkpoints/B_vs_DR_focal_state_pair_summary.tsv

            printf '\n%s\n' \
                '===== OUTPUTS ====='

            ls -lh \
                analysis/figures/B_vs_DR_repolarization
        fi
    fi
fi
