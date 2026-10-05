#!/usr/bin/env bash

project_dir="/scratch/martinlab/joe/GynandroCalling"
repository_dir="${project_dir}/external/gynandromorph_2024"
inventory="${project_dir}/analysis/checkpoints/original_plotting_code_inventory.txt"

if [[ -d "${repository_dir}/.git" ]]
then
    {
        printf 'Repository: %s\n' "${repository_dir}"
        printf 'Commit: '
        git -C "${repository_dir}" rev-parse HEAD

        printf '\nFILES MATCHING GENOTYPE-PLOT TERMS\n'

        git -C "${repository_dir}" grep \
            -I -l -i -E \
            -e 'genotype.?plot|allele.?proportion|allele.?balance|allelic.?depth|posterior|optix|ivory|yellow.?band|Phi|MCC|optix_filtered|ivory_filtered' \
            -- \
            '*.R' '*.Rmd' '*.py' '*.sh' '*.md' '*.txt' '*.ipynb' \
            2>/dev/null \
            | sort \
            | head -150

        printf '\nALL VERSIONED CODE AND DOCUMENTATION FILES\n'

        git -C "${repository_dir}" ls-files \
            | grep -E '\.(R|Rmd|py|sh|md|txt|ipynb)$' \
            | sort \
            | head -250

        printf '\nFILES WITH OPTIX, IVORY, VCF, OR GENOTYPE IN NAME\n'

        git -C "${repository_dir}" ls-files \
            | grep -Ei 'optix|ivory|vcf|genotyp' \
            | sort \
            | head -250

    } > "${inventory}"

    cat "${inventory}"
else
    printf 'ERROR: repository is absent: %s\n' \
        "${repository_dir}" >&2
fi
