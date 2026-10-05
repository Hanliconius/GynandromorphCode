#!/usr/bin/env bash
set -u
set -o pipefail

main() {
    local project_root
    local env_name

    project_root="/scratch/martinlab/joe/GynandroCalling"
    env_name="gynandrocalling"

    mkdir -p \
        "${project_root}/scripts" \
        "${project_root}/logs" \
        "${project_root}/data/raw" \
        "${project_root}/data/derived" \
        "${project_root}/metadata" \
        "${project_root}/analysis/tables" \
        "${project_root}/analysis/figures" \
        "${project_root}/analysis/checkpoints"

    if ! command -v micromamba >/dev/null 2>&1; then
        printf 'ERROR: micromamba is not available in this shell.\n' >&2
        return 1
    fi

    if micromamba run -n "${env_name}" python --version >/dev/null 2>&1; then
        printf 'Environment %s already exists; updating it from environment.yml.\n' \
            "${env_name}"

        micromamba env update \
            -n "${env_name}" \
            --file "${project_root}/environment.yml" \
            --yes
    else
        printf 'Creating environment %s.\n' "${env_name}"

        micromamba create \
            -n "${env_name}" \
            --file "${project_root}/environment.yml" \
            --yes
    fi

    if micromamba run -n "${env_name}" \
        python -c 'import numpy, pandas, scipy, matplotlib, seaborn, openpyxl'
    then
        printf 'Python package check passed.\n'
    else
        printf 'ERROR: Python package check failed.\n' >&2
        return 1
    fi

    if micromamba run -n "${env_name}" \
        Rscript -e 'library(ggplot2); library(dplyr); library(readr); library(patchwork); library(ragg); library(svglite)'
    then
        printf 'R package check passed.\n'
    else
        printf 'ERROR: R package check failed.\n' >&2
        return 1
    fi

    printf '\nInstalled software versions\n'
    printf '%s\n' '---------------------------'

    micromamba run -n "${env_name}" python --version
    micromamba run -n "${env_name}" bcftools --version | head -n 1
    micromamba run -n "${env_name}" Rscript --version

    micromamba env export \
        -n "${env_name}" \
        > "${project_root}/analysis/checkpoints/environment.resolved.yml"

    printf '\nEnvironment setup completed.\n'
    printf 'Project:     %s\n' "${project_root}"
    printf 'Environment: %s\n' "${env_name}"
    printf 'Activate with: micromamba activate %s\n' "${env_name}"
}

main "$@"
