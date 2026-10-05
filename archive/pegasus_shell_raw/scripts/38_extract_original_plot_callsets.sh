#!/usr/bin/env bash

project_dir="/scratch/martinlab/joe/GynandroCalling"
archive="${project_dir}/Genotype_Plots.tar"
output_dir="${project_dir}/external/krista_genotype_plots_original"

mkdir -p \
    "${output_dir}" \
    "${project_dir}/logs" \
    "${project_dir}/analysis/checkpoints"

members=(
    "Genotype_Plots/optix_band-1112651-1114314.vcf"
    "Genotype_Plots/dennis1-1215210-1216124.vcf"
    "Genotype_Plots/dennis2-1218557-1221494.vcf"
    "Genotype_Plots/ray1-1156746-1160220.vcf"
    "Genotype_Plots/ray2-1174547-1187833.vcf"
    "Genotype_Plots/Ivory_Hmr15-1627000-1637000_biallelic.vcf"
    "Genotype_Plots/Optix_metadata.csv"
    "Genotype_Plots/Ivory_metadata.csv"
    "Genotype_Plots/gynandromorph_genotype_plots.Rmd"
    "Genotype_Plots/Outputs"
)

missing_members=()

if [[ ! -s "${archive}" ]]
then
    printf 'ERROR: archive is missing or empty: %s\n' "${archive}"
else
    for member in "${members[@]}"
    do
        target="${output_dir}/${member}"

        if [[ -f "${target}" && -s "${target}" ]]
        then
            printf 'Reusing: %s\n' "${target}"
        elif [[ -d "${target}" ]]
        then
            printf 'Reusing directory: %s\n' "${target}"
        else
            missing_members+=("${member}")
        fi
    done

    if (( ${#missing_members[@]} > 0 ))
    then
        printf 'Extracting %s selected archive members.\n' \
            "${#missing_members[@]}"

        tar \
            --extract \
            --file "${archive}" \
            --directory "${output_dir}" \
            --no-same-owner \
            --no-same-permissions \
            "${missing_members[@]}"

        tar_status=$?

        if (( tar_status == 0 ))
        then
            printf 'Selected extraction completed.\n'
        else
            printf \
                'ERROR: tar returned status %s; inspect extracted files.\n' \
                "${tar_status}"
        fi
    else
        printf 'All selected members already exist.\n'
    fi

    printf '\nExtracted regional VCFs:\n'

    find "${output_dir}/Genotype_Plots" \
        -maxdepth 1 \
        -type f \
        -name '*.vcf' \
        ! -name 'Optix_Hmr18-938195-1239421_biallelic.vcf' \
        -printf '%s\t%p\n' \
        | sort -k2,2

    printf '\nThe 1.997 GB broad Optix VCF was not extracted.\n'
fi
