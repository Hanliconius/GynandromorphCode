#!/usr/bin/env bash

project_dir="/scratch/martinlab/joe/GynandroCalling"
vcf_dir="${project_dir}/data/raw/shared_drive"
source_dir="${project_dir}/external/krista_analysis"
checkpoint_dir="${project_dir}/analysis/checkpoints"

mkdir -p "${vcf_dir}" "${source_dir}" "${checkpoint_dir}"

download_drive_file() {
    file_id="$1"
    output_file="$2"

    if [[ -s "${output_file}" ]]
    then
        printf 'Reusing: %s\n' "${output_file}"
        return
    fi

    temporary_file="${output_file}.part"

    printf 'Downloading: %s\n' "${output_file}"

    if gdown "${file_id}" -O "${temporary_file}" \
        && [[ -s "${temporary_file}" ]]
    then
        mv "${temporary_file}" "${output_file}"
        printf 'Downloaded: %s\n' "${output_file}"
    else
        printf 'ERROR: download failed: %s\n' "${output_file}" >&2
        rm -f "${temporary_file}"
    fi
}

# Missing Ivory VCFs
download_drive_file \
    "1Is-M_RFw19qQeJ2h0gL9gIm_hLfbKlvS" \
    "${vcf_dir}/Ivory_Hmr15-1627000-1637000_biallelic_gynandromorph.vcf"

download_drive_file \
    "1uEU7K0BwTg_cqwav0TGZBMxKOQ1esrX6" \
    "${vcf_dir}/Ivory_Hmr15-1627000-1637000_biallelic_no_miss.vcf"

download_drive_file \
    "1WolFQWMe7aEkrfrqIJosCChPROpQHCLm" \
    "${vcf_dir}/Ivory_biallelic.vcf.gz"

download_drive_file \
    "1Iaa0Eww-zQgV7TGRxaLMwkbQBZRyC44Q" \
    "${vcf_dir}/Ivory_biallelic.vcf.gz.csi"

download_drive_file \
    "1PriB6ACI9hJwiIx9VwXaShYL6xWqD7PZ" \
    "${vcf_dir}/Ivory_Hmr15-1492477-1794947.vcf"

# Original analysis source and rendered reports
download_drive_file \
    "13Yg1mfO_KjfddCBuLM3hfLBu4bURuZD6" \
    "${source_dir}/Optix_and_Ivory_Analysis.Rmd"

download_drive_file \
    "1bv01JVqziCPGYjVuBg0lyHLUf3jZ1eDU" \
    "${source_dir}/Optix_and_Ivory_Analysis.html"

download_drive_file \
    "1QxUpJ63b8VHNZKnW-_Ig_JduFFAGTuW4" \
    "${source_dir}/gynandromorph_genotype_plots.html"

download_drive_file \
    "10lfNq4SQzbAwaurPph3dUQTqbqatF5oT" \
    "${source_dir}/optix_band_genotype_plot-2.pdf"

download_drive_file \
    "1cpAgowLMXScBqltIphDdaJ__uLKDi_Zq" \
    "${source_dir}/ivory_genotype_plot-2.pdf"

download_drive_file \
    "1KOfMnFZYO3tDwllA3Cs3vGis_qTMS5Jt" \
    "${source_dir}/metadata.csv"

inventory="${checkpoint_dir}/missing_shared_analysis_files.tsv"

{
    printf 'category\tfile\tbytes\n'

    for file in "${vcf_dir}"/* "${source_dir}"/*
    do
        if [[ -f "${file}" ]]
        then
            case "${file}" in
                "${vcf_dir}"/*)
                    category="vcf_input"
                    ;;
                *)
                    category="analysis_source"
                    ;;
            esac

            printf '%s\t%s\t%s\n' \
                "${category}" \
                "$(basename "${file}")" \
                "$(stat -c '%s' "${file}")"
        fi
    done
} > "${inventory}"

printf '\nDOWNLOADED FILE INVENTORY\n'
column -t -s $'\t' "${inventory}"

printf '\nInventory written to:\n%s\n' "${inventory}"
