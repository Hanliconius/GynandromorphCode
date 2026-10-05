#!/usr/bin/env bash
set -u
set -o pipefail

main() {
    local project_root
    local raw_root
    local checkpoint_dir
    local sample_dir
    local summary
    local ranges
    local warning_log
    local focal_expected
    local file_count
    local audit_failed
    local vcf
    local relative
    local safe_name
    local sample_file
    local header_temp
    local positions_temp
    local sample_count
    local variant_count
    local focal_present
    local focal_count
    local format_ids
    local first_chrom
    local first_pos
    local last_chrom
    local last_pos
    local warning_count

    project_root="/scratch/martinlab/joe/GynandroCalling"
    raw_root="${project_root}/data/raw"
    checkpoint_dir="${project_root}/analysis/checkpoints"
    sample_dir="${checkpoint_dir}/vcf_sample_lists"

    summary="${checkpoint_dir}/region_vcf_audit.tsv"
    ranges="${checkpoint_dir}/region_vcf_contig_ranges.tsv"
    warning_log="${checkpoint_dir}/region_vcf_audit.bcftools_warnings.log"
    focal_expected="${checkpoint_dir}/expected_focal_samples.txt"

    mkdir -p "${checkpoint_dir}" "${sample_dir}"

    if ! command -v bcftools >/dev/null 2>&1
    then
        printf 'ERROR: bcftools is unavailable.\n' >&2
        return 1
    fi

    cat > "${focal_expected}" <<'SAMPLES'
TN2L_S1
TN2R_S2
TN3L_S3
TN3R_S4
TN_Left
TN_Right
SAMPLES

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "file" \
        "bytes" \
        "samples" \
        "variants" \
        "focal_samples" \
        "focal_n" \
        "first_chrom" \
        "first_pos" \
        "last_chrom" \
        "last_pos" \
        "format_fields" \
        > "${summary}"

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "file" \
        "chrom" \
        "variants" \
        "min_pos" \
        "max_pos" \
        > "${ranges}"

    : > "${warning_log}"

    file_count=0
    audit_failed=0

    while IFS= read -r vcf
    do
        file_count=$((file_count + 1))
        relative="${vcf#${project_root}/}"
        safe_name="$(printf '%s' "${relative}" | tr '/ ' '__')"

        sample_file="${sample_dir}/${safe_name}.samples.txt"
        header_temp="$(mktemp "${checkpoint_dir}/vcf_header.XXXXXX")"
        positions_temp="$(mktemp "${checkpoint_dir}/vcf_positions.XXXXXX")"

        if ! bcftools view -h "${vcf}" \
            > "${header_temp}" \
            2>> "${warning_log}"
        then
            printf 'ERROR: unreadable VCF recorded in warning log: %s\n' \
                "${relative}" >&2

            rm -f "${header_temp}" "${positions_temp}"
            audit_failed=1
            continue
        fi

        if ! bcftools query -l "${vcf}" \
            > "${sample_file}" \
            2>> "${warning_log}"
        then
            printf 'ERROR: sample audit failed: %s\n' \
                "${relative}" >&2

            rm -f "${header_temp}" "${positions_temp}"
            audit_failed=1
            continue
        fi

        if ! bcftools query \
            -f '%CHROM\t%POS\n' \
            "${vcf}" \
            > "${positions_temp}" \
            2>> "${warning_log}"
        then
            printf 'ERROR: positional audit failed: %s\n' \
                "${relative}" >&2

            rm -f "${header_temp}" "${positions_temp}"
            audit_failed=1
            continue
        fi

        sample_count="$(wc -l < "${sample_file}")"
        variant_count="$(wc -l < "${positions_temp}")"

        focal_present="$(
            grep -Fxf "${focal_expected}" "${sample_file}" 2>/dev/null |
                paste -sd, -
        )"

        if [[ -n "${focal_present}" ]]
        then
            focal_count="$(
                grep -Fxf "${focal_expected}" "${sample_file}" |
                    wc -l
            )"
        else
            focal_present="none"
            focal_count=0
        fi

        format_ids="$(
            sed -n \
                's/^##FORMAT=<ID=\([^,>]*\).*/\1/p' \
                "${header_temp}" |
                paste -sd, -
        )"

        if [[ -z "${format_ids}" ]]
        then
            format_ids="none_declared"
        fi

        first_chrom="NA"
        first_pos="NA"
        last_chrom="NA"
        last_pos="NA"

        if [[ -s "${positions_temp}" ]]
        then
            IFS=$'\t' read -r first_chrom first_pos \
                < "${positions_temp}"

            IFS=$'\t' read -r last_chrom last_pos \
                < <(tail -n 1 "${positions_temp}")

            awk -v file="${relative}" '
                {
                    count[$1]++
                    if (!($1 in min) || $2 < min[$1]) {
                        min[$1] = $2
                    }
                    if (!($1 in max) || $2 > max[$1]) {
                        max[$1] = $2
                    }
                }
                END {
                    for (chrom in count) {
                        print file "\t" chrom "\t" count[chrom] \
                            "\t" min[chrom] "\t" max[chrom]
                    }
                }
            ' "${positions_temp}" |
                sort -k2,2V \
                >> "${ranges}"
        fi

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${relative}" \
            "$(stat -c '%s' "${vcf}")" \
            "${sample_count}" \
            "${variant_count}" \
            "${focal_present}" \
            "${focal_count}" \
            "${first_chrom}" \
            "${first_pos}" \
            "${last_chrom}" \
            "${last_pos}" \
            "${format_ids}" \
            >> "${summary}"

        rm -f "${header_temp}" "${positions_temp}"

    done < <(
        find \
            "${raw_root}/region_panel" \
            "${raw_root}/ivory_variants" \
            -type f \
            -name '*.vcf.gz' \
            -print |
            sort
    )

    if [[ "${file_count}" -eq 0 ]]
    then
        printf 'ERROR: no extracted VCFs were found.\n' >&2
        return 1
    fi

    if [[ "${audit_failed}" -ne 0 ]]
    then
        printf 'ERROR: one or more VCFs failed validation.\n' >&2
        printf 'Details: %s\n' "${warning_log}" >&2
        return 1
    fi

    printf '\nREGION VCF AUDIT\n'
    column -t -s $'\t' "${summary}"

    printf '\nCONTIG RANGES\n'
    column -t -s $'\t' "${ranges}"

    warning_count="$(wc -l < "${warning_log}")"

    printf '\nAudit completed for %s VCFs.\n' "${file_count}"
    printf 'Full sample lists: %s\n' "${sample_dir}"
    printf 'bcftools warnings: %s lines saved in %s\n' \
        "${warning_count}" \
        "${warning_log}"
}

main "$@"
