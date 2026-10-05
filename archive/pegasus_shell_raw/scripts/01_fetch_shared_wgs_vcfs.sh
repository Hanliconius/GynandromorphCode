#!/usr/bin/env bash
set -u
set -o pipefail

main() {
    local project_root
    local env_name
    local raw_dir
    local checkpoint_dir
    local optix_id
    local ivory_id
    local optix_vcf
    local ivory_vcf
    local failed

    project_root="/scratch/martinlab/joe/GynandroCalling"
    env_name="gynandrocalling"
    raw_dir="${project_root}/data/raw"
    checkpoint_dir="${project_root}/analysis/checkpoints"

    optix_id="1mFbz4jL51Sv27Wtlt2BksbtgUG9cm7mA"
    ivory_id="1NskQz-a4aUNUm7xJxeHSJWuHFy2P_okn"

    optix_vcf="${raw_dir}/Optix_Hmr18-938195-1239421.vcf"
    ivory_vcf="${raw_dir}/Ivory_Hmr15-1627000-1637000_biallelic.vcf"

    mkdir -p "${raw_dir}" "${checkpoint_dir}"

    if ! micromamba run -n "${env_name}" gdown --version >/dev/null 2>&1
    then
        printf 'ERROR: gdown is unavailable in environment %s.\n' \
            "${env_name}" >&2
        return 1
    fi

    download_and_validate_vcf() {
        local file_id="$1"
        local output="$2"
        local temporary="${output}.part"

        if [[ -s "${output}" ]]
        then
            if micromamba run -n "${env_name}" \
                bcftools view -h "${output}" >/dev/null 2>&1
            then
                printf 'Reusing valid existing input: %s\n' "${output}"
                return 0
            else
                printf 'ERROR: existing file is not a readable VCF: %s\n' \
                    "${output}" >&2
                return 1
            fi
        fi

        rm -f "${temporary}"

        printf 'Downloading Drive file %s\n' "${file_id}"

        if micromamba run -n "${env_name}" \
            gdown \
            "${file_id}" \
            -O "${temporary}"
        then
            if [[ -s "${temporary}" ]] &&
                micromamba run -n "${env_name}" \
                    bcftools view -h "${temporary}" >/dev/null 2>&1
            then
                mv "${temporary}" "${output}"
                printf 'Downloaded and validated: %s\n' "${output}"
                return 0
            else
                printf 'ERROR: downloaded content is not a readable VCF: %s\n' \
                    "${temporary}" >&2
                rm -f "${temporary}"
                return 1
            fi
        else
            printf 'ERROR: Google Drive download failed for file %s\n' \
                "${file_id}" >&2
            rm -f "${temporary}"
            return 1
        fi
    }

    failed=0

    download_and_validate_vcf "${optix_id}" "${optix_vcf}" || failed=1
    download_and_validate_vcf "${ivory_id}" "${ivory_vcf}" || failed=1

    if [[ "${failed}" -ne 0 ]]
    then
        printf '\nOne or more automatic downloads failed.\n' >&2
        printf 'The corresponding files may require browser-authenticated download.\n' >&2
        printf '\nOptix:\n' >&2
        printf '%s\n' \
            'https://drive.google.com/file/d/1mFbz4jL51Sv27Wtlt2BksbtgUG9cm7mA/view' >&2
        printf '\nIvory:\n' >&2
        printf '%s\n' \
            'https://drive.google.com/file/d/1NskQz-a4aUNUm7xJxeHSJWuHFy2P_okn/view' >&2
        return 1
    fi

    printf 'filebytessamplesvariants\n' \
        > "${checkpoint_dir}/shared_wgs_vcf_audit.tsv"

    for vcf in "${optix_vcf}" "${ivory_vcf}"
    do
        printf '%s\t%s\t%s\t%s\n' \
            "$(basename "${vcf}")" \
            "$(stat -c '%s' "${vcf}")" \
            "$(micromamba run -n "${env_name}" bcftools query -l "${vcf}" | wc -l)" \
            "$(micromamba run -n "${env_name}" bcftools view -H "${vcf}" | wc -l)" \
            >> "${checkpoint_dir}/shared_wgs_vcf_audit.tsv"
    done

    sha256sum "${optix_vcf}" "${ivory_vcf}" \
        > "${checkpoint_dir}/shared_wgs_vcf_sha256.txt"

    printf '\nSHARED WGS INPUT AUDIT\n'
    column -t -s $'\t' \
        "${checkpoint_dir}/shared_wgs_vcf_audit.tsv"

    printf '\nChecksums written to:\n%s\n' \
        "${checkpoint_dir}/shared_wgs_vcf_sha256.txt"
}

main "$@"
