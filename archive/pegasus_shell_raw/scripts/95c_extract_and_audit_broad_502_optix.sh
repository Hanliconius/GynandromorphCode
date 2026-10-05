#!/usr/bin/env bash

project_dir="/scratch/martinlab/joe/GynandroCalling"
archive="${project_dir}/Genotype_Plots.tar"

archive_member="Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf"

output_dir="${project_dir}/external/krista_genotype_plots_original/Genotype_Plots"
output_vcf="${output_dir}/Optix_Hmr18-938195-1239421_biallelic.vcf"

checkpoint_dir="${project_dir}/analysis/original_502/checkpoints"
sample_list="${checkpoint_dir}/broad_optix_original502_samples.txt"
summary="${checkpoint_dir}/broad_optix_original502_extraction_audit.tsv"
checksum="${checkpoint_dir}/broad_optix_original502.sha256"

mkdir -p "${output_dir}" \
         "${checkpoint_dir}"

status="ok"
expected_bytes=""

if [[ ! -s "${archive}" ]]
then
    printf 'ERROR: archive is missing or empty:\n%s\n' "${archive}"
    status="failed"
fi

if [[ "${status}" == "ok" ]]
then
    expected_bytes="$(
        tar -tvf "${archive}" "${archive_member}" 2>/dev/null \
            | awk 'NR == 1 {print $3}'
    )"

    if [[ -z "${expected_bytes}" ]]
    then
        printf 'ERROR: archive member was not found:\n%s\n' \
            "${archive_member}"
        status="failed"
    else
        printf 'Archive member: %s\n' "${archive_member}"
        printf 'Expected bytes: %s\n' "${expected_bytes}"
    fi
fi

if [[ "${status}" == "ok" && -s "${output_vcf}" ]]
then
    existing_bytes="$(stat -c '%s' "${output_vcf}")"

    if [[ "${existing_bytes}" == "${expected_bytes}" ]]
    then
        printf 'Reusing existing complete extraction:\n%s\n' \
            "${output_vcf}"
    else
        printf 'ERROR: an existing extraction has the wrong size.\n'
        printf 'File: %s\n' "${output_vcf}"
        printf 'Observed bytes: %s\n' "${existing_bytes}"
        printf 'Expected bytes: %s\n' "${expected_bytes}"
        printf '%s\n' \
            'The existing file was left unchanged for manual inspection.'
        status="failed"
    fi
elif [[ "${status}" == "ok" ]]
then
    temporary_dir="$(
        mktemp -d \
            "${project_dir}/data/incoming/broad_optix_extract.XXXXXX"
    )"

    printf 'Extracting the 502-sample broad Optix VCF.\n'
    printf 'Temporary directory: %s\n' "${temporary_dir}"

    if tar -xf "${archive}" \
        -C "${temporary_dir}" \
        "${archive_member}"
    then
        temporary_vcf="${temporary_dir}/${archive_member}"

        if [[ -s "${temporary_vcf}" ]]
        then
            observed_bytes="$(stat -c '%s' "${temporary_vcf}")"

            if [[ "${observed_bytes}" == "${expected_bytes}" ]]
            then
                mv "${temporary_vcf}" "${output_vcf}"
                printf 'Extraction completed:\n%s\n' "${output_vcf}"
            else
                printf 'ERROR: extracted file size does not match archive metadata.\n'
                printf 'Observed bytes: %s\n' "${observed_bytes}"
                printf 'Expected bytes: %s\n' "${expected_bytes}"
                status="failed"
            fi
        else
            printf '%s\n' 'ERROR: extracted VCF is missing or empty.'
            status="failed"
        fi
    else
        printf '%s\n' 'ERROR: tar extraction failed.'
        status="failed"
    fi

    if [[ -d "${temporary_dir}" ]]
    then
        rm -r "${temporary_dir}"
    fi
fi

if [[ "${status}" == "ok" ]]
then
    printf '\n%s\n' "BROAD ORIGINAL-502 OPTIX VCF AUDIT"

    micromamba run \
        -n gynandrocalling \
        bcftools query \
        -l "${output_vcf}" \
        > "${sample_list}"

    sample_count="$(wc -l < "${sample_list}")"
    file_bytes="$(stat -c '%s' "${output_vcf}")"

    printf 'VCF: %s\n' "${output_vcf}"
    printf 'Bytes: %s\n' "${file_bytes}"
    printf 'Samples: %s\n' "${sample_count}"

    printf '\n%s\n' "FOCAL SAMPLES"

    focal_count=0

    for sample in \
        TN2L_S1 \
        TN3L_S3 \
        TN2R_S2 \
        TN3R_S4 \
        TN_Left \
        TN_Right
    do
        if grep -Fqx "${sample}" "${sample_list}"
        then
            printf '%s\tpresent\n' "${sample}"
            focal_count=$((focal_count + 1))
        else
            printf '%s\tMISSING\n' "${sample}"
        fi
    done

    printf '\n%s\n' "RELEVANT FORMAT DEFINITIONS"

    micromamba run \
        -n gynandrocalling \
        bcftools view \
        -h "${output_vcf}" \
        | grep -E '^##FORMAT=<ID=(GT|AD|DP|GQ|PL|GP),'

    printf '\n%s\n' "VARIANT RANGE"

    range_summary="$(
        micromamba run \
            -n gynandrocalling \
            bcftools query \
            -f '%CHROM\t%POS\n' \
            "${output_vcf}" \
        | awk '
            NR == 1 {
                first_chrom = $1
                first_pos = $2
            }
            {
                last_chrom = $1
                last_pos = $2
                records++
                chroms[$1] = 1
            }
            END {
                chromosome_count = 0
                for (chrom in chroms) {
                    chromosome_count++
                }

                printf "%d\t%d\t%s:%s\t%s:%s",
                    records,
                    chromosome_count,
                    first_chrom,
                    first_pos,
                    last_chrom,
                    last_pos
            }
        '
    )"

    record_count="$(printf '%s\n' "${range_summary}" | cut -f1)"
    chromosome_count="$(printf '%s\n' "${range_summary}" | cut -f2)"
    first_record="$(printf '%s\n' "${range_summary}" | cut -f3)"
    last_record="$(printf '%s\n' "${range_summary}" | cut -f4)"

    printf 'Records: %s\n' "${record_count}"
    printf 'Chromosomes: %s\n' "${chromosome_count}"
    printf 'First record: %s\n' "${first_record}"
    printf 'Last record: %s\n' "${last_record}"

    sha256sum "${output_vcf}" > "${checksum}"

    {
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "file" \
            "bytes" \
            "samples" \
            "focal_samples_present" \
            "records" \
            "chromosomes" \
            "first_record" \
            "last_record"

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${output_vcf}" \
            "${file_bytes}" \
            "${sample_count}" \
            "${focal_count}" \
            "${record_count}" \
            "${chromosome_count}" \
            "${first_record}" \
            "${last_record}"
    } > "${summary}"

    printf '\nAudit summary: %s\n' "${summary}"
    printf 'Sample list: %s\n' "${sample_list}"
    printf 'Checksum: %s\n' "${checksum}"

    if [[ "${sample_count}" -ne 502 ]]
    then
        printf 'WARNING: expected 502 samples but found %s.\n' \
            "${sample_count}"
    fi

    if [[ "${focal_count}" -ne 6 ]]
    then
        printf 'WARNING: expected all 6 focal samples but found %s.\n' \
            "${focal_count}"
    fi
else
    printf '\n%s\n' \
        'Extraction/audit did not complete; no downstream analysis was run.'
fi
