#!/usr/bin/env bash

project_dir="/scratch/martinlab/joe/GynandroCalling"

source_vcf="${project_dir}/data/raw/Optix_Hmr18-938195-1239421.vcf"
output_dir="${project_dir}/data/derived/reconstructed_plot_inputs"
checkpoint_dir="${project_dir}/analysis/checkpoints"

summary="${checkpoint_dir}/reconstructed_plot_specific_optix_vcfs.tsv"
checksums="${checkpoint_dir}/reconstructed_plot_specific_optix_vcfs.sha256"

mkdir -p \
    "${output_dir}" \
    "${checkpoint_dir}"

make_slice() {
    output_name="$1"
    region="$2"
    expected_records="$3"

    output_vcf="${output_dir}/${output_name}"
    temporary_vcf="${output_vcf}.part"

    action="not_created"
    status="error"
    observed_records=""
    sample_count=""
    leg_samples=""
    pooled_samples=""

    if [[ -s "${output_vcf}" ]]
    then
        observed_records="$(
            bcftools view -H "${output_vcf}" 2>/dev/null \
                | wc -l
        )"

        if [[ "${observed_records}" -eq "${expected_records}" ]]
        then
            action="reused"
            status="ok"
        else
            action="existing_file_failed_record_validation"
        fi
    else
        if bcftools view \
            -t "${region}" \
            -Ov \
            -o "${temporary_vcf}" \
            "${source_vcf}"
        then
            observed_records="$(
                bcftools view -H "${temporary_vcf}" 2>/dev/null \
                    | wc -l
            )"

            if [[ "${observed_records}" -eq "${expected_records}" ]]
            then
                mv "${temporary_vcf}" "${output_vcf}"
                action="created"
                status="ok"
            else
                action="new_file_failed_record_validation"
                rm -f "${temporary_vcf}"
            fi
        else
            action="bcftools_failed"
            rm -f "${temporary_vcf}"
        fi
    fi

    if [[ -s "${output_vcf}" ]]
    then
        sample_count="$(
            bcftools query -l "${output_vcf}" 2>/dev/null \
                | wc -l
        )"

        leg_samples="$(
            bcftools query -l "${output_vcf}" 2>/dev/null \
                | awk '
                    $0 == "TN2L_S1" ||
                    $0 == "TN2R_S2" ||
                    $0 == "TN3L_S3" ||
                    $0 == "TN3R_S4" {
                        found++
                    }
                    END {
                        print found + 0
                    }
                '
        )"

        pooled_samples="$(
            bcftools query -l "${output_vcf}" 2>/dev/null \
                | awk '
                    $0 == "TN_Left" ||
                    $0 == "TN_Right" {
                        found++
                    }
                    END {
                        print found + 0
                    }
                '
        )"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${output_name}" \
        "${region}" \
        "${expected_records}" \
        "${observed_records}" \
        "${sample_count}" \
        "${leg_samples}" \
        "${pooled_samples}" \
        "${action}" \
        "${status}"
}

if [[ -s "${source_vcf}" ]]
then
    {
        printf 'file\tregion\texpected_records\tobserved_records\tsamples\tleg_samples\tpooled_samples\taction\tstatus\n'

        make_slice \
            "optix_band-1112651-1114314.vcf" \
            "Hmr18:1112651-1114314" \
            60

        make_slice \
            "dennis1-1215210-1216124.vcf" \
            "Hmr18:1215210-1216124" \
            35

        make_slice \
            "dennis2-1218557-1221494.vcf" \
            "Hmr18:1218557-1221494" \
            252

        make_slice \
            "ray1-1156746-1160220.vcf" \
            "Hmr18:1156746-1160220" \
            172

        make_slice \
            "ray2-1174547-1187833.vcf" \
            "Hmr18:1174547-1187833" \
            716
    } > "${summary}"

    find "${output_dir}" \
        -maxdepth 1 \
        -type f \
        -name '*.vcf' \
        -print0 \
        | sort -z \
        | xargs -0 sha256sum \
        > "${checksums}"

    printf 'RECONSTRUCTED PLOT-SPECIFIC OPTIX VCFS\n\n'
    column -t -s $'\t' "${summary}"

    printf '\nChecksums: %s\n' "${checksums}"
    printf 'Output directory: %s\n' "${output_dir}"
else
    printf 'ERROR: source VCF is missing:\n%s\n' "${source_vcf}"
fi
