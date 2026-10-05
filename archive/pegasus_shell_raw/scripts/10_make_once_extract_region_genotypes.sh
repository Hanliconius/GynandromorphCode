#!/usr/bin/env bash

project_dir="/scratch/martinlab/joe/GynandroCalling"
input_dir="${project_dir}/data/raw/region_panel"
output_dir="${project_dir}/data/derived/region_genotypes"
checkpoint_dir="${project_dir}/analysis/checkpoints"

mkdir -p "${output_dir}" "${checkpoint_dir}"

vcf_files=(
    "optix_filtered.vcf.gz"
    "ivory_filtered.vcf.gz"
    "WntA_filtered.vcf.gz"
    "neutral1_filtered.vcf.gz"
    "neutral2_filtered.vcf.gz"
    "neutral3_filtered.vcf.gz"
    "neutral4_filtered.vcf.gz"
)

summary="${checkpoint_dir}/region_genotype_extraction_summary.tsv"
summary_part="${summary}.part"

printf 'panel\tvariants\tsamples\texpected_rows\tobserved_rows\tindex_action\ttable_action\tstatus\n' \
    > "${summary_part}"

failures=0

for filename in "${vcf_files[@]}"
do
    vcf="${input_dir}/${filename}"
    panel="${filename%.vcf.gz}"
    panel="${panel%_filtered}"

    output="${output_dir}/${panel}.genotypes.tsv.gz"
    output_part="${output}.part"

    variants="NA"
    samples="NA"
    expected_rows="NA"
    observed_rows="NA"
    index_action="not_checked"
    table_action="not_checked"
    status="failed"

    if [[ ! -s "${vcf}" ]]
    then
        printf 'ERROR: missing or empty input: %s\n' "${vcf}" >&2
        failures=$((failures + 1))
    else
        if [[ -s "${vcf}.csi" || -s "${vcf}.tbi" ]]
        then
            index_action="reused"
        else
            index_action="created"

            if ! bcftools index --threads 2 "${vcf}" >/dev/null 2>&1
            then
                index_action="failed"
                printf 'ERROR: indexing failed: %s\n' "${vcf}" >&2
            fi
        fi

        if [[ -s "${vcf}.csi" || -s "${vcf}.tbi" ]]
        then
            variants="$(bcftools index -n "${vcf}" 2>/dev/null)"
            samples="$(bcftools query -l "${vcf}" 2>/dev/null | awk 'END {print NR}')"

            if [[ "${variants}" =~ ^[0-9]+$ && "${samples}" =~ ^[0-9]+$ ]]
            then
                expected_rows=$((variants * samples + 1))
                existing_valid=0

                if [[ -s "${output}" ]] && gzip -t "${output}" 2>/dev/null
                then
                    observed_rows="$(gzip -cd "${output}" 2>/dev/null | wc -l | awk '{print $1}')"

                    if [[ "${observed_rows}" == "${expected_rows}" ]]
                    then
                        existing_valid=1
                    fi
                fi

                if [[ "${existing_valid}" -eq 1 ]]
                then
                    table_action="reused"
                    status="ok"
                else
                    table_action="created"
                    rm -f "${output_part}"

                    if {
                        printf 'chrom\tpos\tref\talt\tqual\tfilter\tsample\tGT\tDP\tGQ\tGP\tPL\tAD\tADF\tADR\tSP\n'
                        bcftools query \
                            -f '[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%SAMPLE\t%GT\t%DP\t%GQ\t%GP\t%PL\t%AD\t%ADF\t%ADR\t%SP\n]' \
                            "${vcf}"
                    } | bgzip -c > "${output_part}"
                    then
                        if [[ -s "${output_part}" ]] && gzip -t "${output_part}" 2>/dev/null
                        then
                            observed_rows="$(gzip -cd "${output_part}" 2>/dev/null | wc -l | awk '{print $1}')"

                            if [[ "${observed_rows}" == "${expected_rows}" ]]
                            then
                                mv "${output_part}" "${output}"
                                status="ok"
                            else
                                printf 'ERROR: %s has %s rows; expected %s\n' \
                                    "${output_part}" "${observed_rows}" "${expected_rows}" >&2
                                rm -f "${output_part}"
                                failures=$((failures + 1))
                            fi
                        else
                            printf 'ERROR: invalid compressed output: %s\n' "${output_part}" >&2
                            rm -f "${output_part}"
                            failures=$((failures + 1))
                        fi
                    else
                        printf 'ERROR: genotype extraction failed: %s\n' "${vcf}" >&2
                        rm -f "${output_part}"
                        failures=$((failures + 1))
                    fi
                fi
            else
                printf 'ERROR: could not determine variant/sample counts: %s\n' "${vcf}" >&2
                failures=$((failures + 1))
            fi
        else
            failures=$((failures + 1))
        fi
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${panel}" \
        "${variants}" \
        "${samples}" \
        "${expected_rows}" \
        "${observed_rows}" \
        "${index_action}" \
        "${table_action}" \
        "${status}" \
        >> "${summary_part}"
done

mv "${summary_part}" "${summary}"

printf '\nREGION GENOTYPE EXTRACTION\n'
column -t -s $'\t' "${summary}"

printf '\nCompressed genotype tables:\n%s\n' "${output_dir}"
printf 'Summary: %s\n' "${summary}"

if [[ "${failures}" -gt 0 ]]
then
    printf '\nCompleted with %s failure(s). See errors above.\n' "${failures}" >&2
    false
else
    printf '\nAll region genotype tables passed validation.\n'
    true
fi
