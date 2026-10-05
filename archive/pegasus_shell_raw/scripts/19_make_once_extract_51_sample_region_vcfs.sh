#!/usr/bin/env bash

project_dir="/scratch/martinlab/joe/GynandroCalling"
input_dir="${project_dir}/data/raw"
output_dir="${project_dir}/data/derived/region_genotypes"
checkpoint_dir="${project_dir}/analysis/checkpoints"

mkdir -p "${output_dir}" "${checkpoint_dir}"

inputs=(
    "optix_51sample|Optix_Hmr18-938195-1239421.vcf"
    "ivory_51sample|Ivory_Hmr15-1627000-1637000_biallelic.vcf"
)

summary="${checkpoint_dir}/region_51sample_extraction_summary.tsv"
summary_part="${summary}.part"

printf 'panel\tinput\tsamples\tvariants\texpected_rows\tobserved_rows\taction\tstatus\n' \
    > "${summary_part}"

failures=0

for specification in "${inputs[@]}"
do
    panel="${specification%%|*}"
    filename="${specification#*|}"

    input="${input_dir}/${filename}"
    output="${output_dir}/${panel}.genotypes.tsv.gz"
    output_part="${output}.part"

    samples="NA"
    variants="NA"
    expected_rows="NA"
    observed_rows="NA"
    action="not_attempted"
    status="failed"

    if [[ ! -s "${input}" ]]
    then
        printf 'ERROR: missing input: %s\n' "${input}" >&2
        failures=$((failures + 1))
    else
        samples="$(
            bcftools query -l "${input}" 2>/dev/null \
                | awk 'END {print NR}'
        )"

        variants="$(
            bcftools view -H "${input}" 2>/dev/null \
                | awk 'END {print NR}'
        )"

        if [[ "${samples}" =~ ^[0-9]+$ \
              && "${variants}" =~ ^[0-9]+$ ]]
        then
            expected_rows=$((samples * variants + 1))
            existing_valid=0

            if [[ -s "${output}" ]] \
                && gzip -t "${output}" 2>/dev/null
            then
                observed_rows="$(
                    gzip -cd "${output}" 2>/dev/null \
                        | awk 'END {print NR}'
                )"

                if [[ "${observed_rows}" == "${expected_rows}" ]]
                then
                    existing_valid=1
                fi
            fi

            if [[ "${existing_valid}" -eq 1 ]]
            then
                action="reused"
                status="ok"
            else
                action="created"
                rm -f "${output_part}"

                if {
                    printf 'chrom\tpos\tref\talt\tqual\tfilter\tsample\tGT\tDP\tGQ\tAD\n'

                    bcftools query \
                        -f '[%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%SAMPLE\t%GT\t%DP\t%GQ\t%AD\n]' \
                        "${input}"
                } | bgzip -c > "${output_part}"
                then
                    if [[ -s "${output_part}" ]] \
                        && gzip -t "${output_part}" 2>/dev/null
                    then
                        observed_rows="$(
                            gzip -cd "${output_part}" 2>/dev/null \
                                | awk 'END {print NR}'
                        )"

                        if [[ "${observed_rows}" == "${expected_rows}" ]]
                        then
                            mv "${output_part}" "${output}"
                            status="ok"
                        else
                            printf \
                                'ERROR: %s rows observed; %s expected for %s\n' \
                                "${observed_rows}" \
                                "${expected_rows}" \
                                "${panel}" >&2

                            rm -f "${output_part}"
                            failures=$((failures + 1))
                        fi
                    else
                        printf 'ERROR: invalid output: %s\n' \
                            "${output_part}" >&2

                        rm -f "${output_part}"
                        failures=$((failures + 1))
                    fi
                else
                    printf 'ERROR: extraction failed: %s\n' \
                        "${input}" >&2

                    rm -f "${output_part}"
                    failures=$((failures + 1))
                fi
            fi
        else
            printf 'ERROR: failed to count samples or variants: %s\n' \
                "${input}" >&2
            failures=$((failures + 1))
        fi
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${panel}" \
        "${filename}" \
        "${samples}" \
        "${variants}" \
        "${expected_rows}" \
        "${observed_rows}" \
        "${action}" \
        "${status}" \
        >> "${summary_part}"
done

mv "${summary_part}" "${summary}"

printf '51-SAMPLE REGION VCF EXTRACTION\n'
column -t -s $'\t' "${summary}"

printf '\nOutputs: %s\n' "${output_dir}"
printf 'Summary: %s\n' "${summary}"

if [[ "${failures}" -gt 0 ]]
then
    printf '\nCompleted with %s failure(s).\n' "${failures}" >&2
    false
else
    printf '\nBoth tables passed validation.\n'
    true
fi
