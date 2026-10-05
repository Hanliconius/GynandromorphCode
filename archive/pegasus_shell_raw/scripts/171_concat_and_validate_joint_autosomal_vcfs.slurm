#!/usr/bin/env bash
#SBATCH --job-name=gyn171ct
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/171_concat_and_validate_joint_autosomal_vcfs.%j.out
#SBATCH --error=logs/171_concat_and_validate_joint_autosomal_vcfs.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin

bcftools="${tool_bin}/bcftools"

manifest="${project}/analysis/genomewide_mixing_input/metadata/chromosome_calling_manifest.tsv"
by_chr="${project}/analysis/genomewide_mixing_input/variant_calls/by_chromosome"
final_dir="${project}/analysis/genomewide_mixing_input/variant_calls/genomewide_chromosomes"
audit="${final_dir}/joint_chromosome_calling.audit.tsv"

mkdir -p \
    "${project}/logs" \
    "${final_dir}"

if [[ ! -s "${manifest}" || ! -x "${bcftools}" ]]
then
    printf 'ERROR: chromosome manifest or bcftools is unavailable.\n' >&2
    false
fi

printf 'dataset\tsamples\tsites\tsite_keys_md5\toutput\n' \
    > "${audit}.new"

for set_name in \
    merged_T2_T3 \
    T2_only \
    T3_only
do
    list_file="${final_dir}/${set_name}.chromosome_bcf.list"

    : > "${list_file}.new"

    while IFS=$'\t' read -r \
        array_id \
        molecule \
        contig \
        length
    do
        [[ "${array_id}" == "array_id" ]] && continue

        chrom_bcf="${by_chr}/chr${molecule}_${contig}/${set_name}.joint_at_control_sites.bcf"
        complete="${by_chr}/chr${molecule}_${contig}/calling.complete.tsv"

        if [[
            ! -s "${chrom_bcf}" ||
            ! -s "${chrom_bcf}.csi" ||
            ! -s "${complete}"
        ]]
        then
            printf 'ERROR: incomplete chromosome output: %s\n' \
                "${chrom_bcf}" >&2
            false
        fi

        printf '%s\n' \
            "${chrom_bcf}" \
            >> "${list_file}.new"
    done < "${manifest}"

    mv \
        "${list_file}.new" \
        "${list_file}"

    final_bcf="${final_dir}/${set_name}.joint_at_control_sites.chromosomes.bcf"
    final_vcf="${final_dir}/${set_name}.joint_at_control_sites.chromosomes.vcf.gz"

    "${bcftools}" concat \
        -f "${list_file}" \
        -Ob \
        --threads 1 \
        -o "${final_bcf}.new"


    "${bcftools}" index \
        -f \
        --threads 1 \
        "${final_bcf}.new"

    "${bcftools}" view \
        -Oz \
        --threads 1 \
        -o "${final_vcf}.new" \
        "${final_bcf}.new"

    "${bcftools}" index \
        -f \
        -t \
        --threads 1 \
        "${final_vcf}.new"

    mv \
        "${final_bcf}.new" \
        "${final_bcf}"

    mv \
        "${final_bcf}.new.csi" \
        "${final_bcf}.csi"

    mv \
        "${final_vcf}.new" \
        "${final_vcf}"

    mv \
        "${final_vcf}.new.tbi" \
        "${final_vcf}.tbi"

    samples=$(
        "${bcftools}" query \
            -l \
            "${final_bcf}" \
            | wc -l
    )

    sites=$(
        "${bcftools}" index \
            -n \
            "${final_bcf}"
    )

    site_md5=$(
        "${bcftools}" query \
            -f '%CHROM\t%POS\t%REF\t%ALT\n' \
            "${final_bcf}" \
            | md5sum \
            | awk '{print $1}'
    )

    if [[
        "${samples}" -ne 20 ||
        "${sites}" -le 0
    ]]
    then
        printf 'ERROR: invalid final %s call: samples=%s sites=%s\n' \
            "${set_name}" \
            "${samples}" \
            "${sites}" >&2
        false
    fi

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${set_name}" \
        "${samples}" \
        "${sites}" \
        "${site_md5}" \
        "${final_vcf}" \
        >> "${audit}.new"
done

unique_site_sets=$(
    awk \
        -F '\t' \
        '
            NR > 1 {
                seen[$4] = 1
            }

            END {
                print length(seen)
            }
        ' \
        "${audit}.new"
)

if [[ "${unique_site_sets}" -ne 1 ]]
then
    printf 'ERROR: final datasets do not have identical frozen site keys.\n' >&2
    false
fi

mv \
    "${audit}.new" \
    "${audit}"

printf 'JOINT NUCLEAR-CHROMOSOME VCF PIPELINE COMPLETE\n'
cat "${audit}"

printf '\nSafeguards:\n'
printf 'site_discovery=18_controls_only\n'
printf 'joint_calling=18_controls_plus_two_focal_samples_per_dataset\n'
printf 'per_sample_vcf_merge=False\n'
printf 'chromosome_operation=concat_nonoverlapping_chromosomes\n'
printf 'sex_chromosomes_included=True\n'
printf 'reference=GCF_056956725.1\n'
