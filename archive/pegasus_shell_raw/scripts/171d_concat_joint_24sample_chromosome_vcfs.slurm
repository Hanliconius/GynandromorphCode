#!/usr/bin/env bash
#SBATCH --job-name=gyn171j24
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/171d_concat_joint_24sample_chromosome_vcfs.%j.out
#SBATCH --error=logs/171d_concat_joint_24sample_chromosome_vcfs.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin

bcftools="${tool_bin}/bcftools"

manifest="${project}/analysis/genomewide_mixing_input/metadata/chromosome_calling_manifest.tsv"
by_chromosome="${project}/analysis/genomewide_mixing_input/variant_calls/joint_24/by_chromosome"
final_dir="${project}/analysis/genomewide_mixing_input/variant_calls/joint_24/genomewide"
list_file="${final_dir}/chromosome_bcf.list"
final_bcf="${final_dir}/all_24_samples.joint_at_control_sites.chromosomes_1_to_21.bcf"
final_vcf="${final_dir}/all_24_samples.joint_at_control_sites.chromosomes_1_to_21.vcf.gz"
audit="${final_dir}/joint_24_sample_genomewide_call.audit.tsv"

mkdir -p "${final_dir}"

: > "${list_file}.new"

expected_sample_md5=""

while IFS=$'\t' read -r \
    array_id \
    molecule \
    contig \
    length \
    chromosome_type
do
    [[ "${array_id}" == "array_id" ]] && continue

    chrom_dir="${by_chromosome}/chr${molecule}_${contig}"
    chrom_bcf="${chrom_dir}/all_24_samples.joint_at_control_sites.bcf"
    complete="${chrom_dir}/calling.complete.tsv"

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

    sample_md5=$(
        awk \
            -F '\t' \
            'NR == 2 {print $8}' \
            "${complete}"
    )

    if [[ -z "${expected_sample_md5}" ]]
    then
        expected_sample_md5="${sample_md5}"
    elif [[ "${sample_md5}" != "${expected_sample_md5}" ]]
    then
        printf 'ERROR: sample order differs on chromosome %s\n' \
            "${molecule}" >&2
        false
    fi

    printf '%s\n' \
        "${chrom_bcf}" \
        >> "${list_file}.new"
done < "${manifest}"

mv \
    "${list_file}.new" \
    "${list_file}"

chromosome_files=$(
    wc -l < "${list_file}"
)

if [[ "${chromosome_files}" -ne 21 ]]
then
    printf 'ERROR: expected 21 chromosome BCFs; found %s\n' \
        "${chromosome_files}" >&2
    false
fi

"${bcftools}" concat \
    -f "${list_file}" \
    -Ob \
    --threads 1 \
    -o "${final_bcf}.new"

"${bcftools}" index \
    -f \
    "${final_bcf}.new"

"${bcftools}" view \
    -Oz \
    --threads 1 \
    -o "${final_vcf}.new" \
    "${final_bcf}.new"

"${bcftools}" index \
    -f \
    -t \
    "${final_vcf}.new"

sample_count=$(
    "${bcftools}" query \
        -l \
        "${final_bcf}.new" \
        | wc -l
)

site_count=$(
    "${bcftools}" index \
        -n \
        "${final_bcf}.new"
)

if [[
    "${sample_count}" -ne 24 ||
    "${site_count}" -le 0
]]
then
    printf 'ERROR: final VCF validation failed: samples=%s sites=%s\n' \
        "${sample_count}" \
        "${site_count}" >&2
    false
fi

mv "${final_bcf}.new" "${final_bcf}"
mv "${final_bcf}.new.csi" "${final_bcf}.csi"
mv "${final_vcf}.new" "${final_vcf}"
mv "${final_vcf}.new.tbi" "${final_vcf}.tbi"

printf 'item\tvalue\n' \
    > "${audit}.new"

printf 'chromosomes\t21\n' \
    >> "${audit}.new"

printf 'autosomes\t20\n' \
    >> "${audit}.new"

printf 'Z_chromosomes\t1\n' \
    >> "${audit}.new"

printf 'samples\t%s\n' \
    "${sample_count}" \
    >> "${audit}.new"

printf 'controls\t18\n' \
    >> "${audit}.new"

printf 'focal_representations\t6\n' \
    >> "${audit}.new"

printf 'sites\t%s\n' \
    "${site_count}" \
    >> "${audit}.new"

printf 'site_discovery\tcontrols_only\n' \
    >> "${audit}.new"

printf 'sample_vcf_merge\tFalse\n' \
    >> "${audit}.new"

printf 'final_vcf\t%s\n' \
    "${final_vcf}" \
    >> "${audit}.new"

mv "${audit}.new" "${audit}"

printf 'JOINT 24-SAMPLE GENOME-WIDE VCF COMPLETE\n'
cat "${audit}"
