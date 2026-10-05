#!/usr/bin/env bash
#SBATCH --job-name=gyn171rc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/171e_concat_recovered_joint_vcfs.%j.out
#SBATCH --error=logs/171e_concat_recovered_joint_vcfs.%j.err

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
temporary_bcf="${final_dir}/all_24_samples.joint_at_control_sites.chromosomes_1_to_21.tmp.bcf"
temporary_vcf="${final_dir}/all_24_samples.joint_at_control_sites.chromosomes_1_to_21.tmp.vcf.gz"
audit="${final_dir}/joint_24_sample_genomewide_call.audit.tsv"
stage_name=concat_recovered_joint_vcfs
start_seconds=${SECONDS}

printf 'STAGE_START\tstage=%s\tjob=%s\thost=%s\ttime=%s\n' \
    "${stage_name}" "${SLURM_JOB_ID:-NA}" "$(hostname)" "$(date --iso-8601=seconds)"
trap 'status=$?; printf "STAGE_END\tstage=%s\tstatus=%s\telapsed_seconds=%s\thost=%s\ttime=%s\n" "${stage_name}" "${status}" "$((SECONDS - start_seconds))" "$(hostname)" "$(date --iso-8601=seconds)"' EXIT

mkdir -p "${final_dir}"
: > "${list_file}.new"
expected_sample_md5=""
expected_total_sites=0
chromosome_files=0

while IFS=$'\t' read -r array_id molecule contig length chromosome_type
do
    [[ "${array_id}" == "array_id" ]] && continue
    chrom_dir="${by_chromosome}/chr${molecule}_${contig}"
    chrom_bcf="${chrom_dir}/all_24_samples.joint_at_control_sites.bcf"
    complete="${chrom_dir}/calling.complete.tsv"

    if [[ ! -s "${chrom_bcf}" || ! -s "${chrom_bcf}.csi" || ! -s "${complete}" ]]
    then
        printf 'ERROR: incomplete chromosome output: %s\n' "${chrom_bcf}" >&2
        false
    fi

    sample_count=$("${bcftools}" query -l "${chrom_bcf}" | awk 'END {print NR + 0}')
    indexed_sites=$("${bcftools}" index -n "${chrom_bcf}")
    recorded_sites=$(awk -F '\t' 'NR == 2 {print $6}' "${complete}")
    sample_md5=$(awk -F '\t' 'NR == 2 {print $8}' "${complete}")

    if [[ "${sample_count}" -ne 24 || "${indexed_sites}" -ne "${recorded_sites}" ]]
    then
        printf 'ERROR: chromosome %s validation mismatch.\n' "${molecule}" >&2
        false
    fi

    if [[ -z "${expected_sample_md5}" ]]
    then
        expected_sample_md5="${sample_md5}"
    elif [[ "${sample_md5}" != "${expected_sample_md5}" ]]
    then
        printf 'ERROR: sample order differs on chromosome %s.\n' "${molecule}" >&2
        false
    fi

    printf '%s\n' "${chrom_bcf}" >> "${list_file}.new"
    expected_total_sites=$((expected_total_sites + indexed_sites))
    chromosome_files=$((chromosome_files + 1))
done < "${manifest}"

if [[ "${chromosome_files}" -ne 21 ]]
then
    printf 'ERROR: expected 21 chromosome BCFs; found %s.\n' "${chromosome_files}" >&2
    false
fi

mv "${list_file}.new" "${list_file}"
rm -f "${temporary_bcf}" "${temporary_bcf}.csi" "${temporary_vcf}" "${temporary_vcf}.tbi"

"${bcftools}" concat -f "${list_file}" -Ob --threads 1 -o "${temporary_bcf}"
"${bcftools}" index -f "${temporary_bcf}"
"${bcftools}" view -Oz --threads 1 -o "${temporary_vcf}" "${temporary_bcf}"
"${bcftools}" index -f -t "${temporary_vcf}"

sample_count=$("${bcftools}" query -l "${temporary_bcf}" | awk 'END {print NR + 0}')
site_count=$("${bcftools}" index -n "${temporary_bcf}")
vcf_sites=$("${bcftools}" index -n "${temporary_vcf}")

if [[ "${sample_count}" -ne 24 || "${site_count}" -ne "${expected_total_sites}" || "${vcf_sites}" -ne "${site_count}" ]]
then
    printf 'ERROR: concatenated output validation failed: samples=%s BCF_sites=%s VCF_sites=%s expected_sites=%s.\n' \
        "${sample_count}" "${site_count}" "${vcf_sites}" "${expected_total_sites}" >&2
    false
fi

mv "${temporary_bcf}" "${final_bcf}"
mv "${temporary_bcf}.csi" "${final_bcf}.csi"
mv "${temporary_vcf}" "${final_vcf}"
mv "${temporary_vcf}.tbi" "${final_vcf}.tbi"

printf 'item\tvalue\n' > "${audit}.new"
printf 'chromosomes\t21\n' >> "${audit}.new"
printf 'autosomes\t20\n' >> "${audit}.new"
printf 'Z_chromosomes\t1\n' >> "${audit}.new"
printf 'samples\t%s\n' "${sample_count}" >> "${audit}.new"
printf 'sites\t%s\n' "${site_count}" >> "${audit}.new"
printf 'sample_md5\t%s\n' "${expected_sample_md5}" >> "${audit}.new"
printf 'final_BCF\t%s\n' "${final_bcf}" >> "${audit}.new"
printf 'final_VCF\t%s\n' "${final_vcf}" >> "${audit}.new"
mv "${audit}.new" "${audit}"

printf 'RECOVERED_GENOMEWIDE_CONCAT_COMPLETE\tsamples=%s\tsites=%s\n' \
    "${sample_count}" "${site_count}"
cat "${audit}"
