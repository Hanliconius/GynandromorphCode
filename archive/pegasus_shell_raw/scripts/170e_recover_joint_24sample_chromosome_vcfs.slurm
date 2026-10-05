#!/usr/bin/env bash
#SBATCH --job-name=gyn170rx
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/170e_recover_joint_24sample_chromosome_vcfs.%A_%a.out
#SBATCH --error=logs/170e_recover_joint_24sample_chromosome_vcfs.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
environment_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin

bcftools="${environment_bin}/bcftools"

manifest="${project}/analysis/genomewide_mixing_input/metadata/chromosome_calling_manifest.tsv"
input_root="${project}/analysis/genomewide_mixing_input"
output_root="${input_root}/variant_calls/joint_24/by_chromosome"

array_task_id="${SLURM_ARRAY_TASK_ID}"

row=$(
    awk \
        -F '\t' \
        -v task_id="${array_task_id}" \
        'NR > 1 && $1 == task_id {print; found++}
         END {
             if (found != 1) {
                 exit 1
             }
         }' \
        "${manifest}"
)

if [[ -z "${row}" ]]
then
    printf 'ERROR: no unique manifest row for array task %s\n' \
        "${array_task_id}" >&2
    false
fi

IFS=$'\t' read -r \
    array_id \
    molecule \
    contig \
    length \
    chromosome_type \
    <<< "${row}"

if [[ "${array_id}" != "${array_task_id}" ]]
then
    printf 'ERROR: manifest array ID mismatch: expected %s, found %s\n' \
        "${array_task_id}" \
        "${array_id}" >&2
    false
fi

chrom_dir="${output_root}/chr${molecule}_${contig}"

control_bcf="${chrom_dir}/controls.discovery.biallelic_snps.bcf"
joint_bcf="${chrom_dir}/all_24_samples.joint_at_control_sites.bcf"
old_temporary_bcf="${joint_bcf}.new"
recovery_bcf="${chrom_dir}/all_24_samples.joint_at_control_sites.recovery.bcf"
completion_file="${chrom_dir}/calling.complete.tsv"

if [[
    -s "${joint_bcf}" &&
    -s "${joint_bcf}.csi" &&
    -s "${completion_file}"
]]
then
    printf 'ALREADY_COMPLETE\tchromosome=%s\tfile=%s\n' \
        "${molecule}" \
        "${joint_bcf}"
else
    if [[ ! -s "${control_bcf}" ]]
    then
        printf 'ERROR: missing control-discovery BCF: %s\n' \
            "${control_bcf}" >&2
        false
    fi

    if [[ ! -s "${control_bcf}.csi" ]]
    then
        printf 'ERROR: missing control-discovery index: %s\n' \
            "${control_bcf}.csi" >&2
        false
    fi

    mapfile -t sample_tables < <(
        find "${input_root}" \
            -type f \
            -name 'all_24_samples.tsv' \
            -print
    )

    if [[ "${#sample_tables[@]}" -ne 1 ]]
    then
        printf 'ERROR: expected exactly one all_24_samples.tsv; found %s\n' \
            "${#sample_tables[@]}" >&2

        printf '%s\n' "${sample_tables[@]}" >&2
        false
    fi

    sample_table="${sample_tables[0]}"

    promote_recovery=0

    if [[ -s "${joint_bcf}" ]]
    then
        candidate_bcf="${joint_bcf}"
    elif [[ -s "${old_temporary_bcf}" ]]
    then
        rm -f \
            "${recovery_bcf}" \
            "${recovery_bcf}.csi"

        ln \
            "${old_temporary_bcf}" \
            "${recovery_bcf}"

        candidate_bcf="${recovery_bcf}"
        promote_recovery=1
    else
        printf 'ERROR: neither final nor temporary joint BCF exists for chromosome %s\n' \
            "${molecule}" >&2
        false
    fi

    expected_bgzf_eof=1f8b08040000000000ff0600424302001b00030000000000000000

    observed_bgzf_eof=$(
        tail -c 28 "${candidate_bcf}" \
            | od -An -tx1 \
            | tr -d '[:space:]'
    )

    if [[ "${observed_bgzf_eof}" != "${expected_bgzf_eof}" ]]
    then
        printf 'ERROR: BCF lacks a complete BGZF EOF marker: %s\n' \
            "${candidate_bcf}" >&2
        false
    fi

    "${bcftools}" view \
        -h \
        "${candidate_bcf}" \
        >/dev/null

    "${bcftools}" index \
        -f \
        "${candidate_bcf}"

    sample_count=$(
        "${bcftools}" query \
            -l \
            "${candidate_bcf}" \
            | awk 'END {print NR + 0}'
    )

    control_sites=$(
        "${bcftools}" index \
            -n \
            "${control_bcf}"
    )

    joint_sites=$(
        "${bcftools}" index \
            -n \
            "${candidate_bcf}"
    )

    expected_samples="${chrom_dir}/expected.samples"
    observed_samples="${chrom_dir}/observed.samples"

    awk \
        -F '\t' \
        'NR > 1 {print $1}' \
        "${sample_table}" \
        > "${expected_samples}.new"

    "${bcftools}" query \
        -l \
        "${candidate_bcf}" \
        > "${observed_samples}.new"

    if ! cmp \
        -s \
        "${expected_samples}.new" \
        "${observed_samples}.new"
    then
        printf 'ERROR: joint BCF sample order is incorrect for chromosome %s.\n' \
            "${molecule}" >&2

        diff \
            -u \
            "${expected_samples}.new" \
            "${observed_samples}.new" >&2 || true

        false
    fi

    if [[
        "${sample_count}" -ne 24 ||
        "${joint_sites}" -ne "${control_sites}"
    ]]
    then
        printf 'ERROR: chromosome %s expected 24 samples and %s sites; found %s samples and %s sites.\n' \
            "${molecule}" \
            "${control_sites}" \
            "${sample_count}" \
            "${joint_sites}" >&2
        false
    fi

    site_md5=$(
        "${bcftools}" query \
            -f '%CHROM\t%POS\t%REF\t%ALT\n' \
            "${candidate_bcf}" \
            | md5sum \
            | awk '{print $1}'
    )

    sample_md5=$(
        "${bcftools}" query \
            -l \
            "${candidate_bcf}" \
            | md5sum \
            | awk '{print $1}'
    )

    if [[ "${promote_recovery}" -eq 1 ]]
    then
        mv \
            "${recovery_bcf}" \
            "${joint_bcf}"

        mv \
            "${recovery_bcf}.csi" \
            "${joint_bcf}.csi"

        rm -f \
            "${old_temporary_bcf}" \
            "${old_temporary_bcf}.csi"
    else
        "${bcftools}" index \
            -f \
            "${joint_bcf}"
    fi

    mv \
        "${expected_samples}.new" \
        "${expected_samples}"

    mv \
        "${observed_samples}.new" \
        "${observed_samples}"

    printf 'array_id\tmolecule\tcontig\tchromosome_type\tsamples\tsites\tsite_md5\tsample_md5\n%s\t%s\t%s\t%s\t24\t%s\t%s\t%s\n' \
        "${array_id}" \
        "${molecule}" \
        "${contig}" \
        "${chromosome_type}" \
        "${joint_sites}" \
        "${site_md5}" \
        "${sample_md5}" \
        > "${completion_file}.new"

    mv \
        "${completion_file}.new" \
        "${completion_file}"

    printf 'RECOVERED_JOINT_24_SAMPLE_CALL\tchromosome=%s\ttype=%s\tsites=%s\n' \
        "${molecule}" \
        "${chromosome_type}" \
        "${joint_sites}"
fi
