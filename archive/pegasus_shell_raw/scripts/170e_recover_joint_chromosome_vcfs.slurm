#!/usr/bin/env bash
#SBATCH --job-name=gyn170rc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/170e_recover_joint_chromosome_vcfs.%A_%a.out
#SBATCH --error=logs/170e_recover_joint_chromosome_vcfs.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
bcftools="${tool_bin}/bcftools"
manifest="${project}/analysis/genomewide_mixing_input/metadata/chromosome_calling_manifest.tsv"
by_chromosome="${project}/analysis/genomewide_mixing_input/variant_calls/joint_24/by_chromosome"
original_script="${project}/scripts/170d_call_joint_24sample_chromosome_vcfs.slurm"
task_id="${SLURM_ARRAY_TASK_ID}"
stage_name=recover_joint_chromosome_vcf
start_seconds=${SECONDS}

printf 'STAGE_START\tstage=%s\tjob=%s\ttask=%s\thost=%s\ttime=%s\n' \
    "${stage_name}" "${SLURM_JOB_ID:-NA}" "${task_id}" \
    "$(hostname)" "$(date --iso-8601=seconds)"
trap 'status=$?; printf "STAGE_END\tstage=%s\tstatus=%s\telapsed_seconds=%s\thost=%s\ttime=%s\n" "${stage_name}" "${status}" "$((SECONDS - start_seconds))" "$(hostname)" "$(date --iso-8601=seconds)"' EXIT

row=$(
    awk -F '\t' -v task_id="${task_id}" \
        'NR > 1 && $1 == task_id {print; found++}
         END {if (found != 1) exit 1}' \
        "${manifest}"
)

IFS=$'\t' read -r array_id molecule contig length chromosome_type <<< "${row}"
chrom_dir="${by_chromosome}/chr${molecule}_${contig}"
control_bcf="${chrom_dir}/controls.discovery.biallelic_snps.bcf"
joint_bcf="${chrom_dir}/all_24_samples.joint_at_control_sites.bcf"
old_candidate="${joint_bcf}.new"
recovery_bcf="${chrom_dir}/all_24_samples.joint_at_control_sites.recovery.bcf"
complete="${chrom_dir}/calling.complete.tsv"

validate_and_record() {
    local candidate="$1"
    local candidate_index="${candidate}.csi"
    local sample_count
    local joint_sites
    local control_sites
    local wrong_contigs
    local site_md5
    local sample_md5

    [[ -s "${candidate}" ]] || return 1
    [[ -s "${control_bcf}" && -s "${control_bcf}.csi" ]] || return 1

    rm -f "${candidate_index}"
    "${bcftools}" index -f "${candidate}" || return 1

    sample_count=$("${bcftools}" query -l "${candidate}" | awk 'END {print NR + 0}') || return 1
    joint_sites=$("${bcftools}" index -n "${candidate}") || return 1
    control_sites=$("${bcftools}" index -n "${control_bcf}") || return 1
    wrong_contigs=$(
        "${bcftools}" query -f '%CHROM\n' "${candidate}" \
            | awk -v expected="${contig}" '$1 != expected {bad++} END {print bad + 0}'
    ) || return 1

    if [[ "${sample_count}" -ne 24 || "${joint_sites}" -ne "${control_sites}" || "${wrong_contigs}" -ne 0 ]]
    then
        return 1
    fi

    site_md5=$(
        "${bcftools}" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${candidate}" \
            | md5sum | awk '{print $1}'
    ) || return 1
    sample_md5=$(
        "${bcftools}" query -l "${candidate}" \
            | md5sum | awk '{print $1}'
    ) || return 1

    if [[ "${candidate}" != "${joint_bcf}" ]]
    then
        mv "${candidate}" "${joint_bcf}"
        mv "${candidate_index}" "${joint_bcf}.csi"
    fi

    printf 'array_id\tmolecule\tcontig\tchromosome_type\tsamples\tsites\tsite_md5\tsample_md5\n%s\t%s\t%s\t%s\t24\t%s\t%s\t%s\n' \
        "${array_id}" "${molecule}" "${contig}" "${chromosome_type}" \
        "${joint_sites}" "${site_md5}" "${sample_md5}" \
        > "${complete}.new"
    mv "${complete}.new" "${complete}"
    return 0
}

if validate_and_record "${joint_bcf}"
then
    printf 'RECOVERY_NOT_NEEDED\tchromosome=%s\tfile=%s\n' "${molecule}" "${joint_bcf}"
elif [[ -s "${old_candidate}" ]]
then
    rm -f "${recovery_bcf}" "${recovery_bcf}.csi"
    ln "${old_candidate}" "${recovery_bcf}"

    if validate_and_record "${recovery_bcf}"
    then
        printf 'RECOVERED_EXISTING_BCF\tchromosome=%s\tsource=%s\n' "${molecule}" "${old_candidate}"
    else
        rm -f "${recovery_bcf}" "${recovery_bcf}.csi" "${joint_bcf}" "${joint_bcf}.csi" "${complete}"
        rerun=1
    fi
else
    rerun=1
fi

if [[ "${rerun:-0}" -eq 1 ]]
then
    [[ -s "${original_script}" ]] || {
        printf 'ERROR: original chromosome caller is missing: %s\n' "${original_script}" >&2
        false
    }

    fixed_script="${SLURM_TMPDIR:-/tmp}/170d_fixed_${SLURM_JOB_ID}_${task_id}.sh"
    sed 's/${joint_bcf}\.new/${joint_bcf}.recovery.bcf/g' \
        "${original_script}" > "${fixed_script}"
    bash -n "${fixed_script}"
    printf 'RERUNNING_INCOMPLETE_CHROMOSOME\tchromosome=%s\n' "${molecule}"
    bash "${fixed_script}"

    if ! validate_and_record "${joint_bcf}"
    then
        printf 'ERROR: rerun output failed final validation for chromosome %s.\n' "${molecule}" >&2
        false
    fi
fi

printf 'RECOVERY_COMPLETE\tchromosome=%s\tcontig=%s\tfile=%s\n' \
    "${molecule}" "${contig}" "${joint_bcf}"
