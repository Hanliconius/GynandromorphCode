#!/usr/bin/env bash
#SBATCH --job-name=gyn194sibpr
#SBATCH --partition=nano
#SBATCH --array=0-20%21
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:30:00
#SBATCH --output=logs/194c_promote_and_validate_full_sib_calls.%A_%a.out
#SBATCH --error=logs/194c_promote_and_validate_full_sib_calls.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
bcftools=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/bcftools

source_root="${project}/analysis/genomewide_mixing_input/variant_calls/joint_24/by_chromosome"
set_dir="${project}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata/bam_set"
output_root="${project}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/variant_calls/by_chromosome"

task_id="${SLURM_ARRAY_TASK_ID}"
chromosome=$((task_id + 1))
started_epoch=$(date +%s)

cd "${project}"

printf 'STAGE_START\tstage=promote_full_sib_BCF\tchromosome=%s\tjob=%s\ttask=%s\thost=%s\ttime=%s\n' \
    "${chromosome}" \
    "${SLURM_ARRAY_JOB_ID}" \
    "${task_id}" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"

shopt -s nullglob

source_directories=(
    "${source_root}/chr${chromosome}_"*
)

if [[ "${#source_directories[@]}" -ne 1 ]]
then
    printf 'ERROR: expected one source directory for chromosome %s; found %s\n' \
        "${chromosome}" \
        "${#source_directories[@]}" >&2
    false
fi

source_directory=${source_directories[0]}
contig=$(basename "${source_directory}")
contig=${contig#chr${chromosome}_}

alleles="${source_directory}/controls.discovery.alleles.tsv.gz"
expected_samples="${set_dir}/full_sib_14_samples.samples"
chrom_dir="${output_root}/chr${chromosome}_${contig}"

joint_bcf="${chrom_dir}/full_sib_14.joint_at_existing_control_sites.bcf"
joint_index="${joint_bcf}.csi"
new_bcf="${joint_bcf}.new"
new_index="${new_bcf}.csi"
completion="${chrom_dir}/call.complete.tsv"
observed_samples="${chrom_dir}/observed.samples"

mkdir -p "${chrom_dir}"

input_status=0

for input in \
    "${alleles}" \
    "${expected_samples}"
do
    if [[ ! -s "${input}" ]]
    then
        printf 'ERROR: missing or empty input: %s\n' \
            "${input}" >&2
        input_status=1
    fi
done

if [[ "${input_status}" -ne 0 ]]
then
    false
fi

expected_sample_count=$(wc -l < "${expected_samples}")

if [[ "${expected_sample_count}" -ne 14 ]]
then
    printf 'ERROR: expected 14 samples; found %s\n' \
        "${expected_sample_count}" >&2
    false
fi

expected_sites=$(
    gzip -dc "${alleles}" \
        | awk '$0 !~ /^#/ {count++} END {print count + 0}'
)

if [[ -s "${joint_bcf}" ]]
then
    candidate="${joint_bcf}"
    candidate_type=already_promoted
elif [[ -s "${new_bcf}" && -s "${new_index}" ]]
then
    candidate="${new_bcf}"
    candidate_type=temporary
else
    printf 'ERROR: no complete BCF/index pair found for chromosome %s\n' \
        "${chromosome}" >&2
    false
fi

printf 'CANDIDATE\tchromosome=%s\ttype=%s\tpath=%s\tbytes=%s\n' \
    "${chromosome}" \
    "${candidate_type}" \
    "${candidate}" \
    "$(wc -c < "${candidate}")"

"${bcftools}" view \
    -h \
    "${candidate}" \
    > /dev/null

sample_tmp="${observed_samples}.tmp.${SLURM_JOB_ID}"

"${bcftools}" query \
    -l \
    "${candidate}" \
    > "${sample_tmp}"

sample_count=$(wc -l < "${sample_tmp}")

if ! cmp \
    -s \
    "${expected_samples}" \
    "${sample_tmp}"
then
    printf 'ERROR: sample order mismatch for chromosome %s.\n' \
        "${chromosome}" >&2

    diff \
        -u \
        "${expected_samples}" \
        "${sample_tmp}" >&2 || true

    false
fi

if [[ "${sample_count}" -ne 14 ]]
then
    printf 'ERROR: expected 14 samples; found %s.\n' \
        "${sample_count}" >&2
    false
fi

if [[ "${candidate_type}" == "temporary" ]]
then
    mv -f \
        "${new_bcf}" \
        "${joint_bcf}"

    mv -f \
        "${new_index}" \
        "${joint_index}"

    printf 'PROMOTED\t%s\n' "${joint_bcf}"
fi

if ! "${bcftools}" index -n "${joint_bcf}" >/dev/null 2>&1
then
    printf 'REBUILDING_FINAL_INDEX\t%s\n' "${joint_index}"

    rm -f "${joint_index}"

    "${bcftools}" index \
        -f \
        "${joint_bcf}"
fi

if [[ ! -s "${joint_index}" ]]
then
    printf 'ERROR: final CSI missing: %s\n' \
        "${joint_index}" >&2
    false
fi

observed_sites=$(
    "${bcftools}" index \
        -n \
        "${joint_bcf}"
)

if [[ "${observed_sites}" -ne "${expected_sites}" ]]
then
    printf 'ERROR: chromosome %s expected %s sites; found %s.\n' \
        "${chromosome}" \
        "${expected_sites}" \
        "${observed_sites}" >&2
    false
fi

mv -f \
    "${sample_tmp}" \
    "${observed_samples}"

completion_tmp="${completion}.tmp.${SLURM_JOB_ID}"

printf 'chromosome\tcontig\tsamples\tsites\tbcf\tbcf_bytes\tcompleted_at\n' \
    > "${completion_tmp}"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${chromosome}" \
    "${contig}" \
    "${sample_count}" \
    "${observed_sites}" \
    "${joint_bcf}" \
    "$(wc -c < "${joint_bcf}")" \
    "$(date --iso-8601=seconds)" \
    >> "${completion_tmp}"

mv -f \
    "${completion_tmp}" \
    "${completion}"

ended_epoch=$(date +%s)

printf 'COMPLETE\tchromosome=%s\tcontig=%s\tsamples=%s\tsites=%s\tbcf_bytes=%s\n' \
    "${chromosome}" \
    "${contig}" \
    "${sample_count}" \
    "${observed_sites}" \
    "$(wc -c < "${joint_bcf}")"

printf 'STAGE_END\tstage=promote_full_sib_BCF\tchromosome=%s\tstatus=0\telapsed_seconds=%s\thost=%s\ttime=%s\n' \
    "${chromosome}" \
    "$((ended_epoch - started_epoch))" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"
