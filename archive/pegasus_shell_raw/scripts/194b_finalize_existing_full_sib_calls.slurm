#!/usr/bin/env bash
#SBATCH --job-name=gyn194sibfx
#SBATCH --partition=nano
#SBATCH --array=0-20%21
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:30:00
#SBATCH --output=logs/194b_finalize_existing_full_sib_calls.%A_%a.out
#SBATCH --error=logs/194b_finalize_existing_full_sib_calls.%A_%a.err

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

printf 'STAGE_START\tstage=finalize_existing_full_sib_call\tchromosome=%s\tjob=%s\ttask=%s\thost=%s\ttime=%s\n' \
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
    printf 'ERROR: expected one source chromosome directory; found %s for chromosome %s\n' \
        "${#source_directories[@]}" \
        "${chromosome}" >&2
    false
fi

source_directory=${source_directories[0]}

if [[ ! -d "${source_directory}" ]]
then
    printf 'ERROR: source path is not a directory: %s\n' \
        "${source_directory}" >&2
    false
fi

contig=$(
    basename "${source_directory}"
)
contig=${contig#chr${chromosome}_}

alleles="${source_directory}/controls.discovery.alleles.tsv.gz"
expected_samples="${set_dir}/full_sib_14_samples.samples"
chrom_dir="${output_root}/chr${chromosome}_${contig}"

joint_bcf="${chrom_dir}/full_sib_14.joint_at_existing_control_sites.bcf"
new_bcf="${joint_bcf}.new"
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

expected_sample_count=$(
    wc -l < "${expected_samples}"
)

if [[ "${expected_sample_count}" -ne 14 ]]
then
    printf 'ERROR: expected 14 sample names; found %s\n' \
        "${expected_sample_count}" >&2
    false
fi

expected_sites=$(
    gzip -dc "${alleles}" \
        | awk '$0 !~ /^#/ {count++} END {print count + 0}'
)

candidate=""

if [[ -s "${new_bcf}" ]]
then
    candidate="${new_bcf}"
    candidate_type=new
elif [[ -s "${joint_bcf}" ]]
then
    candidate="${joint_bcf}"
    candidate_type=final
else
    printf 'ERROR: neither temporary nor final BCF exists for chromosome %s\n' \
        "${chromosome}" >&2
    false
fi

printf 'CANDIDATE\tchromosome=%s\ttype=%s\tbytes=%s\tpath=%s\n' \
    "${chromosome}" \
    "${candidate_type}" \
    "$(wc -c < "${candidate}")" \
    "${candidate}"

"${bcftools}" view \
    -h \
    "${candidate}" \
    > /dev/null

candidate_index="${candidate}.csi"
rm -f "${candidate_index}"

"${bcftools}" index \
    -f \
    -o "${candidate_index}" \
    "${candidate}"

if [[ ! -s "${candidate_index}" ]]
then
    printf 'ERROR: explicit CSI creation failed: %s\n' \
        "${candidate_index}" >&2
    false
fi

observed_tmp="${observed_samples}.tmp.${SLURM_JOB_ID}"

"${bcftools}" query \
    -l \
    "${candidate}" \
    > "${observed_tmp}"

sample_count=$(
    wc -l < "${observed_tmp}"
)

if ! cmp \
    -s \
    "${expected_samples}" \
    "${observed_tmp}"
then
    printf 'ERROR: sample order differs from expected order.\n' >&2

    diff \
        -u \
        "${expected_samples}" \
        "${observed_tmp}" >&2 || true

    false
fi

observed_sites=$(
    "${bcftools}" index \
        -n \
        "${candidate}"
)

if [[ \
    "${sample_count}" -ne 14 || \
    "${observed_sites}" -ne "${expected_sites}" \
]]
then
    printf 'ERROR: chromosome %s expected 14 samples and %s sites; found %s samples and %s sites.\n' \
        "${chromosome}" \
        "${expected_sites}" \
        "${sample_count}" \
        "${observed_sites}" >&2
    false
fi

if [[ "${candidate_type}" == "new" ]]
then
    mv -f \
        "${new_bcf}" \
        "${joint_bcf}"

    mv -f \
        "${new_bcf}.csi" \
        "${joint_bcf}.csi"
fi

mv -f \
    "${observed_tmp}" \
    "${observed_samples}"

"${bcftools}" view \
    -h \
    "${joint_bcf}" \
    > /dev/null

final_sites=$(
    "${bcftools}" index \
        -n \
        "${joint_bcf}"
)

if [[ "${final_sites}" -ne "${expected_sites}" ]]
then
    printf 'ERROR: final BCF site-count validation failed.\n' >&2
    false
fi

completion_tmp="${completion}.tmp.${SLURM_JOB_ID}"

printf 'chromosome\tcontig\tsamples\tsites\tbcf\tbcf_bytes\tcompleted_at\n' \
    > "${completion_tmp}"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${chromosome}" \
    "${contig}" \
    "${sample_count}" \
    "${final_sites}" \
    "${joint_bcf}" \
    "$(wc -c < "${joint_bcf}")" \
    "$(date --iso-8601=seconds)" \
    >> "${completion_tmp}"

mv -f \
    "${completion_tmp}" \
    "${completion}"

ended_epoch=$(date +%s)

printf 'FINALIZED\tchromosome=%s\tcontig=%s\tsamples=%s\tsites=%s\tbcf_bytes=%s\n' \
    "${chromosome}" \
    "${contig}" \
    "${sample_count}" \
    "${final_sites}" \
    "$(wc -c < "${joint_bcf}")"

printf 'STAGE_END\tstage=finalize_existing_full_sib_call\tchromosome=%s\tstatus=0\telapsed_seconds=%s\thost=%s\ttime=%s\n' \
    "${chromosome}" \
    "$((ended_epoch - started_epoch))" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"
