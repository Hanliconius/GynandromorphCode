#!/usr/bin/env bash
#SBATCH --job-name=gyn194sibsb
#SBATCH --partition=nano
#SBATCH --array=0-20%10
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=00:30:00
#SBATCH --output=logs/194d_validate_full_sib_call_subsets.%A_%a.out
#SBATCH --error=logs/194d_validate_full_sib_call_subsets.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
bcftools=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/bcftools

source_root="${project}/analysis/genomewide_mixing_input/variant_calls/joint_24/by_chromosome"
set_dir="${project}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata/bam_set"
output_root="${project}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/variant_calls/by_chromosome"

chromosome=$((SLURM_ARRAY_TASK_ID + 1))
started_epoch=$(date +%s)

cd "${project}"

printf 'STAGE_START\tstage=validate_full_sib_call_subset\tchromosome=%s\tjob=%s\ttask=%s\thost=%s\ttime=%s\n' \
    "${chromosome}" \
    "${SLURM_ARRAY_JOB_ID}" \
    "${SLURM_ARRAY_TASK_ID}" \
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
completion="${chrom_dir}/call.complete.tsv"
observed_samples="${chrom_dir}/observed.samples"

input_status=0

for input in \
    "${alleles}" \
    "${expected_samples}" \
    "${joint_bcf}" \
    "${joint_index}"
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

"${bcftools}" view \
    -h \
    "${joint_bcf}" \
    > /dev/null

observed_tmp="${observed_samples}.tmp.${SLURM_JOB_ID}"

"${bcftools}" query \
    -l \
    "${joint_bcf}" \
    > "${observed_tmp}"

sample_count=$(wc -l < "${observed_tmp}")

if [[ "${sample_count}" -ne 14 ]]
then
    printf 'ERROR: expected 14 samples; found %s.\n' \
        "${sample_count}" >&2
    false
fi

if ! cmp \
    -s \
    "${expected_samples}" \
    "${observed_tmp}"
then
    printf 'ERROR: sample order mismatch for chromosome %s.\n' \
        "${chromosome}" >&2

    diff \
        -u \
        "${expected_samples}" \
        "${observed_tmp}" >&2 || true

    false
fi

expected_sites=$(
    gzip -dc "${alleles}" \
        | awk '$0 !~ /^#/ {count++} END {print count + 0}'
)

observed_sites=$(
    "${bcftools}" index \
        -n \
        "${joint_bcf}"
)

if [[ \
    "${observed_sites}" -le 0 || \
    "${observed_sites}" -gt "${expected_sites}" \
]]
then
    printf 'ERROR: impossible site counts: observed=%s expected=%s\n' \
        "${observed_sites}" \
        "${expected_sites}" >&2
    false
fi

extra_sites=$(
    awk '
        NR == FNR {
            if ($0 !~ /^#/) {
                target[$0] = 1
            }
            next
        }

        !($0 in target) {
            extras++
        }

        END {
            print extras + 0
        }
    ' \
        <(gzip -dc "${alleles}") \
        <(
            "${bcftools}" query \
                -f '%CHROM\t%POS\t%REF\t%ALT\n' \
                "${joint_bcf}"
        )
)

if [[ "${extra_sites}" -ne 0 ]]
then
    printf 'ERROR: found %s emitted sites absent from the target allele set.\n' \
        "${extra_sites}" >&2
    false
fi

retained_fraction=$(
    awk \
        -v observed="${observed_sites}" \
        -v expected="${expected_sites}" \
        'BEGIN {printf "%.10f", observed / expected}'
)

mv -f \
    "${observed_tmp}" \
    "${observed_samples}"

completion_tmp="${completion}.tmp.${SLURM_JOB_ID}"

printf 'chromosome\tcontig\tsamples\ttarget_sites\temitted_sites\tretained_fraction\textra_sites\tbcf\tbcf_bytes\tcompleted_at\n' \
    > "${completion_tmp}"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${chromosome}" \
    "${contig}" \
    "${sample_count}" \
    "${expected_sites}" \
    "${observed_sites}" \
    "${retained_fraction}" \
    "${extra_sites}" \
    "${joint_bcf}" \
    "$(wc -c < "${joint_bcf}")" \
    "$(date --iso-8601=seconds)" \
    >> "${completion_tmp}"

mv -f \
    "${completion_tmp}" \
    "${completion}"

ended_epoch=$(date +%s)

printf 'SUBSET_VALIDATED\tchromosome=%s\tcontig=%s\ttarget=%s\temitted=%s\tretained_fraction=%s\textra=%s\n' \
    "${chromosome}" \
    "${contig}" \
    "${expected_sites}" \
    "${observed_sites}" \
    "${retained_fraction}" \
    "${extra_sites}"

printf 'STAGE_END\tstage=validate_full_sib_call_subset\tchromosome=%s\tstatus=0\telapsed_seconds=%s\thost=%s\ttime=%s\n' \
    "${chromosome}" \
    "$((ended_epoch - started_epoch))" \
    "$(hostname)" \
    "$(date --iso-8601=seconds)"
