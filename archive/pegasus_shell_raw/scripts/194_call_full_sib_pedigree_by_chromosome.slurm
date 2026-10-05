#!/usr/bin/env bash
#SBATCH --job-name=gyn194sibvc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --output=logs/194_call_full_sib_pedigree.%A_%a.out
#SBATCH --error=logs/194_call_full_sib_pedigree.%A_%a.err

set -euo pipefail

project_root=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin

reference=${project_root}/analysis/genomewide_mixing_input/reference/GCF_056956725.1/GCF_056956725.1_ASM5695672v1_genomic.fna

source_root=${project_root}/analysis/genomewide_mixing_input/variant_calls/joint_24/by_chromosome

set_dir=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata/bam_set

output_root=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/variant_calls/by_chromosome

bcftools=${tool_bin}/bcftools

cd "${project_root}"

mkdir -p \
    "${output_root}" \
    logs

task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
chromosome=$((task_id + 1))

shopt -s nullglob
source_directories=(
    "${source_root}/chr${chromosome}_"*
)
shopt -u nullglob

if [[ "${#source_directories[@]}" -ne 1 ]]
then
    printf 'ERROR: expected one source directory for chromosome %s; found %s.\n' \
        "${chromosome}" \
        "${#source_directories[@]}" >&2
    false
fi

source_directory=${source_directories[0]}

contig=$(
    basename "${source_directory}"
)
contig=${contig#chr${chromosome}_}

alleles=${source_directory}/controls.discovery.alleles.tsv.gz
bamlist=${set_dir}/full_sib_14_samples.bamlist
expected_samples=${set_dir}/full_sib_14_samples.samples

chrom_dir=${output_root}/chr${chromosome}_${contig}
joint_bcf=${chrom_dir}/full_sib_14.joint_at_existing_control_sites.bcf
completion=${chrom_dir}/call.complete.tsv

mkdir -p "${chrom_dir}"

printf 'STAGE_START\tstage=call_full_sib_pedigree\tjob=%s\ttask=%s\tchromosome=%s\tcontig=%s\ttime=%s\n' \
    "${SLURM_ARRAY_JOB_ID}" \
    "${task_id}" \
    "${chromosome}" \
    "${contig}" \
    "$(date --iso-8601=seconds)"

input_status=0

for input in \
    "${reference}" \
    "${reference}.fai" \
    "${alleles}" \
    "${bamlist}" \
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

reuse=0

if [[ \
    -s "${completion}" && \
    -s "${joint_bcf}" && \
    -s "${joint_bcf}.csi" \
]]
then
    if "${bcftools}" index -n "${joint_bcf}" >/dev/null
    then
        reuse=1
    fi
fi

if [[ "${reuse}" -eq 1 ]]
then
    printf 'REUSE\tchromosome=%s\tbcf=%s\n' \
        "${chromosome}" \
        "${joint_bcf}"
else
    rm -f \
        "${joint_bcf}.new" \
        "${joint_bcf}.new.csi"

    "${bcftools}" mpileup \
        -Ou \
        --threads 2 \
        -f "${reference}" \
        -r "${contig}" \
        -T "${alleles}" \
        -q 20 \
        -Q 20 \
        -d 1000 \
        -a FORMAT/AD,FORMAT/DP,FORMAT/SP \
        -b "${bamlist}" \
        | "${bcftools}" call \
            -Ob \
            --threads 2 \
            -m \
            -A \
            -C alleles \
            -T "${alleles}" \
            -G - \
            -o "${joint_bcf}.new"

    "${bcftools}" index \
        -f \
        "${joint_bcf}.new"

    observed_samples=${chrom_dir}/observed.samples.new

    "${bcftools}" query \
        -l \
        "${joint_bcf}.new" \
        > "${observed_samples}"

    if ! cmp \
        -s \
        "${expected_samples}" \
        "${observed_samples}"
    then
        printf 'ERROR: sample order differs from the BAM manifest.\n' >&2

        diff \
            -u \
            "${expected_samples}" \
            "${observed_samples}" >&2 || true

        false
    fi

    expected_sites=$(
        gzip -dc "${alleles}" \
            | awk '
                $0 !~ /^#/ {
                    count++
                }

                END {
                    print count + 0
                }
            '
    )

    observed_sites=$(
        "${bcftools}" index \
            -n \
            "${joint_bcf}.new"
    )

    sample_count=$(
        wc -l \
            < "${observed_samples}"
    )

    if [[ "${sample_count}" -ne 14 || "${observed_sites}" -ne "${expected_sites}" ]]
    then
        printf 'ERROR: expected 14 samples and %s sites; found %s samples and %s sites.\n' \
            "${expected_sites}" \
            "${sample_count}" \
            "${observed_sites}" >&2
        false
    fi

    mv \
        "${joint_bcf}.new" \
        "${joint_bcf}"

    mv \
        "${joint_bcf}.new.csi" \
        "${joint_bcf}.csi"

    mv \
        "${observed_samples}" \
        "${chrom_dir}/observed.samples"

    {
        printf 'chromosome\tcontig\tsamples\tsites\tbcf\tbcf_bytes\tcompleted_at\n'

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${chromosome}" \
            "${contig}" \
            "${sample_count}" \
            "${observed_sites}" \
            "${joint_bcf}" \
            "$(wc -c < "${joint_bcf}")" \
            "$(date --iso-8601=seconds)"
    } > "${completion}.new"

    mv \
        "${completion}.new" \
        "${completion}"

    printf 'COMPLETE\tchromosome=%s\tcontig=%s\tsamples=%s\tsites=%s\n' \
        "${chromosome}" \
        "${contig}" \
        "${sample_count}" \
        "${observed_sites}"
fi

printf 'STAGE_END\tstage=call_full_sib_pedigree\tstatus=0\tchromosome=%s\ttime=%s\n' \
    "${chromosome}" \
    "$(date --iso-8601=seconds)"
