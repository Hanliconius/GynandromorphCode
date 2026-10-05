#!/usr/bin/env bash
#SBATCH --job-name=gyn168fq
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --output=logs/168d_merge_focal_fastqs_by_side.%j.out
#SBATCH --error=logs/168d_merge_focal_fastqs_by_side.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
fastq_root="${project}/analysis/genomewide_mixing_input/fastq"

mkdir -p \
    "${fastq_root}/TN_Left_merged" \
    "${fastq_root}/TN_Right_merged"

merge_fastq()
{
    local input_t2=$1
    local input_t3=$2
    local output=$3

    if [[ ! -s "${input_t2}" || ! -s "${input_t3}" ]]
    then
        printf 'ERROR: missing input FASTQ for %s\n' \
            "${output}" >&2
        false
    fi

    gzip -t "${input_t2}"
    gzip -t "${input_t3}"

    if [[ -s "${output}" ]] && gzip -t "${output}"
    then
        printf 'VALID_EXISTING\t%s\n' \
            "${output}"
    else
        cat \
            "${input_t2}" \
            "${input_t3}" \
            > "${output}.new"

        gzip -t "${output}.new"

        mv \
            "${output}.new" \
            "${output}"

        printf 'CREATED\t%s\t%s bytes\n' \
            "${output}" \
            "$(wc -c < "${output}")"
    fi
}

merge_fastq \
    "${fastq_root}/TN_Left_T2/TN_Left_T2_1.fastq.gz" \
    "${fastq_root}/TN_Left_T3/TN_Left_T3_1.fastq.gz" \
    "${fastq_root}/TN_Left_merged/TN_Left_merged_1.fastq.gz"

merge_fastq \
    "${fastq_root}/TN_Left_T2/TN_Left_T2_2.fastq.gz" \
    "${fastq_root}/TN_Left_T3/TN_Left_T3_2.fastq.gz" \
    "${fastq_root}/TN_Left_merged/TN_Left_merged_2.fastq.gz"

merge_fastq \
    "${fastq_root}/TN_Right_T2/TN_Right_T2_1.fastq.gz" \
    "${fastq_root}/TN_Right_T3/TN_Right_T3_1.fastq.gz" \
    "${fastq_root}/TN_Right_merged/TN_Right_merged_1.fastq.gz"

merge_fastq \
    "${fastq_root}/TN_Right_T2/TN_Right_T2_2.fastq.gz" \
    "${fastq_root}/TN_Right_T3/TN_Right_T3_2.fastq.gz" \
    "${fastq_root}/TN_Right_merged/TN_Right_merged_2.fastq.gz"

printf '\nFASTQ-LEVEL SIDE MERGE COMPLETE\n'
printf 'left=T2L+T3L\n'
printf 'right=T2R+T3R\n'
