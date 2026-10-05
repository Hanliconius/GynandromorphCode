#!/usr/bin/env bash
#SBATCH --job-name=gyn168al
#SBATCH --partition=cpu
#SBATCH --array=0-1
#SBATCH --cpus-per-task=3
#SBATCH --mem=20G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/168e_align_fastq_merged_sides.%A_%a.out
#SBATCH --error=logs/168e_align_fastq_merged_sides.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin

bwa="${tool_bin}/bwa"
samtools="${tool_bin}/samtools"

reference_path_file="${project}/analysis/genomewide_mixing_input/metadata/reference_fasta.path"
fastq_root="${project}/analysis/genomewide_mixing_input/fastq"
output_dir="${project}/analysis/genomewide_mixing_input/alignment/libraries"
tmp_root="${project}/analysis/genomewide_mixing_input/tmp/alignment_merged_sides"

mkdir -p \
    "${output_dir}" \
    "${tmp_root}"

case "${SLURM_ARRAY_TASK_ID}" in
    0)
        sample=TN_Left_merged
        ;;
    1)
        sample=TN_Right_merged
        ;;
    *)
        printf 'ERROR: unexpected array task: %s\n' \
            "${SLURM_ARRAY_TASK_ID}" >&2
        false
        ;;
esac

reference=$(head -n 1 "${reference_path_file}")

if [[
    "${reference}" != *'/GCF_056956725.1/'* ||
    ! -s "${reference}" ||
    ! -s "${reference}.fai"
]]
then
    printf 'ERROR: invalid active reference: %s\n' \
        "${reference}" >&2
    false
fi

fastq1="${fastq_root}/${sample}/${sample}_1.fastq.gz"
fastq2="${fastq_root}/${sample}/${sample}_2.fastq.gz"

for fastq in "${fastq1}" "${fastq2}"
do
    if [[ ! -s "${fastq}" ]]
    then
        printf 'ERROR: missing merged FASTQ: %s\n' \
            "${fastq}" >&2
        false
    fi

    gzip -t "${fastq}"
done

final_bam="${output_dir}/${sample}.dedup.bam"
final_bai="${final_bam}.bai"
flagstat="${output_dir}/${sample}.dedup.flagstat.txt"
stats="${output_dir}/${sample}.dedup.stats.txt"
markdup_stats="${output_dir}/${sample}.markdup.stats.txt"
done_file="${output_dir}/${sample}.complete.tsv"

if [[
    -s "${final_bam}" &&
    -s "${final_bai}" &&
    -s "${done_file}"
]] && "${samtools}" quickcheck -q "${final_bam}"
then
    printf 'VALID_EXISTING\t%s\n' \
        "${final_bam}"
else
    work_dir="${tmp_root}/${sample}.${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}"

    mkdir -p "${work_dir}"

    name_bam="${work_dir}/${sample}.name.bam"
    fixmate_bam="${work_dir}/${sample}.fixmate.bam"
    position_bam="${work_dir}/${sample}.position.bam"
    new_bam="${work_dir}/${sample}.dedup.bam"

    rg="@RG\\tID:${sample}\\tSM:${sample}\\tLB:${sample}\\tPL:ILLUMINA\\tPU:T2_plus_T3"

    "${bwa}" mem \
        -M \
        -t 2 \
        -R "${rg}" \
        "${reference}" \
        "${fastq1}" \
        "${fastq2}" \
        | "${samtools}" sort \
            -n \
            -@ 1 \
            -m 3G \
            -o "${name_bam}" \
            -

    "${samtools}" fixmate \
        -m \
        -@ 1 \
        "${name_bam}" \
        "${fixmate_bam}"

    "${samtools}" sort \
        -@ 2 \
        -m 3G \
        -o "${position_bam}" \
        "${fixmate_bam}"

    "${samtools}" markdup \
        -r \
        -@ 2 \
        -s \
        -f "${markdup_stats}.new" \
        "${position_bam}" \
        "${new_bam}"

    "${samtools}" index \
        -@ 2 \
        "${new_bam}"

    "${samtools}" quickcheck \
        -v \
        "${new_bam}"

    "${samtools}" flagstat \
        -@ 2 \
        "${new_bam}" \
        > "${flagstat}.new"

    "${samtools}" stats \
        -@ 2 \
        "${new_bam}" \
        > "${stats}.new"

    mv "${new_bam}" "${final_bam}"
    mv "${new_bam}.bai" "${final_bai}"
    mv "${flagstat}.new" "${flagstat}"
    mv "${stats}.new" "${stats}"
    mv "${markdup_stats}.new" "${markdup_stats}"

    printf 'sample\tbam\treference\tinput\n%s\t%s\t%s\tFASTQ_T2_plus_T3\n' \
        "${sample}" \
        "${final_bam}" \
        "${reference}" \
        > "${done_file}.new"

    mv \
        "${done_file}.new" \
        "${done_file}"

    rm -f \
        "${name_bam}" \
        "${fixmate_bam}" \
        "${position_bam}"

    rmdir "${work_dir}" 2>/dev/null || true
fi

printf 'MERGED_SIDE_ALIGNMENT_COMPLETE\tsample=%s\tbam=%s\n' \
    "${sample}" \
    "${final_bam}"
