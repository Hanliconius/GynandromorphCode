#!/usr/bin/env bash
#SBATCH --job-name=gyn192sibal
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/192_align_full_sib_pedigree.%A_%a.out
#SBATCH --error=logs/192_align_full_sib_pedigree.%A_%a.err

set -euo pipefail

project_root=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

reference=${project_root}/analysis/genomewide_mixing_input/reference/GCF_056956725.1/GCF_056956725.1_ASM5695672v1_genomic.fna
manifest=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/metadata/initial_full_sib_download_manifest.tsv
fastq_root=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/fastq
alignment_root=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/alignment
temporary_root=${project_root}/analysis/genomewide_mixing_corrected_v1/full_sib_validation/tmp/alignment

bwa=${tool_bin}/bwa
samtools=${tool_bin}/samtools

cd "${project_root}"

mkdir -p \
    "${alignment_root}" \
    "${temporary_root}" \
    logs

task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}

sample=$(
    "${python_bin}" - \
        "${manifest}" \
        "${task_id}" <<'PYEOF'
import csv
import sys

with open(
    sys.argv[1],
    newline="",
) as handle:
    rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

matches = [
    row
    for row in rows
    if int(row["array_id"])
    == int(sys.argv[2])
]

if len(matches) != 1:
    raise RuntimeError(
        "Expected one manifest row for "
        f"task {sys.argv[2]}; found "
        f"{len(matches)}"
    )

print(
    matches[0]["sample_accession"]
)
PYEOF
)

fastq1=${fastq_root}/${sample}/${sample}_1.fastq.gz
fastq2=${fastq_root}/${sample}/${sample}_2.fastq.gz

final_bam=${alignment_root}/${sample}.dedup.bam
final_bai=${final_bam}.bai
completion=${alignment_root}/${sample}.complete.tsv

task_tmp=${temporary_root}/${sample}.${SLURM_ARRAY_JOB_ID}.${task_id}

printf 'STAGE_START\tstage=align_full_sib_pedigree\tjob=%s\ttask=%s\tsample=%s\ttime=%s\n' \
    "${SLURM_ARRAY_JOB_ID}" \
    "${task_id}" \
    "${sample}" \
    "$(date --iso-8601=seconds)"

reuse=0

if [[ \
    -s "${completion}" && \
    -s "${final_bam}" && \
    -s "${final_bai}" \
]]
then
    if "${samtools}" quickcheck -v "${final_bam}"
    then
        reuse=1
    fi
fi

if [[ "${reuse}" -eq 1 ]]
then
    printf 'REUSE\tsample=%s\tbam=%s\n' \
        "${sample}" \
        "${final_bam}"
else
    input_status=0

    if [[ ! -s "${fastq1}" || ! -s "${fastq2}" ]]
    then
        printf 'ERROR: missing paired FASTQs for %s\n' \
            "${sample}" >&2
        input_status=1
    fi

    if [[ \
        ! -s "${reference}" || \
        ! -s "${reference}.bwt" || \
        ! -s "${reference}.fai" \
    ]]
    then
        printf 'ERROR: reference or indices are missing.\n' >&2
        input_status=1
    fi

    if [[ "${input_status}" -ne 0 ]]
    then
        false
    fi

    mkdir -p "${task_tmp}"

    rm -f \
        "${final_bam}.new" \
        "${final_bai}.new" \
        "${alignment_root}/${sample}.markdup.stats.txt.new" \
        "${alignment_root}/${sample}.dedup.flagstat.txt.new"

    read_group=$(
        printf '@RG\tID:%s\tSM:%s\tPL:ILLUMINA\tLB:PRJEB11288\tPU:%s' \
            "${sample}" \
            "${sample}" \
            "${sample}"
    )

    "${bwa}" mem \
        -t "${SLURM_CPUS_PER_TASK}" \
        -R "${read_group}" \
        "${reference}" \
        "${fastq1}" \
        "${fastq2}" \
        | "${samtools}" sort \
            -n \
            -@ 2 \
            -m 2G \
            -o "${task_tmp}/${sample}.name.bam" \
            -

    "${samtools}" fixmate \
        -m \
        -@ 2 \
        "${task_tmp}/${sample}.name.bam" \
        "${task_tmp}/${sample}.fixmate.bam"

    "${samtools}" sort \
        -@ 2 \
        -m 2G \
        -o "${task_tmp}/${sample}.position.bam" \
        "${task_tmp}/${sample}.fixmate.bam"

    "${samtools}" markdup \
        -r \
        -@ 2 \
        -s \
        -f "${alignment_root}/${sample}.markdup.stats.txt.new" \
        "${task_tmp}/${sample}.position.bam" \
        "${final_bam}.new"

    "${samtools}" index \
        -@ 2 \
        "${final_bam}.new" \
        "${final_bai}.new"

    "${samtools}" quickcheck -v \
        "${final_bam}.new"

    observed_sample=$(
        "${samtools}" view \
            -H \
            "${final_bam}.new" \
            | awk -F '\t' '
                $1 == "@RG" {
                    for (
                        field = 1;
                        field <= NF;
                        field++
                    ) {
                        if (
                            $field ~ /^SM:/
                        ) {
                            sub(
                                /^SM:/,
                                "",
                                $field
                            )
                            print $field
                        }
                    }
                }
            ' \
            | sort -u
    )

    if [[ "${observed_sample}" != "${sample}" ]]
    then
        printf 'ERROR: BAM sample is %s; expected %s.\n' \
            "${observed_sample}" \
            "${sample}" >&2
        false
    fi

    "${samtools}" flagstat \
        -@ 2 \
        "${final_bam}.new" \
        > "${alignment_root}/${sample}.dedup.flagstat.txt.new"

    mv \
        "${final_bam}.new" \
        "${final_bam}"

    mv \
        "${final_bai}.new" \
        "${final_bai}"

    mv \
        "${alignment_root}/${sample}.markdup.stats.txt.new" \
        "${alignment_root}/${sample}.markdup.stats.txt"

    mv \
        "${alignment_root}/${sample}.dedup.flagstat.txt.new" \
        "${alignment_root}/${sample}.dedup.flagstat.txt"

    {
        printf 'sample\tbam\tbam_bytes\tcompleted_at\n'

        printf '%s\t%s\t%s\t%s\n' \
            "${sample}" \
            "${final_bam}" \
            "$(wc -c < "${final_bam}")" \
            "$(date --iso-8601=seconds)"
    } > "${completion}.new"

    mv \
        "${completion}.new" \
        "${completion}"

    rm -f \
        "${task_tmp}/${sample}.name.bam" \
        "${task_tmp}/${sample}.fixmate.bam" \
        "${task_tmp}/${sample}.position.bam"

    printf 'COMPLETE\tsample=%s\tbam_bytes=%s\n' \
        "${sample}" \
        "$(wc -c < "${final_bam}")"
fi

printf 'STAGE_END\tstage=align_full_sib_pedigree\tstatus=0\tsample=%s\ttime=%s\n' \
    "${sample}" \
    "$(date --iso-8601=seconds)"
