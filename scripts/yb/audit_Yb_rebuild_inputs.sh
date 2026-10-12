#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_ROOT=${GYNANDRO_CODE_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo}
OUT="$PROJECT_ROOT/analysis/yb_legacy_replication/rebuild_inventory"
mkdir -p "$OUT"

READS="$OUT/raw_read_files.tsv"
ARRAYS="$OUT/alignment_call_array_scripts.tsv"
REFERENCES="$OUT/Hmel_reference_candidates.tsv"
BAMS="$OUT/current_BAM_candidates.tsv"

printf 'sample\tbiosample\tside\tleg\tread1\tread2\n' > "$OUT/Yb_samples.tsv"
printf 'TN2L_S1\tSAMN60580458\tLeft\tT2\tTN2L_S1_1.fastq.gz\tTN2L_S1_2.fastq.gz\n' >> "$OUT/Yb_samples.tsv"
printf 'TN2R_S2\tSAMN60580459\tRight\tT2\tTN2R_S2_1.fastq.gz\tTN2R_S2_2.fastq.gz\n' >> "$OUT/Yb_samples.tsv"
printf 'TN3L_S3\tSAMN60580460\tLeft\tT3\tTN3L_S3_1.fastq.gz\tTN3L_S3_2.fastq.gz\n' >> "$OUT/Yb_samples.tsv"
printf 'TN3R_S4\tSAMN60580461\tRight\tT3\tTN3R_S4_1.fastq.gz\tTN3R_S4_2.fastq.gz\n' >> "$OUT/Yb_samples.tsv"

: > "$READS"
: > "$ARRAYS"
: > "$REFERENCES"
: > "$BAMS"

if [[ -d "$PROJECT_ROOT" ]]; then
    find "$PROJECT_ROOT" -maxdepth 9 \
        \( -type d \( -name .git -o -name logs -o -name tmp -o -name cache \) -prune \) -o \
        \( -type f \( \
            -name 'TN2L_S1_1.fastq.gz' -o -name 'TN2L_S1_2.fastq.gz' -o \
            -name 'TN2R_S2_1.fastq.gz' -o -name 'TN2R_S2_2.fastq.gz' -o \
            -name 'TN3L_S3_1.fastq.gz' -o -name 'TN3L_S3_2.fastq.gz' -o \
            -name 'TN3R_S4_1.fastq.gz' -o -name 'TN3R_S4_2.fastq.gz' \
        \) -print \) > "$READS"

    find "$PROJECT_ROOT" -maxdepth 9 -type f \
        \( -iname '*Hmel2.5*.fa' -o -iname '*Hmel2.5*.fasta' -o \
           -iname '*Hmel2.5*.fna' -o -iname '*Hmel*genome*.fa' -o \
           -iname '*Hmel*genome*.fasta' -o -iname '*Hmel*genome*.fna' \) \
        -print > "$REFERENCES"

    BAM_DIR="$PROJECT_ROOT/analysis/genomewide_mixing_input/alignment/libraries"
    if [[ -d "$BAM_DIR" ]]; then
        find "$BAM_DIR" -maxdepth 2 -type f -name '*.bam' \
            \( -iname '*TN*Left*T2*' -o -iname '*TN*Left*T3*' -o \
               -iname '*TN*Right*T2*' -o -iname '*TN*Right*T3*' -o \
               -iname '*TN2L_S1*' -o -iname '*TN2R_S2*' -o \
               -iname '*TN3L_S3*' -o -iname '*TN3R_S4*' \) \
            -print > "$BAMS"
    fi
fi

for root in "$PROJECT_ROOT/scripts" "$CODE_ROOT/scripts"; do
    [[ -d "$root" ]] || continue
    find "$root" -type f \( -name '*.sh' -o -name '*.slurm' -o -name '*.sbatch' \) -print0
done | while IFS= read -r -d '' file; do
    if grep -qiE 'SBATCH.*--array|sbatch.*--array|--array=' "$file" &&
       grep -qiE 'align|mapping|bwa|bowtie|hisat|fastq|variant|call|bcftools|gatk' "$file"; then
        printf '%s\n' "$file"
    fi
done | sort -u > "$ARRAYS"

printf 'YB_REBUILD_INPUT_AUDIT_COMPLETE\n'
printf 'READ_FILES\t%s\t%s\n' "$(wc -l < "$READS" | tr -d ' ')" "$READS"
printf 'ARRAY_SCRIPTS\t%s\t%s\n' "$(wc -l < "$ARRAYS" | tr -d ' ')" "$ARRAYS"
printf 'Hmel_REFERENCES\t%s\t%s\n' "$(wc -l < "$REFERENCES" | tr -d ' ')" "$REFERENCES"
printf 'CURRENT_BAMS\t%s\t%s\n' "$(wc -l < "$BAMS" | tr -d ' ')" "$BAMS"
