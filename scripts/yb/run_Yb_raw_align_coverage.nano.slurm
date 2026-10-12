#!/usr/bin/env bash
#SBATCH --job-name=YbRawMap
#SBATCH --partition=cpu
#SBATCH --array=0-3%4
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_raw_map.%A_%a.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_raw_map.%A_%a.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${YB_CODE_DIR:?YB_CODE_DIR must be exported}
ANALYSIS="$PROJECT_ROOT/analysis/yb_raw_from_scratch"
REF="$ANALYSIS/reference/HmRos_GCF_056956725.1_genomic.fna"
FASTQ="$ANALYSIS/fastq"
OUT="$ANALYSIS/HmRos"
METRICS="$OUT/metrics"
SAMPLES=(TN2L_S1 TN2R_S2 TN3L_S3 TN3R_S4)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

module load BWA/0.7.17
module load samtools
module load Java/1.8.0_60 || true
module load deepTools || module load deeptools || true
module load picard || module load Picard || true
command -v bwa >/dev/null
command -v samtools >/dev/null
command -v bamCoverage >/dev/null || {
    echo "ERROR: bamCoverage is not available in PATH; load the deepTools environment." >&2
    exit 2
}

R1="$FASTQ/${SAMPLE}_R1.fastq.gz"
R2="$FASTQ/${SAMPLE}_R2.fastq.gz"
RAW_BAM="$OUT/${SAMPLE}.raw.bam"
DEDUP_BAM="$OUT/${SAMPLE}.dedup.bam"
BW="$OUT/${SAMPLE}.dedup.50bp.bw"
PICARD_JAR=${PICARD_JAR:-"$PROJECT_ROOT/picard.jar"}

for f in "$REF" "$R1" "$R2"; do
    [[ -s "$f" ]] || { echo "ERROR: required input missing: $f" >&2; exit 2; }
done
mkdir -p "$OUT" "$METRICS"

RG=$(printf '@RG\\tID:%s\\tSM:%s\\tPL:ILLUMINA\\tLB:%s' "$SAMPLE" "$SAMPLE" "$SAMPLE")
bwa mem -M -t "$SLURM_CPUS_PER_TASK" -R "$RG" "$REF" "$R1" "$R2" \
    | samtools sort -@ 4 -m 2G -o "$RAW_BAM" -
samtools index -@ 4 "$RAW_BAM"

if command -v picard >/dev/null 2>&1; then
    picard MarkDuplicates \
        I="$RAW_BAM" O="$DEDUP_BAM" M="$METRICS/${SAMPLE}.picard_metrics.txt" \
        VALIDATION_STRINGENCY=LENIENT
elif [[ -s "$PICARD_JAR" ]]; then
    java -Xmx48g -jar "$PICARD_JAR" MarkDuplicates \
        I="$RAW_BAM" O="$DEDUP_BAM" M="$METRICS/${SAMPLE}.picard_metrics.txt" \
        VALIDATION_STRINGENCY=LENIENT
else
    echo "ERROR: Picard not found; set PICARD_JAR or load the Picard module." >&2
    exit 2
fi
samtools index -@ 4 "$DEDUP_BAM"

# Krista's bigWig recipe uses bamCoverage defaults (no MAPQ filter,
# no normalization, and no duplicate exclusion) with its default 50-bp bins.
bamCoverage --bam "$DEDUP_BAM" --outFileName "$BW" \
    --outFileFormat bigwig --binSize 50 --numberOfProcessors "$SLURM_CPUS_PER_TASK"

samtools flagstat "$DEDUP_BAM" > "$OUT/${SAMPLE}.flagstat.txt"
echo "YB_RAW_SAMPLE_COMPLETE\t$SAMPLE\t$BW"
