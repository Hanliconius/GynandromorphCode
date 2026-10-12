#!/usr/bin/env bash
#SBATCH --job-name=YbHmRosRef
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_hmros_reference.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_hmros_reference.%j.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
OUT="$PROJECT_ROOT/analysis/yb_raw_from_scratch/reference"
ACCESSION="GCF_056956725.1"
ZIP="$OUT/${ACCESSION}.dataset.zip"
TMP="$OUT/${ACCESSION}.dataset"
FASTA="$OUT/HmRos_GCF_056956725.1_genomic.fna"

module load BWA/0.7.17
module load samtools
command -v curl >/dev/null
command -v unzip >/dev/null
command -v bwa >/dev/null
command -v samtools >/dev/null
mkdir -p "$OUT"

if [[ ! -s "$FASTA" ]]; then
    curl --fail --location --retry 4 \
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${ACCESSION}/download?include_annotation_type=GENOME_FASTA" \
        --output "$ZIP"
    rm -rf "$TMP"
    mkdir -p "$TMP"
    unzip -q "$ZIP" -d "$TMP"
    mapfile -t FASTAS < <(find "$TMP" -type f -name '*_genomic.fna')
    [[ ${#FASTAS[@]} -eq 1 ]] || {
        echo "ERROR: expected one genomic FASTA for $ACCESSION, found ${#FASTAS[@]}" >&2
        exit 2
    }
    cp "${FASTAS[0]}" "$FASTA"
fi

[[ -s "$FASTA" ]] || { echo "ERROR: reference FASTA is empty: $FASTA" >&2; exit 2; }
if [[ ! -s "$FASTA.fai" ]]; then samtools faidx "$FASTA"; fi
if [[ ! -s "$FASTA.bwt" ]]; then bwa index "$FASTA"; fi
grep -q '^NC_140403.1[[:space:]]' "$FASTA.fai" || {
    echo "ERROR: expected Yb scaffold NC_140403.1 absent from $FASTA" >&2
    exit 2
}
echo "YB_HMROS_REFERENCE_READY\t$ACCESSION\t$FASTA"
