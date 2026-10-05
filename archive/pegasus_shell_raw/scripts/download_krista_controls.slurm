#!/bin/bash
#SBATCH --job-name=getKrista
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/getKrista_%A_%a.out
#SBATCH --error=logs/getKrista_%A_%a.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
OUTROOT="${ROOT}/analysis/genomewide_mixing_input/fastq"

SRA_BIN="$HOME/.local/share/mamba/envs/colias_w_popgen/bin"
export PATH="${SRA_BIN}:${PATH}"

SAMPLES=(
    ERR1050907
    ERR1050958
    ERR1050994
    ERR1051012
)

s="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
out="${OUTROOT}/${s}"

mkdir -p "$out"

echo "SAMPLE       = $s"
echo "OUT          = $out"
echo "prefetch     = $(command -v prefetch)"
echo "fasterq-dump = $(command -v fasterq-dump)"
echo "pigz         = $(command -v pigz)"
echo

# Idempotent
if [[ -s "${out}/${s}_1.fastq.gz" && -s "${out}/${s}_2.fastq.gz" ]]; then
    echo "Already complete; skipping."
    exit 0
fi

tmp="${out}/sra_tmp"
mkdir -p "$tmp"

echo "===== PREFETCH ====="

prefetch "$s" \
    --output-directory "$tmp"

echo
echo "===== FASTERQ-DUMP ====="

fasterq-dump \
    "$tmp/$s/$s.sra" \
    --split-files \
    --threads 4 \
    --temp "$tmp" \
    --outdir "$out"

test -s "${out}/${s}_1.fastq"
test -s "${out}/${s}_2.fastq"

echo
echo "===== COMPRESS ====="

pigz -p 4 \
    "${out}/${s}_1.fastq" \
    "${out}/${s}_2.fastq"

test -s "${out}/${s}_1.fastq.gz"
test -s "${out}/${s}_2.fastq.gz"

rm -rf "$tmp"

echo
echo "===== COMPLETE ====="

ls -lh \
    "${out}/${s}_1.fastq.gz" \
    "${out}/${s}_2.fastq.gz"
