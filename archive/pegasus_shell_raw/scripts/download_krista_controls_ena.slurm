#!/bin/bash
#SBATCH --job-name=getKristaE
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/getKristaE_%A_%a.out
#SBATCH --error=logs/getKristaE_%A_%a.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
OUTROOT="${ROOT}/analysis/genomewide_mixing_input/fastq"

SAMPLES=(
    ERR1050907
    ERR1050958
    ERR1050994
    ERR1051012
)

s="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
out="${OUTROOT}/${s}"

mkdir -p "$out"

prefix="${s:0:6}"

# ENA layout for these 7-digit ERR accessions:
# ERR1050907 -> ERR105/007/ERR1050907
suffix="00${s: -1}"

base="https://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${suffix}/${s}"

echo "SAMPLE = ${s}"
echo "BASE   = ${base}"
echo "OUT    = ${out}"
echo

if [[ -s "${out}/${s}_1.fastq.gz" && -s "${out}/${s}_2.fastq.gz" ]]; then
    echo "Already complete; skipping."
    exit 0
fi

for mate in 1 2; do
    url="${base}/${s}_${mate}.fastq.gz"
    dest="${out}/${s}_${mate}.fastq.gz"

    echo "Downloading:"
    echo "$url"

    curl \
        --fail \
        --location \
        --retry 5 \
        --retry-delay 10 \
        --continue-at - \
        --output "$dest" \
        "$url"

    test -s "$dest"
done

echo
echo "===== COMPLETE ====="
ls -lh \
    "${out}/${s}_1.fastq.gz" \
    "${out}/${s}_2.fastq.gz"
