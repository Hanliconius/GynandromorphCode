#!/bin/bash
#SBATCH --job-name=rankDepth
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/rankDepth_%j.out
#SBATCH --error=logs/rankDepth_%j.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
BAMLIST="${ROOT}/analysis/genomewide_mixing_input/metadata/bam_sets_24/all_24_samples.bamlist"
OUTDIR="${ROOT}/analysis/smudgeplot_arrangements/control_depth_screen"

mkdir -p "${OUTDIR}"

OUT="${OUTDIR}/all_24_samples.mean_depth.tsv"

printf "sample\tbam\tmean_depth\n" > "${OUT}"

while IFS= read -r bam; do
    [[ -n "$bam" ]] || continue

    base="$(basename "$bam")"

    sample="${base%.bam}"
    sample="${sample%.cram}"

    # mean depth across reference sequences, weighted by covered sequence length
    samtools coverage "$bam" \
      | awk -v sample="$sample" -v bam="$bam" '
        BEGIN {
            num=0
            den=0
        }
        /^#/ {next}
        {
            len=$3-$2+1
            num += $7 * len
            den += len
        }
        END {
            if (den > 0)
                printf "%s\t%s\t%.6f\n", sample, bam, num/den
        }
      ' >> "${OUT}"

done < "${BAMLIST}"

{
    head -1 "${OUT}"
    tail -n +2 "${OUT}" | sort -t $'\t' -k3,3gr
} > "${OUTDIR}/all_24_samples.mean_depth.ranked.tsv"

echo
echo "TOP SAMPLES"
column -t -s $'\t' "${OUTDIR}/all_24_samples.mean_depth.ranked.tsv" | head -15
