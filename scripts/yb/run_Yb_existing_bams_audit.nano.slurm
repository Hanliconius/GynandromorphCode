#!/usr/bin/env bash
#SBATCH --job-name=yb_bam_audit
#SBATCH --partition=nano
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_existing_BAM_audit.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_existing_BAM_audit.%j.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
OUT_DIR="$PROJECT_ROOT/analysis/yb_legacy_replication/rebuild_inventory"
SAMTOOLS=${SAMTOOLS:-$HOME/.local/share/mamba/envs/gynandrocalling/bin/samtools}
CHROM=Hmel215003o
START=1380000
END=1600000

[[ -x "$SAMTOOLS" ]] || { echo "ERROR: samtools not executable: $SAMTOOLS" >&2; exit 2; }
mkdir -p "$OUT_DIR"

printf 'sample\tbam\tindex_stats\treads_overlapping_window\tmean_depth_window\n' \
  > "$OUT_DIR/existing_BAM_Hmel2.5_compatibility.tsv"

for sample in TN_Left_T2 TN_Right_T2 TN_Left_T3 TN_Right_T3; do
    bam="$PROJECT_ROOT/analysis/genomewide_mixing_input/alignment/libraries/${sample}.dedup.bam"
    [[ -s "$bam" ]] || { echo "ERROR: missing BAM: $bam" >&2; exit 3; }

    stats=$("$SAMTOOLS" idxstats "$bam" | awk -F '\t' -v c="$CHROM" \
      '$1==c {printf "length=%s; mapped=%s; unmapped=%s", $2, $3, $4; found=1} END {if (!found) print "CONTIG_ABSENT"}')
    reads=$("$SAMTOOLS" view -c "$bam" "$CHROM:$START-$END")
    depth=$("$SAMTOOLS" depth -aa -r "$CHROM:$START-$END" "$bam" |
      awk 'BEGIN {sum=0; n=0} {sum += $3; n++} END {if (n) printf "%.6f", sum/n; else printf "NA"}')
    printf '%s\t%s\t%s\t%s\t%s\n' "$sample" "$bam" "$stats" "$reads" "$depth" \
      >> "$OUT_DIR/existing_BAM_Hmel2.5_compatibility.tsv"
done

echo "YB_EXISTING_BAM_AUDIT_COMPLETE"
cat "$OUT_DIR/existing_BAM_Hmel2.5_compatibility.tsv"
