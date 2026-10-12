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
SAMTOOLS=${SAMTOOLS:-$(command -v samtools || true)}

[[ -n "$SAMTOOLS" && -x "$SAMTOOLS" ]] || {
    echo "ERROR: samtools not found; load the samtools module before submitting" >&2
    exit 2
}
mkdir -p "$OUT_DIR"

printf 'sample\tbam\tcontig_count\tcontains_Hmel215003o\theader_dictionary\n' \
  > "$OUT_DIR/existing_BAM_reference_dictionary_summary.tsv"

FASTA="$PROJECT_ROOT/analysis/yb_exact_coordinate_transfer/sequences/Hmel2.5_Yb.fa"
if [[ -s "$FASTA" ]]; then
    awk '/^>/ {print "Hmel2.5_Yb_FASTA_HEADER\t" $0; exit}' "$FASTA"
fi

for sample in TN_Left_T2 TN_Right_T2 TN_Left_T3 TN_Right_T3; do
    bam="$PROJECT_ROOT/analysis/genomewide_mixing_input/alignment/libraries/${sample}.dedup.bam"
    [[ -s "$bam" ]] || { echo "ERROR: missing BAM: $bam" >&2; exit 3; }

    dictionary="$OUT_DIR/${sample}.BAM_reference_dictionary.tsv"
    "$SAMTOOLS" view -H "$bam" | awk -F '\t' '
      $1=="@SQ" {
        sn=""; ln="NA"
        for (i=2; i<=NF; i++) {
          if ($i ~ /^SN:/) {sn=substr($i,4)}
          if ($i ~ /^LN:/) {ln=substr($i,4)}
        }
        if (sn!="") print sn "\t" ln
      }' > "$dictionary"

    count=$(wc -l < "$dictionary" | tr -d ' ')
    if awk -F '\t' '$1=="Hmel215003o" {found=1} END {exit !found}' "$dictionary"; then
        has_target=yes
    else
        has_target=no
    fi
    printf '%s\t%s\t%s\t%s\t%s\n' "$sample" "$bam" "$count" "$has_target" "$dictionary" \
      >> "$OUT_DIR/existing_BAM_reference_dictionary_summary.tsv"

    awk -F '\t' -v s="$sample" 'tolower($1) ~ /hmel|215003|yb|chr/ {print "REFERENCE_MATCH\t" s "\t" $1 "\t" $2}' "$dictionary"
    if ! awk -F '\t' 'tolower($1) ~ /hmel|215003|yb|chr/ {found=1} END {exit !found}' "$dictionary"; then
        head -n 8 "$dictionary" | awk -v s="$sample" '{print "REFERENCE_HEAD\t" s "\t" $0}'
    fi
done

echo "YB_EXISTING_BAM_DICTIONARY_AUDIT_COMPLETE"
cat "$OUT_DIR/existing_BAM_reference_dictionary_summary.tsv"
