#!/usr/bin/env bash
#SBATCH --job-name=YbHmRos
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_hmros_normalized.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_hmros_normalized.%j.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${YB_CODE_DIR:?YB_CODE_DIR must be exported}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"
SAMTOOLS="$HOME/.local/share/mamba/envs/colias_w_popgen/bin/samtools"

ROOT="$PROJECT_ROOT/analysis/yb_hmros_normalized_coverage"
META="$ROOT/metadata"
DEPTH="$ROOT/depth_250bp"
RESULTS="$ROOT/results"

REF_DIR="$PROJECT_ROOT/analysis/genomewide_mixing_input/reference/GCF_056956725.1"
REF_PATH_FILE="$PROJECT_ROOT/analysis/genomewide_mixing_input/metadata/reference_fasta.path"
REPORT="$REF_DIR/GCF_056956725.1_ASM5695672v1_assembly_report.txt"
GFF="$REF_DIR/GCF_056956725.1_ASM5695672v1_genomic.gff.gz"

MANIFEST="$PROJECT_ROOT/analysis/genomewide_mixing_input/metadata/resolved_download_manifest.tsv"
LIBDIR="$PROJECT_ROOT/analysis/genomewide_mixing_input/alignment/libraries"

mkdir -p "$PROJECT_ROOT/logs" "$META" "$DEPTH" "$RESULTS"

for f in "$REF_PATH_FILE" "$REPORT" "$MANIFEST"
do
    [[ -s "$f" ]] || {
        echo "ERROR: missing required input: $f" >&2
        exit 1
    }
done

REFERENCE=$(head -n 1 "$REF_PATH_FILE")

[[ -s "$REFERENCE" && -s "$REFERENCE.fai" ]] || {
    echo "ERROR: HmRos reference or index missing: $REFERENCE" >&2
    exit 1
}

if [[ "$REFERENCE" != *"/GCF_056956725.1/"* ]]
then
    echo "ERROR: active reference is not GCF_056956725.1: $REFERENCE" >&2
    exit 1
fi

if [[ ! -s "$GFF" ]]
then
    URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/056/956/725/GCF_056956725.1_ASM5695672v1/GCF_056956725.1_ASM5695672v1_genomic.gff.gz"
    echo "DOWNLOADING_REFSEQ_GFF $URL"
    /usr/bin/curl --fail --location --retry 5 --output "$GFF.part" "$URL"
    gzip -t "$GFF.part"
    mv "$GFF.part" "$GFF"
fi

COORDS="$META/Yb_hmros_coordinates.tsv"

"$PYTHON" "$CODE_DIR/resolve_Yb_hmros_coordinates.py"     --assembly-report "$REPORT"     --gff "$GFF"     --reference-fai "$REFERENCE.fai"     --output "$COORDS"

SEQID=$(awk -F '\t' 'NR==2 {print $4}' "$COORDS")
START=$(awk -F '\t' 'NR==2 {print $9}' "$COORDS")
END=$(awk -F '\t' 'NR==2 {print $10}' "$COORDS")

[[ -n "$SEQID" && "$START" =~ ^[0-9]+$ && "$END" =~ ^[0-9]+$ ]] || {
    echo "ERROR: failed to parse resolved HmRos region" >&2
    cat "$COORDS" >&2
    exit 1
}

echo "RESOLVED_REGION=${SEQID}:${START}-${END}"

SAMPLES="$META/samples.tsv"
printf 'sample\trole\tside\tleg\tbam\tstats\n' > "$SAMPLES.tmp"

while IFS=$'\t' read -r     array_id source_accession resolved_run library_id role side leg merge_group layout fastq1 fastq2
do
    [[ "$array_id" == "array_id" ]] && continue

    bam="$LIBDIR/${library_id}.dedup.bam"
    stats="$LIBDIR/${library_id}.dedup.stats.txt"

    [[ -s "$bam" && -s "$bam.bai" && -s "$stats" ]] || {
        echo "ERROR: missing BAM/index/stats for $library_id" >&2
        exit 1
    }

    printf '%s\t%s\t%s\t%s\t%s\t%s\n'         "$library_id" "$role" "$side" "$leg" "$bam" "$stats"         >> "$SAMPLES.tmp"
done < "$MANIFEST"

mv "$SAMPLES.tmp" "$SAMPLES"

CONTROL_N=$(awk -F '\t' 'NR>1 && $2=="control" {n++} END {print n+0}' "$SAMPLES")
FOCAL_N=$(awk -F '\t' 'NR>1 && $2=="focal" {n++} END {print n+0}' "$SAMPLES")

[[ "$CONTROL_N" -eq 18 && "$FOCAL_N" -eq 4 ]] || {
    echo "ERROR: expected 18 controls + 4 focal libraries; got controls=$CONTROL_N focal=$FOCAL_N" >&2
    exit 1
}

ALL_BINS="$META/Yb_hmros_250bp_raw_bins.tsv"
printf 'sample\tseqid\tstart\tend\tmidpoint\tmean_depth\tmapped_bases\n' > "$ALL_BINS.tmp"

while IFS=$'\t' read -r sample role side leg bam stats
do
    [[ "$sample" == "sample" ]] && continue

    mapped_bases=$(
        awk -F '\t' '
            $1=="SN" && $2 ~ /^bases mapped \(cigar\):/ {
                gsub(/[^0-9]/,"",$3)
                print $3
                exit
            }
        ' "$stats"
    )

    [[ "$mapped_bases" =~ ^[0-9]+$ && "$mapped_bases" -gt 0 ]] || {
        echo "ERROR: could not parse mapped bases from $stats" >&2
        exit 1
    }

    out="$DEPTH/${sample}.250bp.tsv"

    "$SAMTOOLS" depth         -aa         -q 20         -Q 20         -r "${SEQID}:${START}-${END}"         "$bam"     | awk         -v sample="$sample"         -v seq="$SEQID"         -v region_start="$START"         -v region_end="$END"         -v bin=250         -v mapped_bases="$mapped_bases"         'BEGIN {OFS="\t"}
         {
             idx=int(($2-region_start)/bin)
             sum[idx]+=$3
             n[idx]++
         }
         END {
             for (idx=0; ; idx++) {
                 s=region_start + idx*bin
                 if (s > region_end) break
                 e=s+bin-1
                 if (e > region_end) e=region_end
                 expected=e-s+1
                 mean=(n[idx] ? sum[idx]/expected : 0)
                 midpoint=(s+e)/2
                 print sample,seq,s,e,midpoint,mean,mapped_bases
             }
         }'         > "$out.tmp"

    mv "$out.tmp" "$out"
    cat "$out" >> "$ALL_BINS.tmp"

    echo -e "DEPTH_DONE\t$sample\tmapped_bases=$mapped_bases"
done < "$SAMPLES"

mv "$ALL_BINS.tmp" "$ALL_BINS"

"$PYTHON" "$CODE_DIR/plot_Yb_hmros_panel_normalized_coverage.py"     --bins "$ALL_BINS"     --samples "$SAMPLES"     --coordinates "$COORDS"     --output-dir "$RESULTS"

echo "YB_HMROS_WORKFLOW_COMPLETE"
