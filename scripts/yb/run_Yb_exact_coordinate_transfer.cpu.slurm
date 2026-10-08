#!/usr/bin/env bash
#SBATCH --job-name=YbLift
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_exact_liftover.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_exact_liftover.%j.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${YB_CODE_DIR:?YB_CODE_DIR must be exported}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"
SAMTOOLS="$HOME/.local/share/mamba/envs/colias_w_popgen/bin/samtools"

ROOT="$PROJECT_ROOT/analysis/yb_exact_coordinate_transfer"
REFDIR="$ROOT/reference"
SEQDIR="$ROOT/sequences"
RESULTS="$ROOT/results"
HM25_RESULTS="$PROJECT_ROOT/analysis/yb_Hmel2_5_normalized_coverage/results"
HMROS_PATH_FILE="$PROJECT_ROOT/analysis/genomewide_mixing_input/metadata/reference_fasta.path"

mkdir -p "$PROJECT_ROOT/logs" "$REFDIR" "$SEQDIR" "$RESULTS"

"$PYTHON" "$CODE_DIR/analyze_Yb_Hmel2_5_troughs.py" \
  --long "$HM25_RESULTS/Yb_Hmel2.5_Figure8_normalized_50bp.long.tsv.gz" \
  --side-means "$HM25_RESULTS/Yb_Hmel2.5_Figure8_normalized_50bp.side_means.tsv.gz" \
  --output-dir "$RESULTS"

MINIMAP2=""
for x in \
  "$HOME/.local/share/mamba/envs/colias_w_popgen/bin/minimap2" \
  "$HOME/.local/share/mamba/envs/gynandrocalling/bin/minimap2"
do
  [[ -x "$x" ]] && MINIMAP2="$x" && break
done
if [[ -z "$MINIMAP2" ]] && command -v minimap2 >/dev/null 2>&1; then
  MINIMAP2=$(command -v minimap2)
fi
[[ -x "$MINIMAP2" ]] || { echo "ERROR: minimap2 unavailable" >&2; exit 1; }

HM25="$REFDIR/hmel2_5.chromosomes.fa"
if [[ ! -s "$HM25" ]]; then
  LOCAL=""
  while IFS= read -r f; do
    case "$f" in *.gz) continue ;; esac
    if grep -m1 -q '^>Hmel215003o\([[:space:]]\|$\)' "$f"; then
      LOCAL="$f"
      break
    fi
  done < <(
    find "$PROJECT_ROOT" -maxdepth 7 -type f \
      \( -iname '*Hmel2.5*.fa' -o -iname '*Hmel2.5*.fasta' -o \
         -iname '*Hmel2_5*.fa' -o -iname '*Hmel2_5*.fasta' \) \
      -print 2>/dev/null
  )

  if [[ -n "$LOCAL" ]]; then
    echo "REUSING_LOCAL_HMEL2_5=$LOCAL"
    ln -sfn "$LOCAL" "$HM25"
  else
    echo "ERROR: no Hmel2.5 FASTA found under $PROJECT_ROOT" >&2
    echo "TROUGH_ANALYSIS_IS_COMPLETE; coordinate transfer still needs the Hmel2.5 FASTA." >&2
    exit 2
  fi
fi

rm -f "$HM25.fai"
"$SAMTOOLS" faidx "$HM25"
HM25_LEN=$(awk '$1=="Hmel215003o" {print $2; exit}' "$HM25.fai")
[[ "$HM25_LEN" == "10168004" ]] || {
  echo "ERROR: expected Hmel215003o length 10168004; got $HM25_LEN" >&2
  exit 1
}

HMROS=$(head -n1 "$HMROS_PATH_FILE")
[[ -s "$HMROS" && -s "$HMROS.fai" ]] || {
  echo "ERROR: HmRos reference/index missing: $HMROS" >&2
  exit 1
}
[[ "$HMROS" == *"/GCF_056956725.1/"* ]] || {
  echo "ERROR: active reference is not GCF_056956725.1: $HMROS" >&2
  exit 1
}

OLD_START=1350000
OLD_END=1600000
NEW_START=1450000
NEW_END=1900000
OLD_FA="$SEQDIR/Hmel2.5_Yb.fa"
NEW_FA="$SEQDIR/HmRos_Yb.fa"

"$SAMTOOLS" faidx "$HM25" "Hmel215003o:${OLD_START}-${OLD_END}" > "$OLD_FA"
"$SAMTOOLS" faidx "$HMROS" "NC_140403.1:${NEW_START}-${NEW_END}" > "$NEW_FA"

PAF="$RESULTS/Hmel2.5_to_HmRos_Yb.paf"
"$MINIMAP2" -x asm5 -c --cs=long "$NEW_FA" "$OLD_FA" > "$PAF"

[[ -s "$PAF" ]] || { echo "ERROR: minimap2 produced no alignment" >&2; exit 1; }

"$PYTHON" "$CODE_DIR/map_Hmel2_5_to_HmRos.py" \
  --paf "$PAF" \
  --old-region-start "$OLD_START" \
  --new-region-start "$NEW_START" \
  --output-dir "$RESULTS"

echo "YB_EXACT_LIFTOVER_COMPLETE"
