#!/usr/bin/env bash
#SBATCH --job-name=YbNormCov
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_normalized_coverage.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_normalized_coverage.%j.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${YB_CODE_DIR:?YB_CODE_DIR must be exported by submit helper}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

OUT_ROOT="$PROJECT_ROOT/analysis/yb_normalized_coverage"
BW_DIR="$OUT_ROOT/krista_upstream_bigwig"
BG_DIR="$OUT_ROOT/bedgraph"
RESULT_DIR="$OUT_ROOT/results"
TOOLS_DIR="$OUT_ROOT/tools"

mkdir -p     "$PROJECT_ROOT/logs"     "$BW_DIR"     "$BG_DIR"     "$RESULT_DIR"     "$TOOLS_DIR"

UPSTREAM_COMMIT="99376deebc25f4e0bb7aa9927c7969c693591ef4"
UPSTREAM_BASE="https://raw.githubusercontent.com/Krista-Pipho/gynandromorph_2024/$UPSTREAM_COMMIT/R_Markdown"

samples=(
    TN2L_S1
    TN3L_S3
    TN2R_S2
    TN3R_S4
)

echo "YB_NORMALIZED_COVERAGE_START"
echo "job=${SLURM_JOB_ID}"
echo "host=$(hostname)"
echo "upstream_commit=$UPSTREAM_COMMIT"
echo "coordinate_system=Hmel2.5/Hmel215003o"
echo "analysis_region=Hmel215003o:1380000-1600000"
echo

for sample in "${samples[@]}"
do
    bw="$BW_DIR/${sample}.dedup.bw"

    if [[ ! -s "$bw" ]]
    then
        echo "DOWNLOADING $sample"
        curl             --fail             --location             --retry 3             --output "$bw.part"             "$UPSTREAM_BASE/${sample}.dedup.bw"

        mv "$bw.part" "$bw"
    fi

    if [[ ! -s "$bw" ]]
    then
        echo "ERROR: empty bigWig after download: $bw" >&2
        exit 1
    fi

    echo -e "BIGWIG\t$sample\t$(stat -c '%s' "$bw")\t$bw"
done

converter=""

for candidate in     "$(command -v bigWigToBedGraph 2>/dev/null || true)"     "$HOME/.local/share/mamba/envs/gynandrocalling/bin/bigWigToBedGraph"     "$HOME/.local/share/mamba/envs/colias_w_popgen/bin/bigWigToBedGraph"
do
    if [[ -n "$candidate" && -x "$candidate" ]]
    then
        converter="$candidate"
        break
    fi
done

if [[ -z "$converter" ]]
then
    converter="$TOOLS_DIR/bigWigToBedGraph"

    if [[ ! -x "$converter" ]]
    then
        echo "DOWNLOADING_UCSC_BIGWIGTOBEDGRAPH"
        curl             --fail             --location             --retry 3             --output "$converter.part"             "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph"

        chmod +x "$converter.part"
        mv "$converter.part" "$converter"
    fi
fi

echo "BIGWIG_CONVERTER=$converter"
"$converter" 2>&1 | head -n 2 || true

for sample in "${samples[@]}"
do
    bw="$BW_DIR/${sample}.dedup.bw"
    bg="$BG_DIR/${sample}.dedup.Yb_region.bedGraph"

    "$converter"         -chrom=Hmel215003o         -start=1380000         -end=1600000         "$bw"         "$bg.tmp"

    if [[ ! -s "$bg.tmp" ]]
    then
        echo "ERROR: no Yb-region coverage extracted for $sample" >&2
        exit 1
    fi

    mv "$bg.tmp" "$bg"

    echo -e "BEDGRAPH\t$sample\t$(wc -l < "$bg")\t$bg"
done

"$PYTHON"     "$CODE_DIR/plot_Yb_normalized_coverage.py"     --bedgraph-dir "$BG_DIR"     --output-dir "$RESULT_DIR"

echo "YB_NORMALIZED_COVERAGE_FINISHED"
