#!/usr/bin/env bash
#SBATCH --job-name=YbHm25
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:25:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_Hmel2_5_normalized.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Yb_Hmel2_5_normalized.%j.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${YB_CODE_DIR:?YB_CODE_DIR must be exported}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

ROOT="$PROJECT_ROOT/analysis/yb_Hmel2_5_normalized_coverage"
BW_DIR="$ROOT/krista_upstream_bigwig"
BG_DIR="$ROOT/Hmel215003o_bedgraph"
META="$ROOT/metadata"
RESULTS="$ROOT/results"
TOOLS="$ROOT/tools"

mkdir -p "$PROJECT_ROOT/logs" "$BW_DIR" "$BG_DIR" "$META" "$RESULTS" "$TOOLS"

UPSTREAM_COMMIT="99376deebc25f4e0bb7aa9927c7969c693591ef4"
UPSTREAM_BASE="https://raw.githubusercontent.com/Krista-Pipho/gynandromorph_2024/$UPSTREAM_COMMIT/R_Markdown"
CHROM="Hmel215003o"

samples=(
    TN2L_S1
    TN3L_S3
    TN2R_S2
    TN3R_S4
)

get_ucsc_tool() {
    local name=$1
    local path="$TOOLS/$name"

    if [[ ! -x "$path" ]]
    then
        /usr/bin/curl             --fail             --location             --retry 5             --output "$path.part"             "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/$name"

        chmod +x "$path.part"
        mv "$path.part" "$path"
    fi

    printf '%s\n' "$path"
}

BIGWIG_TO_BG=$(get_ucsc_tool bigWigToBedGraph)
BIGWIG_INFO=$(get_ucsc_tool bigWigInfo)

echo "YB_HMEL2_5_NORMALIZED_START"
echo "job=${SLURM_JOB_ID}"
echo "host=$(hostname)"
echo "upstream_commit=$UPSTREAM_COMMIT"
echo "reference=Hmel2.5"
echo "scaffold=$CHROM"
echo "normalization=local 50bp coverage / sample mean full-scaffold coverage"
echo

STATS="$META/Hmel215003o_scaffold_means.tsv"
printf 'sample\tscaffold_length\tscaffold_mean_depth\tbedgraph\n' > "$STATS.tmp"

for sample in "${samples[@]}"
do
    bw="$BW_DIR/${sample}.dedup.bw"

    if [[ ! -s "$bw" ]]
    then
        echo "DOWNLOADING_BIGWIG $sample"
        /usr/bin/curl             --fail             --location             --retry 5             --output "$bw.part"             "$UPSTREAM_BASE/${sample}.dedup.bw"
        mv "$bw.part" "$bw"
    fi

    chrom_length=$(
        "$BIGWIG_INFO" -chroms "$bw"         | awk -v chr="$CHROM" '$1==chr {print $3; exit}'
    )

    if [[ ! "$chrom_length" =~ ^[0-9]+$ || "$chrom_length" -le 0 ]]
    then
        echo "ERROR: could not determine $CHROM length in $bw" >&2
        "$BIGWIG_INFO" -chroms "$bw" >&2
        exit 1
    fi

    bg="$BG_DIR/${sample}.Hmel215003o.bedGraph"

    if [[ ! -s "$bg" ]]
    then
        "$BIGWIG_TO_BG"             -chrom="$CHROM"             -start=0             -end="$chrom_length"             "$bw"             "$bg.tmp"

        [[ -s "$bg.tmp" ]] || {
            echo "ERROR: empty scaffold bedGraph for $sample" >&2
            exit 1
        }

        mv "$bg.tmp" "$bg"
    fi

    scaffold_mean=$(
        awk             -v len="$chrom_length"             'BEGIN {num=0}
             {
                 w=$3-$2
                 if (w > 0) num += w*$4
             }
             END {
                 if (len <= 0) exit 1
                 printf "%.10f", num/len
             }'             "$bg"
    )

    if ! awk -v x="$scaffold_mean" 'BEGIN {exit !(x>0)}'
    then
        echo "ERROR: invalid full-scaffold mean for $sample: $scaffold_mean" >&2
        exit 1
    fi

    printf '%s\t%s\t%s\t%s\n'         "$sample" "$chrom_length" "$scaffold_mean" "$bg"         >> "$STATS.tmp"

    echo -e "SCAFFOLD_MEAN\t$sample\tlength=$chrom_length\tmean=$scaffold_mean"
done

mv "$STATS.tmp" "$STATS"

"$PYTHON" "$CODE_DIR/plot_Yb_Hmel2_5_Figure8_normalized.py"     --bedgraph-dir "$BG_DIR"     --scaffold-stats "$STATS"     --output-dir "$RESULTS"

echo "YB_HMEL2_5_NORMALIZED_FINISHED"
