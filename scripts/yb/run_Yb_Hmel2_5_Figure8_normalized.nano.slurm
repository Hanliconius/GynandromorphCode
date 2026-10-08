#!/usr/bin/env bash
#SBATCH --job-name=YbHm25
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
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
PYBIGWIG_DIR="$TOOLS/pybigwig"

mkdir -p \
    "$PROJECT_ROOT/logs" \
    "$BW_DIR" \
    "$BG_DIR" \
    "$META" \
    "$RESULTS" \
    "$TOOLS"

UPSTREAM_COMMIT="99376deebc25f4e0bb7aa9927c7969c693591ef4"
UPSTREAM_BASE="https://raw.githubusercontent.com/Krista-Pipho/gynandromorph_2024/$UPSTREAM_COMMIT/R_Markdown"
CHROM="Hmel215003o"

samples=(
    TN2L_S1
    TN3L_S3
    TN2R_S2
    TN3R_S4
)

echo "YB_HMEL2_5_NORMALIZED_START"
echo "job=${SLURM_JOB_ID}"
echo "host=$(hostname)"
echo "upstream_commit=$UPSTREAM_COMMIT"
echo "reference=Hmel2.5"
echo "scaffold=$CHROM"
echo "normalization=local 50bp coverage / sample mean full-scaffold coverage"
echo

# Pegasus has an older host glibc than current UCSC precompiled bigWig tools.
# Use pyBigWig instead. Keep it local to this analysis so no shared conda
# environment is modified.
if ! PYTHONPATH="$PYBIGWIG_DIR${PYTHONPATH:+:$PYTHONPATH}" \
    "$PYTHON" -c 'import pyBigWig' >/dev/null 2>&1
then
    echo "INSTALLING_LOCAL_PYBIGWIG"

    rm -rf "$PYBIGWIG_DIR.tmp"
    mkdir -p "$PYBIGWIG_DIR.tmp"

    "$PYTHON" -m pip install \
        --disable-pip-version-check \
        --no-deps \
        --target "$PYBIGWIG_DIR.tmp" \
        "pyBigWig==0.3.24"

    rm -rf "$PYBIGWIG_DIR"
    mv "$PYBIGWIG_DIR.tmp" "$PYBIGWIG_DIR"
fi

export PYTHONPATH="$PYBIGWIG_DIR${PYTHONPATH:+:$PYTHONPATH}"

"$PYTHON" -c 'import pyBigWig; print("PYBIGWIG_VERSION", pyBigWig.__version__)'

STATS="$META/Hmel215003o_scaffold_means.tsv"
printf 'sample\tscaffold_length\tscaffold_mean_depth\tbedgraph\n' > "$STATS.tmp"

for sample in "${samples[@]}"
do
    bw="$BW_DIR/${sample}.dedup.bw"

    if [[ ! -s "$bw" ]]
    then
        echo "DOWNLOADING_BIGWIG $sample"

        /usr/bin/curl \
            --fail \
            --location \
            --retry 5 \
            --output "$bw.part" \
            "$UPSTREAM_BASE/${sample}.dedup.bw"

        mv "$bw.part" "$bw"
    fi

    bg="$BG_DIR/${sample}.Hmel215003o.bedGraph"
    one_stats="$META/${sample}.Hmel215003o.stats.tsv"

    "$PYTHON" \
        "$CODE_DIR/extract_Hmel2_5_bigwig_scaffold.py" \
        --bigwig "$bw" \
        --chrom "$CHROM" \
        --bedgraph "$bg.tmp" \
        --stats "$one_stats.tmp" \
        --sample "$sample"

    mv "$bg.tmp" "$bg"

    # The helper records the temporary bedGraph path while writing; replace it
    # with the final promoted path in the per-sample stats row.
    awk \
        -F '\t' \
        -v OFS='\t' \
        -v final_bg="$bg" \
        'NR==1 {print; next} {$4=final_bg; print}' \
        "$one_stats.tmp" \
        > "$one_stats"

    rm -f "$one_stats.tmp"

    tail -n 1 "$one_stats" >> "$STATS.tmp"

    echo -e "BIGWIG_DONE\t$sample\t$bw"
done

mv "$STATS.tmp" "$STATS"

"$PYTHON" \
    "$CODE_DIR/plot_Yb_Hmel2_5_Figure8_normalized.py" \
    --bedgraph-dir "$BG_DIR" \
    --scaffold-stats "$STATS" \
    --output-dir "$RESULTS"

echo "YB_HMEL2_5_NORMALIZED_FINISHED"
