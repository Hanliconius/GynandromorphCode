#!/bin/bash
#SBATCH --job-name=smudgeKrR
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --array=0-3
#SBATCH --output=logs/smudgeKrR_%A_%a.out
#SBATCH --error=logs/smudgeKrR_%A_%a.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
OUTROOT="${ROOT}/analysis/smudgeplot_arrangements"

export PATH="$HOME/.local/share/mamba/envs/smudgeplot/bin:${PATH}"

SAMPLES=(
    ERR1050907
    ERR1050958
    ERR1050994
    ERR1051012
)

LS=(
    0
    2
    3
    4
)

s="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
L="${LS[$SLURM_ARRAY_TASK_ID]}"
dir="${OUTROOT}/${s}_krista"

cd "$dir"

tmp="${dir}/fastk_tmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
mkdir -p "$tmp"

echo "SAMPLE = $s"
echo "DIR    = $dir"
echo "TMP    = $tmp"

# ERR1050907 needs FastK rerun.
if [[ ! -s FastK_Table.ktab || ! -s FastK_Table.hist ]]; then
    rm -f FastK_Table.ktab FastK_Table.hist

    FastK \
        -v \
        -t4 \
        -k22 \
        -M100 \
        -T60 \
        -P"$tmp" \
        ./*.fastq.gz \
        -NFastK_Table
fi

test -s FastK_Table.ktab
test -s FastK_Table.hist

Histex -A -h1:100 FastK_Table.hist > FastK_Table.hist.tsv

# For task 0, stop after FastK so we can choose L from its histogram.
if [[ "$SLURM_ARRAY_TASK_ID" -eq 0 ]]; then
    echo
    echo "ERR1050907 histogram ready; stopping before hetmers."
    awk '$1 >= 1 && $1 <= 30 {print}' FastK_Table.hist.tsv
    rm -rf "$tmp"
    exit 0
fi

echo
echo "Using L=${L}"

rm -f "kmerpairs_L${L}.smu" output_L${L}*

smudgeplot hetmers \
    -L "$L" \
    -t 60 \
    -o "kmerpairs_L${L}" \
    --verbose \
    FastK_Table

test -s "kmerpairs_L${L}.smu"

smudgeplot all \
    --format pdf \
    -o "output_L${L}" \
    "kmerpairs_L${L}.smu"

echo
ls -lh "kmerpairs_L${L}.smu" output_L${L}* || true

rm -rf "$tmp"
