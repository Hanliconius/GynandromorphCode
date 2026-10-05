#!/bin/bash
#SBATCH --job-name=smudgeDown
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-4
#SBATCH --output=logs/smudgeDown_%A_%a.out
#SBATCH --error=logs/smudgeDown_%A_%a.err

set -euo pipefail

SMUDGE_BIN="$HOME/.local/share/mamba/envs/smudgeplot/bin"
export PATH="${SMUDGE_BIN}:${PATH}"

echo "smudgeplot = $(command -v smudgeplot)"
echo "FastK      = $(command -v FastK)"

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
BASE="${ROOT}/analysis/smudgeplot_arrangements"

DIRS=(
    "${BASE}/TN_Left"
    "${BASE}/TN_Right"
    "${BASE}/TN_T2_merge"
    "${BASE}/ERR2298208_single"
    "${BASE}/pseudoTetraploid_ERR2298208_ERR1307111"
)

dir="${DIRS[$SLURM_ARRAY_TASK_ID]}"
label="$(basename "${dir}")"

echo "TASK  = ${SLURM_ARRAY_TASK_ID}"
echo "LABEL = ${label}"
echo "DIR   = ${dir}"
echo

[[ -s "${dir}/FastK_Table.ktab" ]] || {
    echo "ERROR: missing ${dir}/FastK_Table.ktab"
    exit 1
}

[[ -s "${dir}/FastK_Table.hist" ]] || {
    echo "ERROR: missing ${dir}/FastK_Table.hist"
    exit 1
}

# Critical fix:
# each task works INSIDE its own directory so temporary .trim
# products cannot collide between array jobs.
cd "${dir}"

echo "Running hetmers..."

smudgeplot hetmers \
        -L 12 \
        -t 60 \
        -o kmerpairs \
        --verbose \
        FastK_Table

[[ -s kmerpairs.smu ]] || {
    echo "ERROR: hetmers did not produce kmerpairs.smu"
    exit 1
}

echo
echo "Running smudgeplot all..."

smudgeplot all \
        -o output \
        kmerpairs.smu

[[ -s output.smudge_report.tsv ]] || {
    echo "ERROR: missing output.smudge_report.tsv"
    exit 1
}

echo
echo "SUCCESS: ${label}"
ls -lh \
    kmerpairs.smu \
    output.smudge_report.tsv \
    output_smudgeplot.png \
    output_smudgeplot_log10.png
