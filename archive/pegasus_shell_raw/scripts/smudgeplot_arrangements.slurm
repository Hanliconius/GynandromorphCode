#!/bin/bash
#SBATCH --job-name=smudgeArr
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --array=0-5%3
#SBATCH --output=logs/smudgeArr_%A_%a.out
#SBATCH --error=logs/smudgeArr_%A_%a.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
cd "${ROOT}"

mkdir -p logs

MANIFEST="${ROOT}/analysis/smudgeplot_arrangements/smudgeplot_array_manifest.tsv"

line="$(
    awk -F $'\t' \
        -v id="${SLURM_ARRAY_TASK_ID}" \
        'NR > 1 && $1 == id {print; exit}' \
        "${MANIFEST}"
)"

[[ -n "${line}" ]] || {
    echo "No manifest entry for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
}

IFS=$'\t' read -r array_id label dir <<< "${line}"

echo "ARRAY_ID = ${array_id}"
echo "LABEL    = ${label}"
echo "DIR      = ${dir}"
echo
echo "INPUTS:"
ls -lh "${dir}"/*.fastq.gz
echo


# ------------------------------------------------------------
# Idempotent completed-job guard
# ------------------------------------------------------------

if [[ -s "${dir}/output.smudge_report.tsv" ]] && \
   [[ -s "${dir}/kmerpairs.smu" ]]
then
    echo "Already complete: ${label}"
    exit 0
fi


# ------------------------------------------------------------
# Run exact old parameters inside known environment
# ------------------------------------------------------------

micromamba run -n smudgeplot \
    FastK \
        -v \
        -t4 \
        -k22 \
        -M100 \
        -T60 \
        "${dir}"/*.fastq.gz \
        -N"${dir}/FastK_Table"


micromamba run -n smudgeplot \
    smudgeplot hetmers \
        -L 12 \
        -t 60 \
        -o "${dir}/kmerpairs" \
        --verbose \
        "${dir}/FastK_Table"


micromamba run -n smudgeplot \
    smudgeplot all \
        -o "${dir}/output" \
        "${dir}/kmerpairs.smu"


echo
echo "FINISHED ${label}"
echo

ls -lh \
    "${dir}/kmerpairs.smu" \
    "${dir}"/output*
