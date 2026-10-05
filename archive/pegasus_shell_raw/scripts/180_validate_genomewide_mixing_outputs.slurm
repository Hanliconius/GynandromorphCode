#!/usr/bin/env bash
#SBATCH --job-name=gyn180va
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/180_validate_genomewide_mixing_outputs.%j.out
#SBATCH --error=logs/180_validate_genomewide_mixing_outputs.%j.err

set -euo pipefail

stage_name=validate_genomewide_mixing_outputs
start_seconds=${SECONDS}
printf 'STAGE_START\tstage=%s\tjob=%s\thost=%s\ttime=%s\n' \
    "${stage_name}" "${SLURM_JOB_ID:-NA}" "$(hostname)" "$(date --iso-8601=seconds)"
trap 'status=$?; printf "STAGE_END\tstage=%s\tstatus=%s\telapsed_seconds=%s\thost=%s\ttime=%s\n" "${stage_name}" "${status}" "$((SECONDS - start_seconds))" "$(hostname)" "$(date --iso-8601=seconds)"' EXIT

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
root="${project}/analysis/genomewide_mixing"
audit="${root}/final/genomewide_pipeline.audit.tsv"

mkdir -p "${root}/final"

extracted=$(find "${root}/inputs" -maxdepth 1 -type f -name '*.counts.tsv.gz' | wc -l)
diagnostic_calls=$(find "${root}/observed_fit" -maxdepth 1 -type f -name '*.adjusted_calls.tsv.gz' | wc -l)
final_calls=$(find "${root}/final/adjusted_calls" -maxdepth 1 -type f -name '*.mixing_adjusted_true_states.tsv.gz' | wc -l)
simulations=$(find "${root}/simulations" -maxdepth 1 -type f -name 'simulation_task_*.tsv' | wc -l)
expected_blocks=$(awk 'NR > 1 {count++} END {print count + 0}' "${root}/observed_fit/autosomal_5Mb_block_manifest.tsv")
jackknives=$(find "${root}/block_jackknife" -maxdepth 1 -type f -name 'block_*.tsv' ! -name '*.skip.tsv' | wc -l)

required_files=(
    "${root}/observed_fit/observed_founder_model_fits.tsv"
    "${root}/observed_fit/genomewide_2Mb_window_evidence.tsv"
    "${root}/observed_fit/model_metadata.json"
    "${root}/final/calibrated_founder_and_mixing_summary.tsv"
    "${root}/final/mixing_calibration_curves.tsv"
    "${root}/final/all_simulation_results.tsv.gz"
    "${root}/final/all_5Mb_block_jackknife_results.tsv"
    "${root}/final/jackknife_block_evidence_manhattan.tsv"
    "${root}/final/final_model_selection.json"
    "${root}/figures/genomewide_three_vs_four_founder_mixing_narrative.pdf"
    "${root}/figures/genomewide_founder_mixing_summary.png"
    "${root}/figures/genomewide_narrative.audit.tsv"
)

required_status=0

for path in "${required_files[@]}"
do
    if [[ -s "${path}" ]]
    then
        printf 'PRESENT\t%s\t%s bytes\n' "${path}" "$(wc -c < "${path}")"
    else
        printf 'MISSING_OR_EMPTY\t%s\n' "${path}" >&2
        required_status=1
    fi
done

gzip_status=0
while IFS= read -r path
do
    if ! gzip -t "${path}"
    then
        printf 'INVALID_GZIP\t%s\n' "${path}" >&2
        gzip_status=1
    fi
done < <(
    find "${root}/inputs" "${root}/final/adjusted_calls" \
        -maxdepth 1 -type f -name '*.tsv.gz' -print
)

count_status=0
if [[
    "${extracted}" -ne 21 ||
    "${diagnostic_calls}" -ne 21 ||
    "${final_calls}" -ne 21 ||
    "${simulations}" -ne 98 ||
    "${expected_blocks}" -lt 1 ||
    "${jackknives}" -ne "${expected_blocks}"
]]
then
    count_status=1
fi

printf 'item\tvalue\n' > "${audit}.new"
printf 'extracted_chromosomes\t%s\n' "${extracted}" >> "${audit}.new"
printf 'diagnostic_call_chromosomes\t%s\n' "${diagnostic_calls}" >> "${audit}.new"
printf 'final_adjusted_call_chromosomes\t%s\n' "${final_calls}" >> "${audit}.new"
printf 'simulation_tasks\t%s\n' "${simulations}" >> "${audit}.new"
printf 'expected_jackknife_blocks\t%s\n' "${expected_blocks}" >> "${audit}.new"
printf 'completed_jackknife_blocks\t%s\n' "${jackknives}" >> "${audit}.new"
printf 'required_file_status\t%s\n' "${required_status}" >> "${audit}.new"
printf 'gzip_status\t%s\n' "${gzip_status}" >> "${audit}.new"
printf 'count_status\t%s\n' "${count_status}" >> "${audit}.new"

if [[ "${required_status}" -eq 0 && "${gzip_status}" -eq 0 && "${count_status}" -eq 0 ]]
then
    printf 'pipeline_status\tPASS\n' >> "${audit}.new"
    mv "${audit}.new" "${audit}"
    printf '\nGENOME-WIDE MIXING PIPELINE VALIDATION PASSED\n'
    cat "${audit}"
else
    printf 'pipeline_status\tFAIL\n' >> "${audit}.new"
    mv "${audit}.new" "${audit}"
    printf '\nGENOME-WIDE MIXING PIPELINE VALIDATION FAILED\n' >&2
    cat "${audit}" >&2
    false
fi
