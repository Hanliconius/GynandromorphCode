#!/usr/bin/env bash
#SBATCH --job-name=gyn205sib
#SBATCH --partition=nano
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --chdir=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
#SBATCH --output=logs/205_validate_empirical_sibling_calibration.%j.out
#SBATCH --error=logs/205_validate_empirical_sibling_calibration.%j.err

set -euo pipefail

result_root=analysis/genomewide_mixing_corrected_v1/empirical_full_sib_calibration_v1

required=(
    "${result_root}/final/all_empirical_sibling_simulations.tsv.gz"
    "${result_root}/final/empirical_sibling_simulation_summary.tsv"
    "${result_root}/final/representative_simulations.tsv"
    "${result_root}/final/observed_and_simulated_1Mb_evidence.tsv"
    "${result_root}/final/collection.audit.json"
    "${result_root}/final/window_evidence.audit.json"
    "${result_root}/figures/empirical_full_sib_mixing_narrative.pdf"
    "${result_root}/figures/empirical_full_sib_AF_summary.png"
    "${result_root}/figures/empirical_full_sib_narrative.audit.json"
)

file_status=0

for path in "${required[@]}"
do
    if [[ -s "${path}" ]]
    then
        printf 'PRESENT\t%s\t%s bytes\n' \
            "${path}" \
            "$(wc -c < "${path}")"
    else
        printf 'MISSING_OR_EMPTY\t%s\n' "${path}"
        file_status=1
    fi
done

gzip_status=0

if ! gzip -t \
    "${result_root}/final/all_empirical_sibling_simulations.tsv.gz"
then
    gzip_status=1
fi

row_count=$(
    gzip -dc \
        "${result_root}/final/all_empirical_sibling_simulations.tsv.gz" \
        | awk 'NR > 1 {count++} END {print count + 0}'
)

summary_rows=$(
    awk 'NR > 1 {count++} END {print count + 0}' \
        "${result_root}/final/empirical_sibling_simulation_summary.tsv"
)

representative_rows=$(
    awk 'NR > 1 {count++} END {print count + 0}' \
        "${result_root}/final/representative_simulations.tsv"
)

count_status=0

if [[ "${row_count}" -ne 1320 ]]
then
    count_status=1
fi

if [[ "${summary_rows}" -ne 4 ]]
then
    count_status=1
fi

if [[ "${representative_rows}" -ne 5 ]]
then
    count_status=1
fi

printf '\nVALIDATION COUNTS\n'
printf 'simulation_rows=%s expected=1320\n' "${row_count}"
printf 'summary_rows=%s expected=4\n' "${summary_rows}"
printf 'representative_rows=%s expected=5\n' "${representative_rows}"

pipeline_status=FAIL

if [[ \
    "${file_status}" -eq 0 && \
    "${gzip_status}" -eq 0 && \
    "${count_status}" -eq 0 \
]]
then
    pipeline_status=PASS
fi

audit="${result_root}/final/pipeline_validation.tsv"

{
    printf 'item\tvalue\n'
    printf 'required_file_status\t%s\n' "${file_status}"
    printf 'gzip_status\t%s\n' "${gzip_status}"
    printf 'simulation_rows\t%s\n' "${row_count}"
    printf 'summary_rows\t%s\n' "${summary_rows}"
    printf 'representative_rows\t%s\n' "${representative_rows}"
    printf 'pipeline_status\t%s\n' "${pipeline_status}"
} > "${audit}"

printf '\nPIPELINE STATUS\t%s\n' "${pipeline_status}"

if [[ "${pipeline_status}" == PASS ]]
then
    printf 'ENDPOINT=EMPIRICAL_FULL_SIBLING_CALIBRATION_READY\n'

    printf 'NARRATIVE_PDF\t%s\n' \
        "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/${result_root}/figures/empirical_full_sib_mixing_narrative.pdf"
else
    false
fi
