#!/usr/bin/env bash

set -u

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
python_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
calling_array_job="${1:-}"

if cd "${project}"
then
    :
else
    printf 'FAIL\tcannot enter project directory: %s\n' "${project}"
    exit 1
fi

mkdir -p logs

scripts=(
    scripts/170e_recover_joint_chromosome_vcfs.slurm
    scripts/171e_concat_recovered_joint_vcfs.slurm
    scripts/173_extract_joint_counts.slurm
    scripts/174_fit_observed_genomewide_mixing.slurm
    scripts/175_simulate_genomewide_mixing.slurm
    scripts/176_jackknife_genomewide_mixing.slurm
    scripts/177_collect_genomewide_mixing.slurm
    scripts/178_finalize_mixing_adjusted_calls.slurm
    scripts/179_plot_genomewide_mixing.slurm
    scripts/180_validate_genomewide_mixing_outputs.slurm
)

preflight_status=0
printf '%s\n' 'GENOME-WIDE MIXING PIPELINE PREFLIGHT'

if [[ "${calling_array_job}" =~ ^[0-9]+$ ]]
then
    printf 'PASS\texisting joint-calling array job: %s\n' "${calling_array_job}"
else
    printf 'FAIL\tprovide the existing joint-calling array job ID; received: %s\n' "${calling_array_job}"
    preflight_status=1
fi

for script in "${scripts[@]}"
do
    if [[ -s "${script}" ]] && bash -n "${script}"
    then
        printf 'PASS\t%s\n' "${script}"
    else
        printf 'FAIL\t%s\n' "${script}"
        preflight_status=1
    fi
done

for tool in bcftools bgzip tabix
do
    if [[ -x "${tool_bin}/${tool}" ]]
    then
        printf 'PASS\ttool: %s\n' "${tool_bin}/${tool}"
    else
        printf 'FAIL\ttool: %s\n' "${tool_bin}/${tool}"
        preflight_status=1
    fi
done

original_caller=scripts/170d_call_joint_24sample_chromosome_vcfs.slurm
if [[ -s "${original_caller}" ]] \
    && bash -n "${original_caller}" \
    && grep -Fq '${joint_bcf}.new' "${original_caller}"
then
    printf 'PASS\toriginal caller can be patched safely at recovery time\n'
else
    printf 'FAIL\toriginal caller is missing, invalid, or lacks the expected temporary-output token\n'
    preflight_status=1
fi

if [[ -s scripts/173_genomewide_mixing_models.py ]] \
    && "${python_bin}" -m py_compile scripts/173_genomewide_mixing_models.py \
    && "${python_bin}" -c 'import numpy, pandas, scipy, matplotlib'
then
    printf 'PASS\tPython model program and dependencies\n'
else
    printf 'FAIL\tPython model program or dependencies\n'
    preflight_status=1
fi

manifest=analysis/genomewide_mixing_input/metadata/chromosome_calling_manifest.tsv
if [[ -s "${manifest}" ]]
then
    manifest_rows=$(awk 'NR > 1 {count++} END {print count + 0}' "${manifest}")
else
    manifest_rows=0
fi

if [[ "${manifest_rows}" -eq 21 ]]
then
    printf 'PASS\tchromosome manifest rows: %s\n' "${manifest_rows}"
else
    printf 'FAIL\tchromosome manifest rows: %s\n' "${manifest_rows}"
    preflight_status=1
fi

printf 'preflight_status=%s\n' "${preflight_status}"

recovery_status=127
concat_status=127
extract_status=127
fit_status=127
simulation_status=127
jackknife_status=127
collect_status=127
adjust_status=127
plot_status=127
validation_status=127

recovery_job=""
concat_job=""
extract_job=""
fit_job=""
simulation_job=""
jackknife_job=""
collect_job=""
adjust_job=""
plot_job=""
validation_job=""

if [[ "${preflight_status}" -eq 0 ]]
then
    recovery_submission=$(sbatch --parsable --dependency="afterany:${calling_array_job}" --array=0-20%21 scripts/170e_recover_joint_chromosome_vcfs.slurm)
    recovery_status=$?
fi

if [[ "${recovery_status}" -eq 0 ]]
then
    recovery_job=${recovery_submission%%;*}
    concat_submission=$(sbatch --parsable --dependency="afterok:${recovery_job}" scripts/171e_concat_recovered_joint_vcfs.slurm)
    concat_status=$?
fi

if [[ "${concat_status}" -eq 0 ]]
then
    concat_job=${concat_submission%%;*}
    extract_submission=$(sbatch --parsable --dependency="afterok:${concat_job}" --array=0-20 scripts/173_extract_joint_counts.slurm)
    extract_status=$?
fi

if [[ "${extract_status}" -eq 0 ]]
then
    extract_job=${extract_submission%%;*}
    fit_submission=$(sbatch --parsable --dependency="afterok:${extract_job}" scripts/174_fit_observed_genomewide_mixing.slurm)
    fit_status=$?
fi

if [[ "${fit_status}" -eq 0 ]]
then
    fit_job=${fit_submission%%;*}
    simulation_submission=$(sbatch --parsable --dependency="afterok:${fit_job}" --array=0-97%40 scripts/175_simulate_genomewide_mixing.slurm)
    simulation_status=$?
    jackknife_submission=$(sbatch --parsable --dependency="afterok:${fit_job}" --array=0-99%30 scripts/176_jackknife_genomewide_mixing.slurm)
    jackknife_status=$?
fi

if [[ "${simulation_status}" -eq 0 && "${jackknife_status}" -eq 0 ]]
then
    simulation_job=${simulation_submission%%;*}
    jackknife_job=${jackknife_submission%%;*}
    collect_submission=$(sbatch --parsable --dependency="afterok:${simulation_job}:${jackknife_job}" scripts/177_collect_genomewide_mixing.slurm)
    collect_status=$?
fi

if [[ "${collect_status}" -eq 0 ]]
then
    collect_job=${collect_submission%%;*}
    adjust_submission=$(sbatch --parsable --dependency="afterok:${collect_job}" --array=0-20 scripts/178_finalize_mixing_adjusted_calls.slurm)
    adjust_status=$?
fi

if [[ "${adjust_status}" -eq 0 ]]
then
    adjust_job=${adjust_submission%%;*}
    plot_submission=$(sbatch --parsable --dependency="afterok:${adjust_job}" scripts/179_plot_genomewide_mixing.slurm)
    plot_status=$?
fi

if [[ "${plot_status}" -eq 0 ]]
then
    plot_job=${plot_submission%%;*}
    validation_submission=$(sbatch --parsable --dependency="afterok:${plot_job}" scripts/180_validate_genomewide_mixing_outputs.slurm)
    validation_status=$?
fi

if [[ "${validation_status}" -eq 0 ]]
then
    validation_job=${validation_submission%%;*}
fi

printf '\nSUBMISSION STATUS\n'
printf 'recovery_status=%s\n' "${recovery_status}"
printf 'concat_status=%s\n' "${concat_status}"
printf 'extract_status=%s\n' "${extract_status}"
printf 'fit_status=%s\n' "${fit_status}"
printf 'simulation_status=%s\n' "${simulation_status}"
printf 'jackknife_status=%s\n' "${jackknife_status}"
printf 'collect_status=%s\n' "${collect_status}"
printf 'adjust_status=%s\n' "${adjust_status}"
printf 'plot_status=%s\n' "${plot_status}"
printf 'validation_status=%s\n' "${validation_status}"

if [[ "${validation_status}" -eq 0 ]]
then
    printf '\nSUBMITTED JOB CHAIN\n'
    printf 'existing_joint_call_array\t%s\n' "${calling_array_job}"
    printf 'chromosome_recovery_array\t%s\tarray=0-20%%21\n' "${recovery_job}"
    printf 'corrected_genomewide_concat\t%s\n' "${concat_job}"
    printf 'count_extraction_array\t%s\tarray=0-20\n' "${extract_job}"
    printf 'observed_fit\t%s\n' "${fit_job}"
    printf 'simulation_array\t%s\tarray=0-97%%40\n' "${simulation_job}"
    printf 'jackknife_array\t%s\tarray=0-99%%30\n' "${jackknife_job}"
    printf 'calibration_collection\t%s\n' "${collect_job}"
    printf 'final_adjusted_call_array\t%s\tarray=0-20\n' "${adjust_job}"
    printf 'narrative_plotting\t%s\n' "${plot_job}"
    printf 'final_validation\t%s\n' "${validation_job}"

    printf '\nCURRENT QUEUE\n'
    squeue -j "${recovery_job},${concat_job},${extract_job},${fit_job},${simulation_job},${jackknife_job},${collect_job},${adjust_job},${plot_job},${validation_job}"
else
    printf '\nThe complete downstream chain was not submitted.\n'
    false
fi
