#!/usr/bin/env bash
#SBATCH --job-name=optixRD
#SBATCH --partition=nano
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/martinlab/joe/GynandroCalling/logs/96a_scan_broad_optix_Ray_Dennis.nano.%A_%a.out
#SBATCH --error=/scratch/martinlab/joe/GynandroCalling/logs/96a_scan_broad_optix_Ray_Dennis.nano.%A_%a.err

project_dir="/scratch/martinlab/joe/GynandroCalling"

scripts=(
    "${project_dir}/scripts/96a_scan_broad_optix_Ray.py"
    "${project_dir}/scripts/96a_scan_broad_optix_Dennis.py"
)

labels=(
    "Ray"
    "Dennis"
)

task_index="${SLURM_ARRAY_TASK_ID}"
script="${scripts[${task_index}]}"
label="${labels[${task_index}]}"

cd "${project_dir}" || {
    printf 'ERROR: could not enter:\n%s\n' "${project_dir}"
}

if [[ -s "${script}" ]]
then
    printf 'Trait: %s\n' "${label}"
    printf 'Host: %s\n' "$(hostname)"
    printf 'Started: %s\n' "$(date -Is)"
    printf 'Script: %s\n\n' "${script}"

    micromamba run \
        -n gynandrocalling \
        python "${script}"

    command_status="$?"

    printf '\nFinished: %s\n' "$(date -Is)"
    printf 'Python command status: %s\n' "${command_status}"
else
    printf 'ERROR: script is missing:\n%s\n' "${script}"
fi
