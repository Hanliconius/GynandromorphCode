#!/usr/bin/env bash
#SBATCH --job-name=optixDennis
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/martinlab/joe/GynandroCalling/logs/96b_scan_broad_optix_Dennis.nano.%j.out
#SBATCH --error=/scratch/martinlab/joe/GynandroCalling/logs/96b_scan_broad_optix_Dennis.nano.%j.err

project_dir="/scratch/martinlab/joe/GynandroCalling"
script="${project_dir}/scripts/96a_scan_broad_optix_Dennis.py"

cd "${project_dir}" || {
    printf 'ERROR: could not enter:\n%s\n' "${project_dir}"
}

if [[ -x "${GYNANDRO_PYTHON}" && -s "${script}" ]]
then
    printf 'Host: %s\n' "$(hostname)"
    printf 'Started: %s\n' "$(date -Is)"
    printf 'Python: %s\n' "${GYNANDRO_PYTHON}"
    printf 'Script: %s\n\n' "${script}"

    "${GYNANDRO_PYTHON}" "${script}"

    command_status="$?"

    printf '\nFinished: %s\n' "$(date -Is)"
    printf 'Python command status: %s\n' "${command_status}"
else
    printf '%s\n' 'ERROR: Python executable or script is missing.'
    printf 'Python: %s\n' "${GYNANDRO_PYTHON}"
    printf 'Script: %s\n' "${script}"
fi
