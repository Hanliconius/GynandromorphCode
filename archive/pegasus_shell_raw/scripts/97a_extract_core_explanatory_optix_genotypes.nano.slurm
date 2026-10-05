#!/usr/bin/env bash
#SBATCH --job-name=optixCoreGT
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/martinlab/joe/GynandroCalling/logs/97a_extract_core_explanatory_optix_genotypes.nano.%j.out
#SBATCH --error=/scratch/martinlab/joe/GynandroCalling/logs/97a_extract_core_explanatory_optix_genotypes.nano.%j.err

project_dir="/scratch/martinlab/joe/GynandroCalling"
script="${project_dir}/scripts/97a_extract_core_explanatory_optix_genotypes.py"

cd "${project_dir}" || {
    printf 'ERROR: could not enter:\n%s\n' "${project_dir}"
}

if [[ -x "${GYNANDRO_PYTHON}" && -s "${script}" ]]
then
    printf 'Host: %s\n' "$(hostname)"
    printf 'Started: %s\n' "$(date -Is)"

    "${GYNANDRO_PYTHON}" "${script}"

    status="$?"

    printf '\nFinished: %s\n' "$(date -Is)"
    printf 'Python command status: %s\n' "${status}"
else
    printf '%s\n' 'ERROR: Python executable or script is missing.'
fi
