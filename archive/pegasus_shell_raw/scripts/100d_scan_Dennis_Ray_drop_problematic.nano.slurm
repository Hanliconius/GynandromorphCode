#!/usr/bin/env bash
#SBATCH --job-name=gynDRdrop6
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=0-1
#SBATCH --output=/scratch/martinlab/joe/GynandroCalling/logs/100d_scan_Dennis_Ray_drop_problematic.nano.%A_%a.out
#SBATCH --error=/scratch/martinlab/joe/GynandroCalling/logs/100d_scan_Dennis_Ray_drop_problematic.nano.%A_%a.err

project_dir="/scratch/martinlab/joe/GynandroCalling"
env_python="/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python"

cd "${project_dir}"

if [[ "${SLURM_ARRAY_TASK_ID}" == "0" ]]
then
    analysis_script="scripts/100d_scan_Dennis_drop_problematic.py"
else
    analysis_script="scripts/100d_scan_Ray_drop_problematic.py"
fi

if [[ -x "${env_python}" && -s "${analysis_script}" ]]
then
    "${env_python}" "${analysis_script}"
else
    printf 'ERROR: required Python or analysis script is missing.\n' >&2
    printf 'Python: %s\n' "${env_python}" >&2
    printf 'Script: %s\n' "${analysis_script}" >&2
    false
fi
