#!/usr/bin/env bash
#SBATCH --job-name=gyn_band502
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/martinlab/joe/GynandroCalling/logs/95d_scan_long_band_original502_WGS_vs_WGSplusCS.%j.out
#SBATCH --error=/scratch/martinlab/joe/GynandroCalling/logs/95d_scan_long_band_original502_WGS_vs_WGSplusCS.%j.err

project_dir="/scratch/martinlab/joe/GynandroCalling"
script="${project_dir}/scripts/95d_scan_long_band_original502_WGS_vs_WGSplusCS.py"
vcf="${project_dir}/external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf"

cd "${project_dir}" || {
    printf 'ERROR: could not enter project directory:\n%s\n' \
        "${project_dir}"
}

if [[ -s "${script}" && -s "${vcf}" ]]
then
    printf 'Host: %s\n' "$(hostname)"
    printf 'Started: %s\n' "$(date -Is)"
    printf 'Input VCF: %s\n' "${vcf}"
    printf 'Input bytes: %s\n\n' "$(stat -c '%s' "${vcf}")"

    micromamba run \
        -n gynandrocalling \
        python "${script}"

    command_status="$?"

    printf '\nFinished: %s\n' "$(date -Is)"
    printf 'Python command status: %s\n' "${command_status}"
else
    printf '%s\n' 'ERROR: required script or VCF is missing.'
    printf 'Script: %s\n' "${script}"
    printf 'VCF: %s\n' "${vcf}"
fi
