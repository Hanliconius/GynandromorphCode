#!/usr/bin/env bash
#SBATCH --job-name=gyn169bs
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=logs/169e_build_24sample_calling_set.%j.out
#SBATCH --error=logs/169e_build_24sample_calling_set.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin

samtools="${tool_bin}/samtools"

download_manifest="${project}/analysis/genomewide_mixing_input/metadata/resolved_download_manifest.tsv"
bam_dir="${project}/analysis/genomewide_mixing_input/alignment/libraries"
set_dir="${project}/analysis/genomewide_mixing_input/metadata/bam_sets_24"

mkdir -p "${set_dir}"

controls_tsv="${set_dir}/controls.tsv"
all_tsv="${set_dir}/all_24_samples.tsv"
groups_tsv="${set_dir}/all_24_samples.groups.tsv"

printf 'sample\trole\tbam\n' \
    > "${controls_tsv}.new"

awk \
    -F '\t' \
    -v OFS='\t' \
    -v directory="${bam_dir}" \
    '
        NR > 1 && $5 == "control" {
            print $4, "control",
                directory "/" $4 ".dedup.bam"
        }
    ' \
    "${download_manifest}" \
    >> "${controls_tsv}.new"

cp \
    "${controls_tsv}.new" \
    "${all_tsv}.new"

for sample in \
    TN_Left_T2 \
    TN_Left_T3 \
    TN_Right_T2 \
    TN_Right_T3 \
    TN_Left_merged \
    TN_Right_merged
do
    case "${sample}" in
        TN_Left*)
            role=focal_left
            ;;
        TN_Right*)
            role=focal_right
            ;;
    esac

    printf '%s\t%s\t%s/%s.dedup.bam\n' \
        "${sample}" \
        "${role}" \
        "${bam_dir}" \
        "${sample}" \
        >> "${all_tsv}.new"
done

missing=0

while IFS=$'\t' read -r sample role bam
do
    [[ "${sample}" == "sample" ]] && continue

    if [[ ! -s "${bam}" || ! -s "${bam}.bai" ]] ||
        ! "${samtools}" quickcheck -q "${bam}"
    then
        printf 'MISSING_OR_INVALID\t%s\t%s\n' \
            "${sample}" \
            "${bam}" >&2
        missing=1
    fi
done < "${all_tsv}.new"

controls_n=$(
    awk 'END {print NR - 1}' \
        "${controls_tsv}.new"
)

all_n=$(
    awk 'END {print NR - 1}' \
        "${all_tsv}.new"
)

if [[
    "${missing}" -ne 0 ||
    "${controls_n}" -ne 18 ||
    "${all_n}" -ne 24
]]
then
    printf 'ERROR: calling-set validation failed: controls=%s all=%s missing=%s\n' \
        "${controls_n}" \
        "${all_n}" \
        "${missing}" >&2
    false
fi

mv \
    "${controls_tsv}.new" \
    "${controls_tsv}"

mv \
    "${all_tsv}.new" \
    "${all_tsv}"

awk \
    -F '\t' \
    'NR > 1 {print $3}' \
    "${controls_tsv}" \
    > "${set_dir}/controls.bamlist.new"

mv \
    "${set_dir}/controls.bamlist.new" \
    "${set_dir}/controls.bamlist"

awk \
    -F '\t' \
    'NR > 1 {print $3}' \
    "${all_tsv}" \
    > "${set_dir}/all_24_samples.bamlist.new"

mv \
    "${set_dir}/all_24_samples.bamlist.new" \
    "${set_dir}/all_24_samples.bamlist"

: > "${groups_tsv}.new"

awk \
    -F '\t' \
    '
        NR > 1 {
            if ($2 == "control") {
                print $1 "\tcontrols"
            } else {
                print $1 "\t" $1
            }
        }
    ' \
    "${all_tsv}" \
    >> "${groups_tsv}.new"

mv \
    "${groups_tsv}.new" \
    "${groups_tsv}"

printf 'CALLING SET COMPLETE\n'
printf 'controls=%s\n' "${controls_n}"
printf 'focal_representations=6\n'
printf 'total_samples=%s\n' "${all_n}"
printf '\nSAMPLES\n'
cat "${all_tsv}"
