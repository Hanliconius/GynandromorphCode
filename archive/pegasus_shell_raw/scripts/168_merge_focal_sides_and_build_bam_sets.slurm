#!/usr/bin/env bash
#SBATCH --job-name=gyn168mg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=3
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/168_merge_focal_sides_and_build_bam_sets.%j.out
#SBATCH --error=logs/168_merge_focal_sides_and_build_bam_sets.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
samtools="${tool_bin}/samtools"

manifest="${project}/analysis/genomewide_mixing_input/metadata/resolved_download_manifest.tsv"
library_dir="${project}/analysis/genomewide_mixing_input/alignment/libraries"
focal_dir="${project}/analysis/genomewide_mixing_input/alignment/focal_sets"
set_dir="${project}/analysis/genomewide_mixing_input/metadata/bam_sets"
tmp_dir="${project}/analysis/genomewide_mixing_input/tmp/focal_merge.${SLURM_JOB_ID}"

mkdir -p \
    "${project}/logs" \
    "${focal_dir}" \
    "${set_dir}" \
    "${tmp_dir}"

if [[ ! -s "${manifest}" || ! -x "${samtools}" ]]
then
    printf 'ERROR: manifest or samtools is unavailable.\n' >&2
    false
fi

missing=0

while IFS=$'\t' read -r \
    array_id \
    source_accession \
    resolved_run \
    library_id \
    role \
    side \
    leg \
    merge_group \
    layout \
    fastq1 \
    fastq2
do
    [[ "${array_id}" == "array_id" ]] && continue

    bam="${library_dir}/${library_id}.dedup.bam"

    if [[ ! -s "${bam}" || ! -s "${bam}.bai" ]] ||
        ! "${samtools}" quickcheck -q "${bam}"
    then
        printf 'MISSING_OR_INVALID\t%s\n' \
            "${bam}" >&2
        missing=1
    fi
done < "${manifest}"

if [[ "${missing}" -ne 0 ]]
then
    printf 'ERROR: not all 22 library BAMs are complete.\n' >&2
    false
fi

retag_one()
{
    local input_bam=$1
    local output_bam=$2
    local rg_id=$3
    local sample_id=$4

    "${samtools}" addreplacerg \
        -@ 2 \
        -m overwrite_all \
        -r "ID:${rg_id}" \
        -r "SM:${sample_id}" \
        -r "LB:${rg_id}" \
        -r "PL:ILLUMINA" \
        -o "${output_bam}" \
        "${input_bam}"
}

build_side()
{
    local side_name=$1
    local t2_library=$2
    local t3_library=$3
    local final_bam="${focal_dir}/${side_name}.dedup_merged.bam"

    if [[
        -s "${final_bam}" &&
        -s "${final_bam}.bai"
    ]] && "${samtools}" quickcheck -q "${final_bam}"
    then
        printf 'Already complete and valid: %s\n' \
            "${final_bam}"
    else
        local t2_tmp="${tmp_dir}/${side_name}.T2.retagged.bam"
        local t3_tmp="${tmp_dir}/${side_name}.T3.retagged.bam"
        local merged_tmp="${tmp_dir}/${side_name}.merged.bam"

        retag_one \
            "${library_dir}/${t2_library}.dedup.bam" \
            "${t2_tmp}" \
            "${side_name}_T2" \
            "${side_name}"

        retag_one \
            "${library_dir}/${t3_library}.dedup.bam" \
            "${t3_tmp}" \
            "${side_name}_T3" \
            "${side_name}"

        "${samtools}" merge \
            -@ 2 \
            -f \
            "${merged_tmp}" \
            "${t2_tmp}" \
            "${t3_tmp}"

        "${samtools}" index \
            -@ 2 \
            "${merged_tmp}"

        "${samtools}" quickcheck \
            -v \
            "${merged_tmp}"

        mv "${merged_tmp}" "${final_bam}"
        mv "${merged_tmp}.bai" "${final_bam}.bai"

        rm -f \
            "${t2_tmp}" \
            "${t3_tmp}"
    fi
}

build_side \
    TN_Left_merged \
    TN_Left_T2 \
    TN_Left_T3

build_side \
    TN_Right_merged \
    TN_Right_T2 \
    TN_Right_T3

controls_tsv="${set_dir}/controls.tsv"
primary_tsv="${set_dir}/merged_T2_T3.tsv"
t2_tsv="${set_dir}/T2_only.tsv"
t3_tsv="${set_dir}/T3_only.tsv"

printf 'sample\trole\tbam\n' \
    > "${controls_tsv}.new"

awk \
    -F '\t' \
    -v OFS='\t' \
    -v directory="${library_dir}" \
    '
        NR > 1 && $5 == "control" {
            print $4, "control", directory "/" $4 ".dedup.bam"
        }
    ' \
    "${manifest}" \
    >> "${controls_tsv}.new"

cp "${controls_tsv}.new" "${primary_tsv}.new"

printf 'TN_Left_merged\tfocal_left\t%s/TN_Left_merged.dedup_merged.bam\n' \
    "${focal_dir}" \
    >> "${primary_tsv}.new"

printf 'TN_Right_merged\tfocal_right\t%s/TN_Right_merged.dedup_merged.bam\n' \
    "${focal_dir}" \
    >> "${primary_tsv}.new"

cp "${controls_tsv}.new" "${t2_tsv}.new"

printf 'TN_Left_T2\tfocal_left\t%s/TN_Left_T2.dedup.bam\n' \
    "${library_dir}" \
    >> "${t2_tsv}.new"

printf 'TN_Right_T2\tfocal_right\t%s/TN_Right_T2.dedup.bam\n' \
    "${library_dir}" \
    >> "${t2_tsv}.new"

cp "${controls_tsv}.new" "${t3_tsv}.new"

printf 'TN_Left_T3\tfocal_left\t%s/TN_Left_T3.dedup.bam\n' \
    "${library_dir}" \
    >> "${t3_tsv}.new"

printf 'TN_Right_T3\tfocal_right\t%s/TN_Right_T3.dedup.bam\n' \
    "${library_dir}" \
    >> "${t3_tsv}.new"

for set_name in \
    controls \
    merged_T2_T3 \
    T2_only \
    T3_only
do
    tsv="${set_dir}/${set_name}.tsv"

    mv "${tsv}.new" "${tsv}"

    awk \
        -F '\t' \
        'NR > 1 {print $3}' \
        "${tsv}" \
        > "${set_dir}/${set_name}.bamlist.new"

    mv \
        "${set_dir}/${set_name}.bamlist.new" \
        "${set_dir}/${set_name}.bamlist"
done

controls_n=$(
    awk 'END {print NR - 1}' \
        "${controls_tsv}"
)

primary_n=$(
    awk 'END {print NR - 1}' \
        "${primary_tsv}"
)

if [[
    "${controls_n}" -ne 18 ||
    "${primary_n}" -ne 20
]]
then
    printf 'ERROR: expected 18 controls and 20 primary samples; observed %s and %s.\n' \
        "${controls_n}" \
        "${primary_n}" >&2
    false
fi

rm -f "${tmp_dir}"/*.bam
rmdir "${tmp_dir}" 2>/dev/null || true

printf 'FOCAL_MERGE_AND_BAM_SETS_COMPLETE\n'
printf 'controls=%s\n' "${controls_n}"
printf 'primary_samples=%s\n' "${primary_n}"
