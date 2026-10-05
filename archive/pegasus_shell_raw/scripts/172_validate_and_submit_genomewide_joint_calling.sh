#!/usr/bin/env bash

set -u

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
python_bin="${tool_bin}/python"
metadata="${project}/analysis/genomewide_mixing_input/metadata"
resolved_manifest="${metadata}/resolved_download_manifest.tsv"
reference_path_file="${metadata}/reference_fasta.path"

validation_status=0

if ! cd "${project}"
then
    printf 'FAIL\tcannot enter project directory: %s\n' \
        "${project}"
    validation_status=1
fi

mkdir -p logs

printf '%s\n' 'SCRIPT SYNTAX'

for script in \
    scripts/167_align_genomewide_libraries.slurm \
    scripts/168_merge_focal_sides_and_build_bam_sets.slurm \
    scripts/170_call_joint_autosomal_vcfs.slurm \
    scripts/171_concat_and_validate_joint_autosomal_vcfs.slurm
do
    if bash -n "${script}"
    then
        printf 'PASS\t%s\n' "${script}"
    else
        printf 'FAIL\t%s\n' "${script}"
        validation_status=1
    fi
done

if "${python_bin}" -m py_compile \
    scripts/169_prepare_autosome_calling_manifest.py
then
    printf 'PASS\t%s\n' \
        scripts/169_prepare_autosome_calling_manifest.py
else
    printf 'FAIL\t%s\n' \
        scripts/169_prepare_autosome_calling_manifest.py
    validation_status=1
fi

printf '\n%s\n' 'REQUIRED TOOLS'

for tool in \
    bwa \
    samtools \
    bcftools \
    bgzip \
    tabix \
    python
do
    if [[ -x "${tool_bin}/${tool}" ]]
    then
        printf 'PASS\t%s\t%s\n' \
            "${tool}" \
            "${tool_bin}/${tool}"
    else
        printf 'FAIL\t%s\t%s\n' \
            "${tool}" \
            "${tool_bin}/${tool}"
        validation_status=1
    fi
done

if [[ -s "${resolved_manifest}" ]]
then
    manifest_rows=$(
        awk 'END {print NR - 1}' \
            "${resolved_manifest}"
    )
else
    manifest_rows=0
fi

printf '\nresolved_manifest_rows=%s\n' \
    "${manifest_rows}"

if [[ "${manifest_rows}" -ne 22 ]]
then
    validation_status=1
fi

if [[ -s "${reference_path_file}" ]]
then
    reference=$(head -n 1 "${reference_path_file}")
else
    reference=""
fi

printf 'active_reference=%s\n' \
    "${reference}"

if [[
    "${reference}" != *'/GCF_056956725.1/'* ||
    ! -s "${reference}" ||
    ! -s "${reference}.fai"
]]
then
    printf 'FAIL\tactive reference must be indexed GCF_056956725.1\n'
    validation_status=1
else
    printf 'PASS\tactive reference is indexed GCF_056956725.1\n'
fi

if [[ "${validation_status}" -eq 0 ]]
then
    "${python_bin}" \
        scripts/169_prepare_autosome_calling_manifest.py \
        > logs/169_prepare_autosome_calling_manifest.log \
        2>&1

    manifest_status=$?

    cat logs/169_prepare_autosome_calling_manifest.log
else
    manifest_status=127
fi

autosome_manifest="${metadata}/autosome_calling_manifest.tsv"

if [[
    "${manifest_status}" -eq 0 &&
    -s "${autosome_manifest}"
]]
then
    autosome_rows=$(
        awk 'END {print NR - 1}' \
            "${autosome_manifest}"
    )
else
    autosome_rows=0
fi

printf 'autosome_manifest_rows=%s\n' \
    "${autosome_rows}"

if [[ "${autosome_rows}" -ne 20 ]]
then
    validation_status=1
fi

if grep \
    -Eq \
    'bcftools.*merge|\$\{bcftools\}.*merge' \
    scripts/170_call_joint_autosomal_vcfs.slurm \
    scripts/171_concat_and_validate_joint_autosomal_vcfs.slurm
then
    printf 'FAIL\tbcftools merge detected in joint-calling workflow\n'
    validation_status=1
else
    printf 'PASS\tno per-sample bcftools merge operation\n'
fi

if grep \
    -q \
    'controls.bamlist' \
    scripts/170_call_joint_autosomal_vcfs.slurm \
    && grep \
        -q \
        -- '-C alleles' \
        scripts/170_call_joint_autosomal_vcfs.slurm \
    && grep \
        -q \
        'concat' \
        scripts/171_concat_and_validate_joint_autosomal_vcfs.slurm
then
    printf 'PASS\tcontrol-only discovery, constrained joint calling, and chromosome concatenation detected\n'
else
    printf 'FAIL\tstatic joint-calling audit failed\n'
    validation_status=1
fi

printf '\nFINAL PREFLIGHT STATUS\n'
printf 'validation_status=%s\n' \
    "${validation_status}"
printf 'manifest_status=%s\n' \
    "${manifest_status}"
printf 'autosome_rows=%s\n' \
    "${autosome_rows}"

if [[
    "${validation_status}" -eq 0 &&
    "${manifest_status}" -eq 0 &&
    "${autosome_rows}" -eq 20
]]
then
    align_submission=$(
        sbatch \
            --parsable \
            --array=0-21 \
            scripts/167_align_genomewide_libraries.slurm
    )
    align_job=${align_submission%%;*}

    merge_submission=$(
        sbatch \
            --parsable \
            --dependency="afterok:${align_job}" \
            scripts/168_merge_focal_sides_and_build_bam_sets.slurm
    )
    merge_job=${merge_submission%%;*}

    call_submission=$(
        sbatch \
            --parsable \
            --dependency="afterok:${merge_job}" \
            --array=0-19 \
            scripts/170_call_joint_autosomal_vcfs.slurm
    )
    call_job=${call_submission%%;*}

    concat_submission=$(
        sbatch \
            --parsable \
            --dependency="afterok:${call_job}" \
            scripts/171_concat_and_validate_joint_autosomal_vcfs.slurm
    )
    concat_job=${concat_submission%%;*}

    printf '\nSUBMITTED JOB CHAIN\n'
    printf 'alignment_array\t%s\t22 simultaneous tasks\n' \
        "${align_job}"
    printf 'focal_merge\t%s\tafterok:%s\n' \
        "${merge_job}" \
        "${align_job}"
    printf 'autosome_call_array\t%s\t20 simultaneous tasks; afterok:%s\n' \
        "${call_job}" \
        "${merge_job}"
    printf 'genomewide_concat\t%s\tafterok:%s\n' \
        "${concat_job}" \
        "${call_job}"

    printf '\nQueue command:\n'
    printf 'squeue -j %s,%s,%s,%s\n' \
        "${align_job}" \
        "${merge_job}" \
        "${call_job}" \
        "${concat_job}"
else
    printf '\nNo jobs submitted because preflight failed.\n'
fi
