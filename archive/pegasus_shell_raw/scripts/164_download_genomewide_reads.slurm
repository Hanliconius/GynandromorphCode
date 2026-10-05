#!/usr/bin/env bash
#SBATCH --job-name=gyn164dl
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/164_download_genomewide_reads.%A_%a.out
#SBATCH --error=logs/164_download_genomewide_reads.%A_%a.err

set -euo pipefail

tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
export PATH="${tool_bin}:/usr/bin:/bin:${PATH}"

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
root="${project}/analysis/genomewide_mixing_input"
manifest="${root}/metadata/download_manifest.tsv"
task_id=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is unset}

cd "${project}"

for tool in prefetch fasterq-dump pigz gzip awk find
do
    command -v "${tool}" >/dev/null 2>&1 || {
        printf 'Required command is unavailable: %s\n' \
            "${tool}" >&2
        false
    }
done

field_value() {
    local wanted=$1

    awk -F '\t' \
        -v target="${task_id}" \
        -v wanted="${wanted}" '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "array_id") {
                    array_col = i
                }

                if ($i == wanted) {
                    wanted_col = i
                }
            }

            next
        }

        $array_col == target {
            print $wanted_col
        }
    ' "${manifest}"
}

accession=$(field_value accession)
library_id=$(field_value library_id)
role=$(field_value role)
side=$(field_value side)
leg=$(field_value leg)
merge_group=$(field_value merge_group)

if [[ ! "${accession}" =~ ^(SRX|SRR|ERR|DRR)[0-9]+$ ]]
then
    printf 'Invalid accession for task %s: %s\n' \
        "${task_id}" \
        "${accession}" >&2
    false
fi

if [[ ! "${library_id}" =~ ^[A-Za-z0-9_.-]+$ ]]
then
    printf 'Unsafe library ID: %s\n' \
        "${library_id}" >&2
    false
fi

sra_dir="${root}/sra/${library_id}"
fastq_dir="${root}/fastq/${library_id}"
tmp_dir="${root}/tmp/fasterq_${library_id}.${SLURM_JOB_ID}.${task_id}"
mapping_dir="${root}/metadata/resolved_runs"

complete="${fastq_dir}/download.complete"
mapping="${mapping_dir}/$(printf '%03d' "${task_id}").${library_id}.tsv"

mkdir -p \
    "${sra_dir}" \
    "${fastq_dir}" \
    "${tmp_dir}/scratch" \
    "${mapping_dir}"

already_complete=0

if [[ -s "${complete}" && -s "${mapping}" ]]
then
    shopt -s nullglob
    existing_fastq=("${fastq_dir}"/*.fastq.gz)
    shopt -u nullglob

    if [[ ${#existing_fastq[@]} -gt 0 ]]
    then
        gzip -t "${existing_fastq[@]}"
        already_complete=1
    fi
fi

if [[ "${already_complete}" -eq 1 ]]
then
    printf 'ALREADY COMPLETE\n'
    cat "${mapping}"
else
    prefetch \
        --max-size u \
        --output-directory "${sra_dir}" \
        "${accession}"

    mapfile -t sra_files < <(
        find "${sra_dir}" \
            -type f \
            \( \
                -name 'SRR*.sra' \
                -o -name 'ERR*.sra' \
                -o -name 'DRR*.sra' \
            \) \
            -print \
            | sort
    )

    if [[ ${#sra_files[@]} -ne 1 ]]
    then
        printf 'Expected one run for %s; found %s\n' \
            "${accession}" \
            "${#sra_files[@]}" >&2

        printf '%s\n' "${sra_files[@]}" >&2
        false
    fi

    sra_file=${sra_files[0]}
    resolved_run=$(basename "${sra_file}" .sra)

    if [[ ! "${resolved_run}" =~ ^(SRR|ERR|DRR)[0-9]+$ ]]
    then
        printf 'Invalid resolved run: %s\n' \
            "${resolved_run}" >&2
        false
    fi

    if command -v vdb-validate >/dev/null 2>&1
    then
        vdb-validate "${sra_file}"
    fi

    fasterq-dump \
        "${sra_file}" \
        --split-files \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --temp "${tmp_dir}/scratch" \
        --outdir "${tmp_dir}"

    mapfile -t mate1_candidates < <(
        find "${tmp_dir}" \
            -maxdepth 1 \
            -type f \
            -name '*_1.fastq' \
            -print \
            | sort
    )

    mapfile -t mate2_candidates < <(
        find "${tmp_dir}" \
            -maxdepth 1 \
            -type f \
            -name '*_2.fastq' \
            -print \
            | sort
    )

    mapfile -t single_candidates < <(
        find "${tmp_dir}" \
            -maxdepth 1 \
            -type f \
            -name '*.fastq' \
            ! -name '*_1.fastq' \
            ! -name '*_2.fastq' \
            -print \
            | sort
    )

    if [[ \
        ${#mate1_candidates[@]} -eq 1 \
        && ${#mate2_candidates[@]} -eq 1 \
        && ${#single_candidates[@]} -eq 0 \
    ]]
    then
        layout=PAIRED
        raw1=${mate1_candidates[0]}
        raw2=${mate2_candidates[0]}

        pigz -p "${SLURM_CPUS_PER_TASK}" \
            "${raw1}" \
            "${raw2}"

        gzip -t \
            "${raw1}.gz" \
            "${raw2}.gz"

        fastq1="${fastq_dir}/${library_id}_1.fastq.gz"
        fastq2="${fastq_dir}/${library_id}_2.fastq.gz"

        mv "${raw1}.gz" "${fastq1}"
        mv "${raw2}.gz" "${fastq2}"

    elif [[ \
        ${#single_candidates[@]} -eq 1 \
        && ${#mate1_candidates[@]} -eq 0 \
        && ${#mate2_candidates[@]} -eq 0 \
    ]]
    then
        layout=SINGLE
        raw1=${single_candidates[0]}

        pigz -p "${SLURM_CPUS_PER_TASK}" \
            "${raw1}"

        gzip -t "${raw1}.gz"

        fastq1="${fastq_dir}/${library_id}.fastq.gz"
        fastq2=

        mv "${raw1}.gz" "${fastq1}"

    elif [[ \
        ${#mate1_candidates[@]} -eq 1 \
        && ${#mate2_candidates[@]} -eq 0 \
        && ${#single_candidates[@]} -eq 0 \
    ]]
    then
        layout=SINGLE
        raw1=${mate1_candidates[0]}

        pigz -p "${SLURM_CPUS_PER_TASK}" \
            "${raw1}"

        gzip -t "${raw1}.gz"

        fastq1="${fastq_dir}/${library_id}.fastq.gz"
        fastq2=

        mv "${raw1}.gz" "${fastq1}"
    else
        printf 'Ambiguous FASTQ output for %s (%s)\n' \
            "${accession}" \
            "${resolved_run}" >&2

        find "${tmp_dir}" \
            -maxdepth 1 \
            -type f \
            -printf '%f\n' >&2

        false
    fi

    if [[ "${layout}" == "PAIRED" ]]
    then
        gzip -t "${fastq1}" "${fastq2}"
    else
        gzip -t "${fastq1}"
    fi

    mapping_tmp="${mapping}.part.${SLURM_JOB_ID}"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
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
        fastq2 \
        > "${mapping_tmp}"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${task_id}" \
        "${accession}" \
        "${resolved_run}" \
        "${library_id}" \
        "${role}" \
        "${side}" \
        "${leg}" \
        "${merge_group}" \
        "${layout}" \
        "${fastq1}" \
        "${fastq2}" \
        >> "${mapping_tmp}"

    mv "${mapping_tmp}" "${mapping}"

    date --iso-8601=seconds \
        > "${complete}"

    find "${tmp_dir}" \
        -type f \
        -delete

    find "${tmp_dir}" \
        -depth \
        -type d \
        -empty \
        -delete

    printf 'DOWNLOAD COMPLETE\n'
    cat "${mapping}"
fi
