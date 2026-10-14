#!/usr/bin/env bash
#SBATCH --job-name=pixyPi
#SBATCH --partition=nano
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=00:30:00
#SBATCH --output=logs/232_pixy_pi.%A_%a.out
#SBATCH --error=logs/232_pixy_pi.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
root="${project}/analysis/pixy_pi_17controls_allsites"
metadata="${root}/metadata"

pixy="${HOME}/.local/share/mamba/envs/pixy/bin/pixy"

chunks="${metadata}/chunks.tsv"
populations="${metadata}/populations.tsv"

row=$(
    awk \
        -F '\t' \
        -v id="${SLURM_ARRAY_TASK_ID}" \
        '
        NR > 1 && $1 == id {
            print
            n++
        }
        END {
            if (n != 1)
                exit 1
        }
        ' \
        "${chunks}"
)

IFS=$'\t' read -r \
    task_id \
    molecule \
    contig \
    chromosome_length \
    start \
    end \
    region \
    <<< "${row}"

call_dir="${root}/calls/chr${molecule}_${contig}"
vcf="${call_dir}/chunk_${start}_${end}.allsites.dp10.vcf.gz"
call_complete="${call_dir}/chunk_${start}_${end}.complete.tsv"

if [[
    ! -s "${vcf}" ||
    ! -s "${vcf}.tbi" ||
    ! -s "${call_complete}"
]]
then
    echo "ERROR: validated all-sites input missing for task ${task_id}: ${vcf}" >&2
    exit 1
fi

invariant=$(
    awk \
        -F '\t' \
        'NR == 2 {print $8}' \
        "${call_complete}"
)

if [[ -z "${invariant}" || "${invariant}" -le 0 ]]
then
    echo "ERROR: task ${task_id} has no validated invariant records" >&2
    exit 1
fi

printf -v task_tag '%04d' "${task_id}"

outdir="${root}/pixy/by_chunk/task_${task_tag}"
prefix="task_${task_tag}"
pi_file="${outdir}/${prefix}_pi.txt"
complete="${outdir}/COMPLETE.tsv"

mkdir -p "${outdir}"

if [[ -s "${pi_file}" && -s "${complete}" ]]
then
    echo "ALREADY_COMPLETE task=${task_id} region=${region}"
    exit 0
fi

rm -f \
    "${pi_file}" \
    "${complete}.new"

echo "PIXY_START task=${task_id} region=${region} host=$(hostname) time=$(date --iso-8601=seconds)"
"${pixy}" --version

###############################################################################
# Deliberately NO --bypass_invariant_check.
###############################################################################

"${pixy}" \
    --stats pi \
    --vcf "${vcf}" \
    --populations "${populations}" \
    --chromosomes "${contig}" \
    --interval_start "${start}" \
    --interval_end "${end}" \
    --window_size 100000 \
    --n_cores "${SLURM_CPUS_PER_TASK}" \
    --output_folder "${outdir}" \
    --output_prefix "${prefix}"

if [[ ! -s "${pi_file}" ]]
then
    echo "ERROR: pixy produced no pi output for task ${task_id}" >&2
    exit 1
fi

header=$(
    head -n 1 "${pi_file}"
)

if [[
    "${header}" != *"count_diffs"* ||
    "${header}" != *"count_comparisons"*
]]
then
    echo "ERROR: unexpected pixy pi header:" >&2
    echo "${header}" >&2
    exit 1
fi

rows=$(
    awk 'END {print NR - 1}' \
        "${pi_file}"
)

if [[ "${rows}" -le 0 ]]
then
    echo "ERROR: no pixy result rows for task ${task_id}" >&2
    exit 1
fi

printf \
    'task_id\tmolecule\tcontig\tstart\tend\tresult_rows\tpi_file\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${task_id}" \
    "${molecule}" \
    "${contig}" \
    "${start}" \
    "${end}" \
    "${rows}" \
    "${pi_file}" \
    > "${complete}.new"

mv "${complete}.new" "${complete}"

echo "PIXY_COMPLETE task=${task_id} region=${region} rows=${rows}"
