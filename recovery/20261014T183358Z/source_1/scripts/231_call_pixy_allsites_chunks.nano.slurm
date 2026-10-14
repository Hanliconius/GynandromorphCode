#!/usr/bin/env bash
#SBATCH --job-name=piAllSite
#SBATCH --partition=nano
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/231_pixy_allsites.%A_%a.out
#SBATCH --error=logs/231_pixy_allsites.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
root="${project}/analysis/pixy_pi_17controls_allsites"
metadata="${root}/metadata"

tool_bin="${HOME}/.local/share/mamba/envs/colias_w_popgen/bin"
bcftools="${tool_bin}/bcftools"
tabix="${tool_bin}/tabix"

reference=$(head -n 1 \
    "${project}/analysis/genomewide_mixing_input/metadata/reference_fasta.path")

chunks="${metadata}/chunks.tsv"
bamlist="${metadata}/controls17.bamlist"
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

outdir="${root}/calls/chr${molecule}_${contig}"
mkdir -p "${outdir}"

vcf="${outdir}/chunk_${start}_${end}.allsites.dp10.vcf.gz"
complete="${outdir}/chunk_${start}_${end}.complete.tsv"

if [[
    -s "${vcf}" &&
    -s "${vcf}.tbi" &&
    -s "${complete}"
]] && "${bcftools}" view -h "${vcf}" >/dev/null 2>&1
then
    echo "ALREADY_COMPLETE task=${task_id} region=${region}"
    cat "${complete}"
    exit 0
fi

tmp="${outdir}/.chunk_${start}_${end}.${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}.tmp.vcf.gz"

rm -f \
    "${tmp}" \
    "${tmp}.tbi" \
    "${complete}.new"

echo "CALL_START task=${task_id} region=${region} host=$(hostname) time=$(date --iso-8601=seconds)"

###############################################################################
# Important:
#   - joint call across all 17 individuals
#   - autosomes explicitly diploid
#   - -I skips indel discovery: we want invariant sites + SNPs
#   - no -v: homozygous-reference/invariant records are retained
#   - population groups prevent one HWE prior being imposed across all taxa
#   - GT with DP < 10 becomes missing
#   - keep only GT and DP to keep the all-sites files manageable
###############################################################################

"${bcftools}" mpileup \
    -Ou \
    --threads 1 \
    -f "${reference}" \
    -r "${region}" \
    -q 20 \
    -Q 20 \
    -d 1000 \
    -I \
    -a FORMAT/AD,FORMAT/DP \
    -b "${bamlist}" \
| "${bcftools}" call \
    -Ou \
    --threads 1 \
    -m \
    --ploidy 2 \
    -G "${populations}" \
| "${bcftools}" view \
    -Ou \
    -m1 \
    -M2 \
| "${bcftools}" filter \
    -Ou \
    -S . \
    -e 'FMT/DP<10' \
| "${bcftools}" annotate \
    -Oz \
    --threads 2 \
    -x INFO,FORMAT/PL,FORMAT/AD \
    -o "${tmp}"

"${tabix}" \
    -f \
    -p vcf \
    "${tmp}"

expected_samples="${outdir}/.expected_samples.${SLURM_ARRAY_TASK_ID}.txt"
observed_samples="${outdir}/.observed_samples.${SLURM_ARRAY_TASK_ID}.txt"

cut -f1 "${populations}" \
    > "${expected_samples}"

"${bcftools}" query -l "${tmp}" \
    > "${observed_samples}"

if ! cmp -s "${expected_samples}" "${observed_samples}"
then
    echo "ERROR: sample order mismatch in ${region}" >&2
    diff -u "${expected_samples}" "${observed_samples}" >&2 || true
    rm -f "${expected_samples}" "${observed_samples}"
    exit 1
fi

rm -f "${expected_samples}" "${observed_samples}"

read -r records invariant variant < <(
    "${bcftools}" query \
        -f '%ALT\n' \
        "${tmp}" \
    | awk '
        {
            n++
            if ($1 == ".")
                inv++
            else
                var++
        }
        END {
            print n + 0, inv + 0, var + 0
        }
    '
)

if [[ "${records}" -le 0 ]]
then
    echo "ERROR: zero records in ${region}" >&2
    exit 1
fi

if [[ "${invariant}" -le 0 ]]
then
    echo "ERROR: no invariant records in ${region}; refusing SNP-only callset" >&2
    exit 1
fi

chunk_length=$((end - start + 1))

fraction=$(
    awk \
        -v n="${records}" \
        -v d="${chunk_length}" \
        'BEGIN {printf "%.6f", n / d}'
)

mv "${tmp}" "${vcf}"
mv "${tmp}.tbi" "${vcf}.tbi"

printf \
    'task_id\tmolecule\tcontig\tstart\tend\tchunk_bp\trecords\tinvariant_records\tvariant_records\trecord_fraction\tsamples\tmin_DP\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t17\t10\n' \
    "${task_id}" \
    "${molecule}" \
    "${contig}" \
    "${start}" \
    "${end}" \
    "${chunk_length}" \
    "${records}" \
    "${invariant}" \
    "${variant}" \
    "${fraction}" \
    > "${complete}.new"

mv "${complete}.new" "${complete}"

echo "CALL_COMPLETE task=${task_id} region=${region} records=${records} invariant=${invariant} variant=${variant} record_fraction=${fraction}"
