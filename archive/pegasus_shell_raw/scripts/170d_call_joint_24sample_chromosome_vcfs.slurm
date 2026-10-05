#!/usr/bin/env bash
#SBATCH --job-name=gyn170j24
#SBATCH --partition=cpu
#SBATCH --array=0-20
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/170d_call_joint_24sample_chromosome_vcfs.%A_%a.out
#SBATCH --error=logs/170d_call_joint_24sample_chromosome_vcfs.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin

bcftools="${tool_bin}/bcftools"
bgzip="${tool_bin}/bgzip"
tabix="${tool_bin}/tabix"

reference_path_file="${project}/analysis/genomewide_mixing_input/metadata/reference_fasta.path"
chrom_manifest="${project}/analysis/genomewide_mixing_input/metadata/chromosome_calling_manifest.tsv"
set_dir="${project}/analysis/genomewide_mixing_input/metadata/bam_sets_24"
output_root="${project}/analysis/genomewide_mixing_input/variant_calls/joint_24/by_chromosome"

mkdir -p "${output_root}"

reference=$(head -n 1 "${reference_path_file}")

if [[
    "${reference}" != *'/GCF_056956725.1/'* ||
    ! -s "${reference}" ||
    ! -s "${reference}.fai"
]]
then
    printf 'ERROR: invalid active reference: %s\n' \
        "${reference}" >&2
    false
fi

row=$(
    awk \
        -F '\t' \
        -v id="${SLURM_ARRAY_TASK_ID}" \
        '
            NR > 1 && $1 == id {
                print
                found++
            }

            END {
                if (found != 1) {
                    exit 1
                }
            }
        ' \
        "${chrom_manifest}"
)

IFS=$'\t' read -r \
    array_id \
    molecule \
    contig \
    length \
    chromosome_type \
    <<< "${row}"

chrom_dir="${output_root}/chr${molecule}_${contig}"

mkdir -p "${chrom_dir}"

control_bcf="${chrom_dir}/controls.discovery.biallelic_snps.bcf"
alleles="${chrom_dir}/controls.discovery.alleles.tsv.gz"
joint_bcf="${chrom_dir}/all_24_samples.joint_at_control_sites.bcf"

if [[
    ! -s "${control_bcf}" ||
    ! -s "${control_bcf}.csi"
]]
then
    "${bcftools}" mpileup \
        -Ou \
        --threads 1 \
        -f "${reference}" \
        -r "${contig}" \
        -q 20 \
        -Q 20 \
        -d 1000 \
        -a FORMAT/AD,FORMAT/DP,FORMAT/SP \
        -b "${set_dir}/controls.bamlist" \
        | "${bcftools}" call \
            -Ou \
            -m \
            -v \
        | "${bcftools}" norm \
            -Ou \
            -f "${reference}" \
        | "${bcftools}" view \
            -Ob \
            --threads 1 \
            -m2 \
            -M2 \
            -v snps \
            -o "${control_bcf}.new"

    "${bcftools}" index \
        -f \
        "${control_bcf}.new"

    mv "${control_bcf}.new" "${control_bcf}"
    mv "${control_bcf}.new.csi" "${control_bcf}.csi"
fi

"${bcftools}" query \
    -f '%CHROM\t%POS\t%REF,%ALT\n' \
    "${control_bcf}" \
    | "${bgzip}" -c \
    > "${alleles}.new"

"${tabix}" \
    -f \
    -s 1 \
    -b 2 \
    -e 2 \
    "${alleles}.new"

mv "${alleles}.new" "${alleles}"
mv "${alleles}.new.tbi" "${alleles}.tbi"

control_sites=$(
    "${bcftools}" index \
        -n \
        "${control_bcf}"
)

if [[ "${control_sites}" -le 0 ]]
then
    printf 'ERROR: no control-discovered SNPs on %s\n' \
        "${contig}" >&2
    false
fi

"${bcftools}" mpileup \
    -Ou \
    --threads 1 \
    -f "${reference}" \
    -r "${contig}" \
    -T "${alleles}" \
    -q 20 \
    -Q 20 \
    -d 1000 \
    -a FORMAT/AD,FORMAT/DP,FORMAT/SP \
    -b "${set_dir}/all_24_samples.bamlist" \
    | "${bcftools}" call \
        -Ob \
        --threads 1 \
        -m \
        -A \
        -C alleles \
        -T "${alleles}" \
        -G "${set_dir}/all_24_samples.groups.tsv" \
        -o "${joint_bcf}.new"

"${bcftools}" index \
    -f \
    "${joint_bcf}.new"

sample_count=$(
    "${bcftools}" query \
        -l \
        "${joint_bcf}.new" \
        | wc -l
)

joint_sites=$(
    "${bcftools}" index \
        -n \
        "${joint_bcf}.new"
)

expected_samples="${chrom_dir}/expected.samples"
observed_samples="${chrom_dir}/observed.samples"

awk \
    -F '\t' \
    'NR > 1 {print $1}' \
    "${set_dir}/all_24_samples.tsv" \
    > "${expected_samples}"

"${bcftools}" query \
    -l \
    "${joint_bcf}.new" \
    > "${observed_samples}"

if ! cmp \
    -s \
    "${expected_samples}" \
    "${observed_samples}"
then
    printf 'ERROR: joint VCF sample order is incorrect.\n' >&2
    diff \
        -u \
        "${expected_samples}" \
        "${observed_samples}" >&2 || true
    false
fi

if [[
    "${sample_count}" -ne 24 ||
    "${joint_sites}" -ne "${control_sites}"
]]
then
    printf 'ERROR: expected 24 samples and %s sites; found %s samples and %s sites.\n' \
        "${control_sites}" \
        "${sample_count}" \
        "${joint_sites}" >&2
    false
fi

mv "${joint_bcf}.new" "${joint_bcf}"
mv "${joint_bcf}.new.csi" "${joint_bcf}.csi"

site_md5=$(
    "${bcftools}" query \
        -f '%CHROM\t%POS\t%REF\t%ALT\n' \
        "${joint_bcf}" \
        | md5sum \
        | awk '{print $1}'
)

sample_md5=$(
    "${bcftools}" query \
        -l \
        "${joint_bcf}" \
        | md5sum \
        | awk '{print $1}'
)

printf 'array_id\tmolecule\tcontig\tchromosome_type\tsamples\tsites\tsite_md5\tsample_md5\n%s\t%s\t%s\t%s\t24\t%s\t%s\t%s\n' \
    "${array_id}" \
    "${molecule}" \
    "${contig}" \
    "${chromosome_type}" \
    "${joint_sites}" \
    "${site_md5}" \
    "${sample_md5}" \
    > "${chrom_dir}/calling.complete.tsv.new"

mv \
    "${chrom_dir}/calling.complete.tsv.new" \
    "${chrom_dir}/calling.complete.tsv"

printf 'JOINT_24_SAMPLE_CALL_COMPLETE\tchromosome=%s\ttype=%s\tsites=%s\n' \
    "${molecule}" \
    "${chromosome_type}" \
    "${joint_sites}"
