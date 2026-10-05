#!/usr/bin/env bash
#SBATCH --job-name=gyn170vc
#SBATCH --partition=cpu
#SBATCH --array=0-20
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=logs/170_call_joint_autosomal_vcfs.%A_%a.out
#SBATCH --error=logs/170_call_joint_autosomal_vcfs.%A_%a.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin

bcftools="${tool_bin}/bcftools"
bgzip="${tool_bin}/bgzip"
tabix="${tool_bin}/tabix"

reference_path_file="${project}/analysis/genomewide_mixing_input/metadata/reference_fasta.path"
chrom_manifest="${project}/analysis/genomewide_mixing_input/metadata/chromosome_calling_manifest.tsv"
set_dir="${project}/analysis/genomewide_mixing_input/metadata/bam_sets"
output_root="${project}/analysis/genomewide_mixing_input/variant_calls/by_chromosome"

mkdir -p \
    "${project}/logs" \
    "${output_root}"

for required in \
    "${bcftools}" \
    "${bgzip}" \
    "${tabix}" \
    "${reference_path_file}" \
    "${chrom_manifest}" \
    "${set_dir}/controls.bamlist" \
    "${set_dir}/merged_T2_T3.bamlist" \
    "${set_dir}/T2_only.bamlist" \
    "${set_dir}/T3_only.bamlist"
do
    if [[ ! -s "${required}" ]]
    then
        printf 'ERROR: required input is missing or empty: %s\n' \
            "${required}" >&2
        false
    fi
done

reference=$(head -n 1 "${reference_path_file}")

if [[
    "${reference}" != *'/GCF_056956725.1/'* ||
    ! -s "${reference}" ||
    ! -s "${reference}.fai"
]]
then
    printf 'ERROR: wrong or unindexed active reference: %s\n' \
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
    <<< "${row}"

chrom_dir="${output_root}/chr${molecule}_${contig}"

mkdir -p "${chrom_dir}"

control_bcf="${chrom_dir}/controls.discovery.biallelic_snps.bcf"
alleles="${chrom_dir}/controls.discovery.alleles.tsv.gz"

if [[
    ! -s "${control_bcf}" ||
    ! -s "${control_bcf}.csi"
]] || ! "${bcftools}" view \
    -h \
    "${control_bcf}" \
    >/dev/null 2>&1
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
        --threads 1 \
        "${control_bcf}.new"

    mv \
        "${control_bcf}.new" \
        "${control_bcf}"

    mv \
        "${control_bcf}.new.csi" \
        "${control_bcf}.csi"
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

mv \
    "${alleles}.new" \
    "${alleles}"

mv \
    "${alleles}.new.tbi" \
    "${alleles}.tbi"

discovery_sites=$(
    "${bcftools}" index \
        -n \
        "${control_bcf}"
)

if [[ "${discovery_sites}" -le 0 ]]
then
    printf 'ERROR: no control-discovered biallelic SNPs on %s.\n' \
        "${contig}" >&2
    false
fi

joint_call()
{
    local set_name=$1
    local bamlist="${set_dir}/${set_name}.bamlist"
    local output_bcf="${chrom_dir}/${set_name}.joint_at_control_sites.bcf"

    if [[
        ! -s "${output_bcf}" ||
        ! -s "${output_bcf}.csi"
    ]] || ! "${bcftools}" view \
        -h \
        "${output_bcf}" \
        >/dev/null 2>&1
    then
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
            -b "${bamlist}" \
            | "${bcftools}" call \
                -Ob \
                --threads 1 \
                -m \
                -A \
                -C alleles \
                -T "${alleles}" \
                -o "${output_bcf}.new"

        "${bcftools}" index \
            -f \
            --threads 1 \
            "${output_bcf}.new"

        mv \
            "${output_bcf}.new" \
            "${output_bcf}"

        mv \
            "${output_bcf}.new.csi" \
            "${output_bcf}.csi"
    fi

    local observed_samples
    local expected_samples
    local observed_sites

    observed_samples=$(
        "${bcftools}" query \
            -l \
            "${output_bcf}" \
            | wc -l
    )

    expected_samples=$(
        wc -l < "${bamlist}"
    )

    observed_sites=$(
        "${bcftools}" index \
            -n \
            "${output_bcf}"
    )

    if [[
        "${observed_samples}" -ne "${expected_samples}" ||
        "${observed_sites}" -ne "${discovery_sites}"
    ]]
    then
        printf 'ERROR: %s validation failed: samples %s/%s; sites %s/%s.\n' \
            "${set_name}" \
            "${observed_samples}" \
            "${expected_samples}" \
            "${observed_sites}" \
            "${discovery_sites}" >&2
        false
    fi

    "${bcftools}" query \
        -f '%CHROM\t%POS\t%REF\t%ALT\n' \
        "${output_bcf}" \
        | md5sum \
        | awk '{print $1}' \
        > "${output_bcf}.site_keys.md5"
}

joint_call merged_T2_T3
joint_call T2_only
joint_call T3_only

primary_md5=$(
    cat \
        "${chrom_dir}/merged_T2_T3.joint_at_control_sites.bcf.site_keys.md5"
)

t2_md5=$(
    cat \
        "${chrom_dir}/T2_only.joint_at_control_sites.bcf.site_keys.md5"
)

t3_md5=$(
    cat \
        "${chrom_dir}/T3_only.joint_at_control_sites.bcf.site_keys.md5"
)

if [[
    "${primary_md5}" != "${t2_md5}" ||
    "${primary_md5}" != "${t3_md5}"
]]
then
    printf 'ERROR: frozen site keys differ among joint datasets on %s.\n' \
        "${contig}" >&2
    false
fi

printf 'array_id\tmolecule\tcontig\tlength\tdiscovery_sites\tprimary_site_md5\treference\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${array_id}" \
    "${molecule}" \
    "${contig}" \
    "${length}" \
    "${discovery_sites}" \
    "${primary_md5}" \
    "${reference}" \
    > "${chrom_dir}/calling.complete.tsv.new"

mv \
    "${chrom_dir}/calling.complete.tsv.new" \
    "${chrom_dir}/calling.complete.tsv"

printf 'JOINT_CHROMOSOME_CALL_COMPLETE\tchr=%s\tcontig=%s\tsites=%s\n' \
    "${molecule}" \
    "${contig}" \
    "${discovery_sites}"
