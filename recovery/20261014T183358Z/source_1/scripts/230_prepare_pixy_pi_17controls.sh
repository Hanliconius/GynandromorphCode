#!/usr/bin/env bash

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
root="${project}/analysis/pixy_pi_17controls_allsites"
metadata="${root}/metadata"

source_controls="${project}/analysis/genomewide_mixing_input/metadata/bam_sets_24/controls.tsv"
chrom_manifest="${project}/analysis/genomewide_mixing_input/metadata/chromosome_calling_manifest.tsv"
reference_path_file="${project}/analysis/genomewide_mixing_input/metadata/reference_fasta.path"

tool_bin="${HOME}/.local/share/mamba/envs/colias_w_popgen/bin"
bcftools="${tool_bin}/bcftools"
samtools="${tool_bin}/samtools"
pixy="${HOME}/.local/share/mamba/envs/pixy/bin/pixy"

mkdir -p \
    "${metadata}" \
    "${root}/calls" \
    "${root}/pixy/by_chunk" \
    "${root}/summary" \
    "${project}/logs"

for f in \
    "${source_controls}" \
    "${chrom_manifest}" \
    "${reference_path_file}"
do
    if [[ ! -s "${f}" ]]
    then
        echo "ERROR: missing required input: ${f}" >&2
        exit 1
    fi
done

for x in \
    "${bcftools}" \
    "${samtools}" \
    "${pixy}"
do
    if [[ ! -x "${x}" ]]
    then
        echo "ERROR: missing executable: ${x}" >&2
        exit 1
    fi
done

reference=$(head -n 1 "${reference_path_file}")

if [[
    "${reference}" != *'/GCF_056956725.1/'* ||
    ! -s "${reference}" ||
    ! -s "${reference}.fai"
]]
then
    echo "ERROR: invalid reference: ${reference}" >&2
    exit 1
fi

###############################################################################
# Exact 17 biological controls retained for population pi.
# Excludes:
#   SRR11243339 = H. melpomene x H. cydno hybrid
#   all TN gynandromorph representations
###############################################################################

{
    printf '%s\t%s\n' ERR2298201 H_m_melpomene
    printf '%s\t%s\n' ERR2298216 H_m_melpomene
    printf '%s\t%s\n' ERR3790122 H_m_bellula
    printf '%s\t%s\n' ERR3790123 H_m_bellula
    printf '%s\t%s\n' ERR3790118 H_m_bellula
    printf '%s\t%s\n' ERR3790124 H_m_bellula
    printf '%s\t%s\n' ERR2298257 H_m_amaryllis
    printf '%s\t%s\n' ERR2298260 H_m_amaryllis
    printf '%s\t%s\n' ERR1307113 H_m_malleti
    printf '%s\t%s\n' ERR1307112 H_m_malleti
    printf '%s\t%s\n' ERR1307111 H_m_malleti
    printf '%s\t%s\n' ERR2298208 H_m_rosina
    printf '%s\t%s\n' ERR2298221 H_m_rosina
    printf '%s\t%s\n' ERR2298265 H_timareta_thelxinoe
    printf '%s\t%s\n' ERR2298263 H_timareta_thelxinoe
    printf '%s\t%s\n' ERR2298235 H_timareta_NG2008
    printf '%s\t%s\n' ERR2298236 H_timareta_NG2008
} > "${metadata}/populations.tsv.new"

mv \
    "${metadata}/populations.tsv.new" \
    "${metadata}/populations.tsv"

###############################################################################
# Resolve BAMs from the already validated 18-control calling set.
###############################################################################

printf 'sample\tpopulation\tbam\n' \
    > "${metadata}/samples.tsv.new"

: > "${metadata}/controls17.bamlist.new"

while IFS=$'\t' read -r sample population
do
    matches=$(
        awk \
            -F '\t' \
            -v s="${sample}" \
            'NR > 1 && $1 == s {n++} END {print n + 0}' \
            "${source_controls}"
    )

    if [[ "${matches}" -ne 1 ]]
    then
        echo "ERROR: sample ${sample} occurs ${matches} times in ${source_controls}" >&2
        exit 1
    fi

    bam=$(
        awk \
            -F '\t' \
            -v s="${sample}" \
            'NR > 1 && $1 == s {print $3}' \
            "${source_controls}"
    )

    if [[
        ! -s "${bam}" ||
        ! -s "${bam}.bai"
    ]] || ! "${samtools}" quickcheck -q "${bam}"
    then
        echo "ERROR: invalid BAM for ${sample}: ${bam}" >&2
        exit 1
    fi

    printf '%s\t%s\t%s\n' \
        "${sample}" \
        "${population}" \
        "${bam}" \
        >> "${metadata}/samples.tsv.new"

    printf '%s\n' \
        "${bam}" \
        >> "${metadata}/controls17.bamlist.new"

done < "${metadata}/populations.tsv"

mv \
    "${metadata}/samples.tsv.new" \
    "${metadata}/samples.tsv"

mv \
    "${metadata}/controls17.bamlist.new" \
    "${metadata}/controls17.bamlist"

sample_n=$(
    wc -l < "${metadata}/populations.tsv"
)

if [[ "${sample_n}" -ne 17 ]]
then
    echo "ERROR: expected 17 controls, found ${sample_n}" >&2
    exit 1
fi

if grep -q '^SRR11243339' "${metadata}/populations.tsv"
then
    echo "ERROR: hybrid SRR11243339 was not excluded" >&2
    exit 1
fi

###############################################################################
# Validate 20 autosomes and construct 1-Mb non-overlapping chunks.
# Chunk size is exactly divisible by the later 100-kb pixy window size.
###############################################################################

autosome_n=$(
    awk \
        -F '\t' \
        'NR > 1 && $2 ~ /^[0-9]+$/ && $2 >= 1 && $2 <= 20 {n++}
         END {print n + 0}' \
        "${chrom_manifest}"
)

autosome_bp=$(
    awk \
        -F '\t' \
        'NR > 1 && $2 ~ /^[0-9]+$/ && $2 >= 1 && $2 <= 20 {s += $4}
         END {printf "%.0f\n", s}' \
        "${chrom_manifest}"
)

if [[ "${autosome_n}" -ne 20 ]]
then
    echo "ERROR: expected 20 autosomes; found ${autosome_n}" >&2
    exit 1
fi

if [[ "${autosome_bp}" -ne 275202915 ]]
then
    echo "ERROR: autosomal bp=${autosome_bp}; expected 275202915" >&2
    exit 1
fi

awk \
    -F '\t' \
    -v OFS='\t' \
    -v chunk=1000000 \
    '
    BEGIN {
        print "task_id","molecule","contig","chromosome_length","start","end","region"
        task = 0
    }

    NR > 1 && $2 ~ /^[0-9]+$/ && $2 >= 1 && $2 <= 20 {
        molecule = $2
        contig = $3
        chrom_len = $4

        for (start = 1; start <= chrom_len; start += chunk) {
            end = start + chunk - 1
            if (end > chrom_len)
                end = chrom_len

            print task,molecule,contig,chrom_len,start,end,contig ":" start "-" end
            task++
        }
    }
    ' \
    "${chrom_manifest}" \
    > "${metadata}/chunks.tsv.new"

mv \
    "${metadata}/chunks.tsv.new" \
    "${metadata}/chunks.tsv"

chunk_n=$(
    awk 'END {print NR - 1}' \
        "${metadata}/chunks.tsv"
)

chunk_bp=$(
    awk \
        -F '\t' \
        'NR > 1 {s += $6 - $5 + 1}
         END {printf "%.0f\n", s}' \
        "${metadata}/chunks.tsv"
)

if [[ "${chunk_bp}" -ne "${autosome_bp}" ]]
then
    echo "ERROR: chunk coverage ${chunk_bp} != autosomal bp ${autosome_bp}" >&2
    exit 1
fi

###############################################################################
# Tool/version preflight.
###############################################################################

echo
echo "===== TOOL VERSIONS ====="
"${bcftools}" --version | head -n 1
"${pixy}" --version

if ! "${pixy}" --help 2>&1 | grep -q -- '--interval_start'
then
    echo "ERROR: pixy lacks --interval_start" >&2
    exit 1
fi

if ! "${pixy}" --help 2>&1 | grep -q -- '--interval_end'
then
    echo "ERROR: pixy lacks --interval_end" >&2
    exit 1
fi

echo
echo "===== PIXY PI PREPARATION COMPLETE ====="
echo "controls=${sample_n}"
echo "populations=$(cut -f2 "${metadata}/populations.tsv" | sort -u | wc -l)"
echo "autosomes=${autosome_n}"
echo "autosomal_bp=${autosome_bp}"
echo "chunk_size_bp=1000000"
echo "chunk_tasks=${chunk_n}"
echo
echo "POPULATIONS"
awk -F '\t' '{n[$2]++} END {for (p in n) print p "\t" n[p]}' \
    "${metadata}/populations.tsv" \
    | sort
