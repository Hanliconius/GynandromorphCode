#!/usr/bin/env bash

set -u
set -o pipefail

project="/scratch/martinlab/joe/GynandroCalling"

printf '\n============================================================\n'
printf '1. BAND 51-SAMPLE / ORIGINAL-502 ANALYSIS FILES\n'
printf '============================================================\n\n'

find "$project" \
    -type f \
    \( \
        -iname '*51sample*' \
        -o -iname '*51_sample*' \
        -o -iname '*original_502*' \
        -o -iname '*502*metadata*' \
        -o -iname '*phenotype*' \
        -o -iname '*sample*metadata*' \
        -o -iname '*sample*info*' \
        -o -iname '*taxon*' \
    \) \
    2>/dev/null \
    | grep -Ei \
        'optix|band|metadata|phenotype|sample|502|51sample' \
    | sort

printf '\n============================================================\n'
printf '2. TABULAR FILES WHOSE HEADER LOOKS LIKE SAMPLE METADATA\n'
printf '============================================================\n\n'

find "$project" \
    -type f \
    \( \
        -name '*.tsv' \
        -o -name '*.csv' \
        -o -name '*.txt' \
    \) \
    2>/dev/null \
    | while IFS= read -r f
do
    if [[ ! -s "$f" ]]
    then
        continue
    fi

    header="$(
        head -n 1 "$f" 2>/dev/null
    )"

    if printf '%s\n' "$header" \
        | grep -Eqi \
        'sample|accession|species|subspecies|taxon|phenotype|locality|location|country|population|race|morph'
    then
        printf '%s\t%s\n' \
            "$f" \
            "$header"
    fi
done

printf '\n============================================================\n'
printf '3. VCF/BCF FILES WITH BAND REGION + SAMPLE COUNTS\n'
printf '============================================================\n\n'

find "$project" \
    -type f \
    \( \
        -name '*.vcf' \
        -o -name '*.vcf.gz' \
        -o -name '*.bcf' \
    \) \
    2>/dev/null \
    | while IFS= read -r f
do
    # Avoid printing all samples here; just count.
    n="$(
        bcftools query -l "$f" \
            2>/dev/null \
        | awk 'END {print NR+0}'
    )"

    if [[ "$n" -gt 0 ]]
    then
        first="$(
            bcftools query \
                -r Hmr18:1112000-1114000 \
                -f '%CHROM\t%POS\n' \
                "$f" \
                2>/dev/null \
            | head -n 1
        )"

        if [[ -n "$first" ]]
        then
            printf '%s\t%s samples\t%s\n' \
                "$f" \
                "$n" \
                "$first"
        fi
    fi
done

printf '\n============================================================\n'
printf '4. SAMPLE NAMES FROM LIKELY BAND VCFs\n'
printf '============================================================\n\n'

find "$project" \
    -type f \
    \( \
        -iname '*band*.vcf' \
        -o -iname '*band*.vcf.gz' \
        -o -iname '*optix*.vcf' \
        -o -iname '*optix*.vcf.gz' \
    \) \
    2>/dev/null \
    | while IFS= read -r f
do
    n="$(
        bcftools query -l "$f" \
            2>/dev/null \
        | awk 'END {print NR+0}'
    )"

    if [[ "$n" -gt 0 ]]
    then
        printf '\n--- %s (%s samples) ---\n' \
            "$f" \
            "$n"

        bcftools query -l "$f" \
            2>/dev/null \
        | head -n 70
    fi
done
