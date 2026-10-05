#!/usr/bin/env bash
set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
FASTQ="${ROOT}/analysis/genomewide_mixing_input/fastq"
OUT="${ROOT}/analysis/smudgeplot_arrangements"

mkdir -p "${OUT}"

link_reads () {
    local sample="$1"
    local prefix="$2"
    local dest="$3"

    local r1="${FASTQ}/${sample}/${sample}_1.fastq.gz"
    local r2="${FASTQ}/${sample}/${sample}_2.fastq.gz"

    [[ -s "${r1}" ]] || { echo "MISSING ${r1}"; return 1; }
    [[ -s "${r2}" ]] || { echo "MISSING ${r2}"; return 1; }

    #
    # Deliberately name links *.fastq.gz.
    #
    ln -sfn "${r1}" "${dest}/${prefix}_R1.fastq.gz"
    ln -sfn "${r2}" "${dest}/${prefix}_R2.fastq.gz"
}


make_group () {
    local name="$1"
    shift

    local dest="${OUT}/${name}"

    mkdir -p "${dest}"

    local i=1

    for sample in "$@"; do
        link_reads \
            "${sample}" \
            "$(printf '%02d_%s' "${i}" "${sample}")" \
            "${dest}"

        i=$((i + 1))
    done
}


# Left alone = both Left sequencing libraries
make_group \
    TN_Left \
    TN_Left_T2 \
    TN_Left_T3

# Right alone = both Right sequencing libraries
make_group \
    TN_Right \
    TN_Right_T2 \
    TN_Right_T3

# T2 merge = Left + Right from T2
make_group \
    TN_T2_merge \
    TN_Left_T2 \
    TN_Right_T2

# T3 merge = Left + Right from T3
make_group \
    TN_T3_merge \
    TN_Left_T3 \
    TN_Right_T3

# Normal diploid candidate
make_group \
    ERR2298208_single \
    ERR2298208

# Artificial two-individual merge = pseudo-tetraploid benchmark
make_group \
    pseudoTetraploid_ERR2298208_ERR1307111 \
    ERR2298208 \
    ERR1307111


MANIFEST="${OUT}/smudgeplot_array_manifest.tsv"

{
    printf "array_id\tlabel\tdirectory\n"

    printf "0\tTN_Left\t%s\n" \
        "${OUT}/TN_Left"

    printf "1\tTN_Right\t%s\n" \
        "${OUT}/TN_Right"

    printf "2\tTN_T2_merge\t%s\n" \
        "${OUT}/TN_T2_merge"

    printf "3\tTN_T3_merge\t%s\n" \
        "${OUT}/TN_T3_merge"

    printf "4\tERR2298208_single\t%s\n" \
        "${OUT}/ERR2298208_single"

    printf "5\tpseudoTetraploid_ERR2298208_ERR1307111\t%s\n" \
        "${OUT}/pseudoTetraploid_ERR2298208_ERR1307111"

} > "${MANIFEST}"


echo
echo "MANIFEST"
column -t -s $'\t' "${MANIFEST}"

echo
echo "INPUT COUNTS"

while IFS=$'\t' read -r id label dir; do

    [[ "${id}" == "array_id" ]] && continue

    printf "%-45s " "${label}"

    find "${dir}" \
        -maxdepth 1 \
        -type l \
        -name '*.fastq.gz' \
        | wc -l

done < "${MANIFEST}"
