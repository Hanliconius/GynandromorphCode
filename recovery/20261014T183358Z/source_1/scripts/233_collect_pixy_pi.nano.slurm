#!/usr/bin/env bash
#SBATCH --job-name=pixyPiSum
#SBATCH --partition=nano
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/233_collect_pixy_pi.%j.out
#SBATCH --error=logs/233_collect_pixy_pi.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
root="${project}/analysis/pixy_pi_17controls_allsites"
metadata="${root}/metadata"
pixy_root="${root}/pixy/by_chunk"
summary="${root}/summary"

mkdir -p "${summary}"

expected=$(
    awk 'END {print NR - 1}' \
        "${metadata}/chunks.tsv"
)

observed=$(
    find "${pixy_root}" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name 'task_*_pi.txt' \
        | wc -l
)

if [[ "${observed}" -ne "${expected}" ]]
then
    echo "ERROR: expected ${expected} pixy chunk outputs; found ${observed}" >&2
    exit 1
fi

combined="${summary}/pi_100kb_all_windows.tsv"
tmp="${combined}.new"

rm -f "${tmp}"

first=1

for task_id in $(seq 0 $((expected - 1)))
do
    printf -v tag '%04d' "${task_id}"

    file="${pixy_root}/task_${tag}/task_${tag}_pi.txt"

    if [[ ! -s "${file}" ]]
    then
        echo "ERROR: missing ${file}" >&2
        exit 1
    fi

    if [[ "${first}" -eq 1 ]]
    then
        cat "${file}" \
            > "${tmp}"
        first=0
    else
        tail -n +2 "${file}" \
            >> "${tmp}"
    fi
done

mv "${tmp}" "${combined}"

###############################################################################
# IMPORTANT:
# Genome-wide pi is NOT mean(window avg_pi).
# Reconstruct it as:
#
#       sum(count_diffs) / sum(count_comparisons)
#
###############################################################################

genomewide="${summary}/pi_genomewide_by_population.tsv"

awk \
    -F '\t' \
    -v OFS='\t' \
    '
    NR == 1 {
        for (i = 1; i <= NF; i++)
            col[$i] = i

        required[1] = "pop"
        required[2] = "count_diffs"
        required[3] = "count_comparisons"

        for (j = 1; j <= 3; j++) {
            if (!(required[j] in col)) {
                print "ERROR missing column " required[j] > "/dev/stderr"
                exit 1
            }
        }

        order[1] = "H_m_melpomene"
        order[2] = "H_m_bellula"
        order[3] = "H_m_amaryllis"
        order[4] = "H_m_malleti"
        order[5] = "H_m_rosina"
        order[6] = "H_timareta_thelxinoe"
        order[7] = "H_timareta_NG2008"

        next
    }

    {
        pop = $(col["pop"])

        diffs[pop] += $(col["count_diffs"])
        comparisons[pop] += $(col["count_comparisons"])
        windows[pop]++

        if ("no_sites" in col)
            sites[pop] += $(col["no_sites"])
    }

    END {
        print \
            "population", \
            "windows", \
            "sum_no_sites", \
            "sum_count_diffs", \
            "sum_count_comparisons", \
            "pi", \
            "pi_percent"

        for (i = 1; i <= 7; i++) {
            pop = order[i]

            if (comparisons[pop] <= 0) {
                pi = "NA"
                pct = "NA"
            } else {
                pi_value = diffs[pop] / comparisons[pop]
                pi = sprintf("%.8f", pi_value)
                pct = sprintf("%.4f", 100 * pi_value)
            }

            printf "%s\t%d\t%.0f\t%.6f\t%.6f\t%s\t%s\n", \
                pop, \
                windows[pop], \
                sites[pop], \
                diffs[pop], \
                comparisons[pop], \
                pi, \
                pct
        }
    }
    ' \
    "${combined}" \
    > "${genomewide}.new"

mv \
    "${genomewide}.new" \
    "${genomewide}"

echo
echo "===== GENOME-WIDE AUTOSOMAL PIXY PI ====="
column -t -s $'\t' "${genomewide}"

echo
echo "WINDOW TABLE:"
echo "${combined}"
echo
echo "GENOME-WIDE SUMMARY:"
echo "${genomewide}"
