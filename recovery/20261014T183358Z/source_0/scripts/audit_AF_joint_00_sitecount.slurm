#!/bin/bash
#SBATCH --job-name=AFjoint00
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/AF_joint_00_sitecount.%j.out
#SBATCH --error=logs/AF_joint_00_sitecount.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
bcftools=/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/bcftools

bcf="${project}/analysis/genomewide_mixing_input/variant_calls/joint_24/genomewide/all_24_samples.joint_at_control_sites.chromosomes_1_to_21.bcf"

outdir="${project}/analysis/figures/pairwise_AF_cli/site_count_audit"
outfile="${outdir}/TN_joint_original_GT_site_counts.tsv"

mkdir -p "${outdir}" "${project}/logs"

"${bcftools}" query \
    -s TN_Left_merged,TN_Right_merged \
    -f '%CHROM\t%POS\t%INFO/AC[\t%GT\t%AD]\n' \
    "${bcf}" \
| awk -F '\t' '
function ad_ok(s, a, n) {
    n = split(s, a, ",")
    if (n != 2) return 0
    if (a[1] == "." || a[2] == ".") return 0
    if ((a[1] + 0) < 0 || (a[2] + 0) < 0) return 0
    return ((a[1] + 0) + (a[2] + 0) >= 10)
}

function gtcat(g) {
    gsub(/\|/, "/", g)

    if (g == "1/0")
        g = "0/1"

    if (g == "0/0" || g == "0/1" || g == "1/1")
        return g

    return "missing_or_other"
}

{
    records_read++

    # Fields:
    # 1 CHROM
    # 2 POS
    # 3 INFO/AC
    # 4 Left GT
    # 5 Left AD
    # 6 Right GT
    # 7 Right AD

    if (!(ad_ok($5) && ad_ok($7))) {
        excluded_depth++
        next
    }

    retained++

    left_gt  = gtcat($4)
    right_gt = gtcat($6)

    pair[left_gt "|" right_gt]++

    if (left_gt == "0/0" && right_gt == "0/0") {
        joint00++

        split($5, lad, ",")
        split($7, rad, ",")

        if ((lad[2] + 0) == 0 && (rad[2] + 0) == 0)
            joint00_AF00++
        else
            joint00_with_ALT_reads++

        # Because both merged TN genotypes are 0/0 here,
        # AC > 0 means some other sample in the joint BCF
        # carries a called ALT allele.
        if ($3 != "." && ($3 + 0) > 0)
            joint00_other_sample_ALT_GT++
        else
            joint00_joint_AC_zero++
    }
}

END {
    print "metric\tcount"

    print "records_read\t" records_read
    print "retained_AD10_both\t" retained
    print "excluded_below_AD10_either_side\t" excluded_depth

    print "joint_original_GT_0/0|0/0\t" joint00
    print "joint_0/0|0/0_and_AF_exactly_0|0\t" joint00_AF00
    print "joint_0/0|0/0_with_ALT_reads_in_TN\t" joint00_with_ALT_reads
    print "joint_0/0|0/0_with_INFO_AC_gt_0\t" joint00_other_sample_ALT_GT
    print "joint_0/0|0/0_with_INFO_AC_eq_0\t" joint00_joint_AC_zero

    for (key in pair)
        print "GT_pair_" key "\t" pair[key]
}
' > "${outfile}.tmp"

mv "${outfile}.tmp" "${outfile}"

records=$(
    awk -F '\t' '$1=="records_read" {print $2}' "${outfile}"
)

retained=$(
    awk -F '\t' '$1=="retained_AD10_both" {print $2}' "${outfile}"
)

if [[ "${records}" -ne 35195121 ]]; then
    echo "ERROR: expected 35,195,121 records; observed ${records}" >&2
    exit 1
fi

if [[ "${retained}" -ne 31126758 ]]; then
    echo "ERROR: expected 31,126,758 retained sites; observed ${retained}" >&2
    exit 1
fi

echo
echo "===== AF SITE-COUNT AUDIT ====="
column -t -s $'\t' "${outfile}"

echo
echo "WROTE:"
echo "${outfile}"
