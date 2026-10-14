#!/bin/bash
#SBATCH --job-name=taxonVar
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=03:00:00
#SBATCH --output=logs/variant_sites_by_taxon.%j.out
#SBATCH --error=logs/variant_sites_by_taxon.%j.err

set -euo pipefail

bcftools=/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/bcftools

bcf=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/genomewide_mixing_input/variant_calls/joint_24/genomewide/all_24_samples.joint_at_control_sites.chromosomes_1_to_21.bcf

outdir=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/genomewide_mixing_input/variant_calls/joint_24/audits
outfile="${outdir}/variant_sites_by_subspecies_and_gynandromorph.tsv"

mkdir -p "$outdir" logs

# Order is important: the awk block below maps these columns to groups.
samples="ERR3790122,ERR3790123,ERR3790118,ERR3790124,ERR2298257,ERR2298260,ERR2298201,ERR2298216,ERR1307113,ERR1307112,ERR1307111,ERR2298208,ERR2298221,ERR2298265,ERR2298263,ERR2298235,ERR2298236,TN_Left_merged,TN_Right_merged"

"$bcftools" query \
    -s "$samples" \
    -f '%CHROM\t%POS[\t%GT]\n' \
    "$bcf" \
| awk -F '\t' '
BEGIN {
    name[1] = "H_m_melpomene"
    name[2] = "H_m_bellula"
    name[3] = "H_m_amaryllis"
    name[4] = "H_m_malleti"
    name[5] = "H_m_rosina"
    name[6] = "H_timareta_thelxinoe"
    name[7] = "H_timareta_NG2008"
    name[8] = "gynandromorph_merged_LR"

    n_columns[1] = 2
    n_columns[2] = 4
    n_columns[3] = 2
    n_columns[4] = 3
    n_columns[5] = 2
    n_columns[6] = 2
    n_columns[7] = 2
    n_columns[8] = 2

    # bcftools query columns:
    # 1 CHROM
    # 2 POS
    #
    # 3-6   bellula
    # 7-8   amaryllis
    # 9-10  melpomene
    # 11-13 malleti
    # 14-15 rosina
    # 16-17 timareta thelxinoe
    # 18-19 timareta NG-2008
    # 20-21 gynandromorph merged L/R

    for (i=3;  i<=6;  i++) group[i]=2
    for (i=7;  i<=8;  i++) group[i]=3
    for (i=9;  i<=10; i++) group[i]=1
    for (i=11; i<=13; i++) group[i]=4
    for (i=14; i<=15; i++) group[i]=5
    for (i=16; i<=17; i++) group[i]=6
    for (i=18; i<=19; i++) group[i]=7
    for (i=20; i<=21; i++) group[i]=8
}

function inspect_gt(gt, g,    n,a,j,valid,allele) {
    gsub(/\|/, "/", gt)

    n = split(gt, a, "/")

    if (n < 1 || n > 2)
        return

    valid = 1

    for (j=1; j<=n; j++) {
        if (a[j] != "0" && a[j] != "1") {
            valid = 0
            break
        }
    }

    if (!valid)
        return

    called[g]++

    for (j=1; j<=n; j++) {
        allele = a[j]

        if (allele == "0")
            has_ref[g] = 1
        else if (allele == "1")
            has_alt[g] = 1
    }
}

{
    total++

    delete called
    delete has_ref
    delete has_alt

    for (i=3; i<=NF; i++)
        inspect_gt($i, group[i])

    for (g=1; g<=8; g++) {
        if (called[g] > 0)
            callable[g]++

        if (has_ref[g] && has_alt[g])
            variant[g]++
    }
}

END {
    expected = 35195121

    if (total != expected) {
        printf "ERROR: expected %d records but saw %d\n", \
            expected, total > "/dev/stderr"
        exit 1
    }

    print "group\tn_sample_columns\tvariant_sites\tpct_all_discovery_sites\tcallable_sites\tpct_variant_among_callable"

    for (g=1; g<=8; g++) {
        pct_all = 100 * variant[g] / total
        pct_callable = (callable[g] > 0) \
            ? 100 * variant[g] / callable[g] \
            : 0

        printf "%s\t%d\t%d\t%.4f\t%d\t%.4f\n", \
            name[g], \
            n_columns[g], \
            variant[g], \
            pct_all, \
            callable[g], \
            pct_callable
    }
}
' > "${outfile}.tmp"

mv "${outfile}.tmp" "$outfile"

echo
echo "===== VARIANT SITES BY TAXON ====="
column -t -s $'\t' "$outfile"

echo
echo "Denominator for pct_all_discovery_sites = 35,195,121 control-discovered SNP records."
echo "Variant = both REF and ALT represented among called GT alleles in that group."
echo "Gynandromorph uses TN_Left_merged + TN_Right_merged jointly."
echo
echo "WROTE: $outfile"
