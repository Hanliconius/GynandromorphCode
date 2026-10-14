#!/bin/bash
#SBATCH --job-name=taxonPi
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=03:00:00
#SBATCH --output=logs/autosomal_pi_by_taxon.%j.out
#SBATCH --error=logs/autosomal_pi_by_taxon.%j.err

set -euo pipefail

bcftools=/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/bcftools

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

bcf="${project}/analysis/genomewide_mixing_input/variant_calls/joint_24/genomewide/all_24_samples.joint_at_control_sites.chromosomes_1_to_21.bcf"

manifest="${project}/analysis/genomewide_mixing_input/metadata/chromosome_calling_manifest.tsv"

outdir="${project}/analysis/genomewide_mixing_input/variant_calls/joint_24/audits"
outfile="${outdir}/autosomal_pi_by_taxon.tsv"

mkdir -p "$outdir" logs

autosome_bp=$(
    awk -F $'\t' '
        NR > 1 && $3 != "NC_140409.1" {
            sum += $4
        }
        END {
            printf "%.0f\n", sum
        }
    ' "$manifest"
)

echo "AUTOSOMAL_BP=${autosome_bp}"
echo "EXPECTED_FROM_TOTAL_MINUS_Z=$((290171155 - 14968240))"

if [[ "$autosome_bp" -ne 275202915 ]]; then
    echo "ERROR: unexpected autosomal genome length: $autosome_bp" >&2
    exit 1
fi

samples="ERR3790122,ERR3790123,ERR3790118,ERR3790124,ERR2298257,ERR2298260,ERR2298201,ERR2298216,ERR1307113,ERR1307112,ERR1307111,ERR2298208,ERR2298221,ERR2298265,ERR2298263,ERR2298235,ERR2298236"

"$bcftools" query \
    -s "$samples" \
    -f '%CHROM\t%POS[\t%GT]\n' \
    "$bcf" \
| awk -F '\t' -v genome_bp="$autosome_bp" '
BEGIN {
    name[1] = "H_m_melpomene"
    name[2] = "H_m_bellula"
    name[3] = "H_m_amaryllis"
    name[4] = "H_m_malleti"
    name[5] = "H_m_rosina"
    name[6] = "H_timareta_thelxinoe"
    name[7] = "H_timareta_NG2008"

    individuals[1] = 2
    individuals[2] = 4
    individuals[3] = 2
    individuals[4] = 3
    individuals[5] = 2
    individuals[6] = 2
    individuals[7] = 2

    for (i=3;  i<=6;  i++) group[i]=2
    for (i=7;  i<=8;  i++) group[i]=3
    for (i=9;  i<=10; i++) group[i]=1
    for (i=11; i<=13; i++) group[i]=4
    for (i=14; i<=15; i++) group[i]=5
    for (i=16; i<=17; i++) group[i]=6
    for (i=18; i<=19; i++) group[i]=7
}

function add_gt(gt, g, n, a) {
    gsub(/\|/, "/", gt)

    n = split(gt, a, "/")

    if (n != 2)
        return

    if ((a[1] != "0" && a[1] != "1") || (a[2] != "0" && a[2] != "1"))
        return

    called_alleles[g] += 2
    alt_alleles[g] += (a[1] + 0) + (a[2] + 0)
}

{
    if ($1 == "NC_140409.1")
        next

    autosomal_records++

    delete called_alleles
    delete alt_alleles

    for (i=3; i<=NF; i++)
        add_gt($i, group[i])

    for (g=1; g<=7; g++) {
        n = called_alleles[g]
        k = alt_alleles[g]

        if (n < 2)
            continue

        callable_snp_records[g]++

        site_pi = 2 * k * (n - k) / (n * (n - 1))
        pi_sum[g] += site_pi

        if (k > 0 && k < n)
            segregating_sites[g]++
    }
}

END {
    print "group\tn_individuals\tsegregating_sites\tcallable_discovery_records\tpi_sum\tpi_per_bp\tpi_percent"

    for (g=1; g<=7; g++) {
        pi = pi_sum[g] / genome_bp

        printf "%s\t%d\t%d\t%d\t%.6f\t%.8f\t%.4f\n", \
            name[g], \
            individuals[g], \
            segregating_sites[g], \
            callable_snp_records[g], \
            pi_sum[g], \
            pi, \
            100 * pi
    }

    print ""
    print "autosomal_VCF_records\t" autosomal_records
    print "autosomal_genome_bp\t" genome_bp
}
' > "${outfile}.tmp"

mv "${outfile}.tmp" "$outfile"

echo
echo "===== AUTOSOMAL PI BY TAXON ====="
column -t -s $'\t' "$outfile"

echo
echo "Definition: pi_i = 2*k*(n-k) / [n*(n-1)]"
echo "Genome-wide pi = sum(pi_i) / autosomal genome bp"
echo
echo "WROTE: $outfile"
