# Stacked heterozygosity and phenotype BA

Worker: `scripts/figures/figure3/stack_LR_heterozygosity_BA.py --submit` on Pegasus in `gynandrocalling`.

Output: `/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/LR_Ho_BA_stacked_v1`.

Reuses completed `LR_window_heterozygosity_v1` per-window statistics without new calling, counting, parameter fitting or changes to previous figures. All input files are SHA256 recorded in run.json.

Dennis reuses `Optix_Arnaud_review_nonPK_v2/merged/new_comparative_SNP_scores.tsv.gz` from Rosina recall `20261014T223341Z_f956c8bc`. All Dennis module scores are displayed, without a BA cutoff or filtering by focal genotype. Identity and BA values at the frozen selected SNPs must agree with `Dennis_BA90_intervals_sites.tsv`. The phenotype cohorts must remain n=142 positive and n=281 negative.

Yb(v) computes NEW EXPLORATORY scores from all SNPs in `ybv_REF_ALT_rebuild_v1/SNPs.tsv` and `individual_ALT_dosages.tsv.gz`, restricted to the individuals in the frozen `ybv_short_ALT_geography_v7/membership_geography.tsv`. Unknown phenotype entries are not scored; duplicate representatives or potential PK/focal entries cause an error. No new phenotype assignments or genotype quality rules are introduced; the existing quality-filtered dosage matrix is reused.

Balanced accuracy is one half sensitivity plus one half specificity, calculated among callable individuals. Each phenotype class requires at least 80% callability. The allele orientation with higher BA is used; ties remain unscored, matching the historical Dennis scoring rule. Two explicitly separate Yb versions are plotted:

* carrier: phenotype positive predicted by at least one ALT copy, or by at least one REF copy;
* homozygote: phenotype positive predicted by two ALT copies, or by two REF copies.

This explores two inheritance classifications; it does not establish dominance from these data. These are descriptive in-cohort scores, not cross-validated predictive performance, and retain the geographical/race exclusions of v7. Do not present this as a range-wide association test. BA uses comparative individuals, never the focal corrected calls.

The expanded Yb cached genotype window is only about 3 kb wide. Flanks with no expanded-cohort BA data are hatched, never assigned zero. This limits comparison with the upstream heterozygosity peak. Dennis and Yb use different cohorts and historical call sources, so BA should not be treated as identically ascertained across the two loci.

Outputs:

* `Dennis_Ho_BA_stacked.pdf`: four pages (wide/local views, 1/5 kb Ho windows).
* `Ybv_Ho_BA_stacked.pdf`: eight pages (wide/local views, 1/2 kb Ho windows, both classification rules).
* Exported Dennis and Yb score tables, Yb individual membership, input hashes, jobs.json, complete.json.

Each page aligns interval annotations, LR heterozygosity (posterior .95/.99), per-SNP BA, Right-minus-Left heterozygosity and discovery/accepted SNP counts. BA>=.9 points are highlighted; lower scores remain visible. Both curve axes are common between loci (Ho 0-.6; difference -.5 to .5). Actual data rendering must be checked after the nano job. Numerical checks cover known carrier/homozygote scores, REF orientation, class callability and orientation ties.
