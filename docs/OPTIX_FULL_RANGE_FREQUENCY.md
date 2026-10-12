# Full-range optix frequency track

Use `scripts/figures/figure3/submit_Fig3_optix_full_range_frequency.sh` on Pegasus.

The track spans Hmr18:1,060,000–1,240,000 and uses all eligible biallelic SNPs
with comparative BA >0.99 for Band, Dennis or Ray. It is not intersected with
the five legacy locus tables, except that the existing Band table supplies
the established comparative Band phenotype assignments. Multiple trait
assignments for one SNP are retained.

The BA sources and comparative manifest are the same as the existing BA099
workflow. Band uses analysis142 WGS_plus_CS / DP2_GQ20 / score_eligible;
Dennis and Ray use the Wallbank dominant-BA table. Calls are retrieved from
the existing original optix VCF and filtered by DP>=2 and GQ>=20. Both
comparative phenotype groups must be >=80% callable.

Population rows show the actual pattern-associated allele frequency among
called diploid individuals. Focal rows show reciprocal-mixture recalled
genotype dosage (0, 0.5 or 1), not the mixed-read allele fraction. Focal
posterior must be >=0.95. Neither focal genotype state nor focal callability
selects sites: unavailable focal values remain grey. There is no raw-GT
fallback and no phasing interpretation.

The figure has three aligned four-row tracks, with SNP ticks at true genomic
positions. Pale blue/gold/purple indicates frequency 0/0.5/1. Outputs include
PDF, PNG, candidate-site callability, full plotting values, and a source/cohort
audit. The audit records BA-selected variant keys absent from the original VCF;
its named interval ends at 1,239,421, slightly before the requested 1.24 Mb axis.

Source code is frozen into each dated output directory before submission.
This stage does not alter the running Yb rebuild.
