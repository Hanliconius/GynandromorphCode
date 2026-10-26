# Final Figure 3 morph-frequency supplements and reviewer checks

## Correcting the old 80-column interpretation

The old `plot_optix_uniform_genoplots.py` at commit `f132ae9` has
`--max-sites=80`. Its `display_sites()` function sorts sites by genomic identity
and, when necessary, takes evenly spaced indices using linspace. Thus the
screenshot's 174 callable SNPs / 80 displayed means actual display thinning,
not the top 80 by BA, and not only 80 passing the threshold.

## New rendering

Run on Pegasus after pulling main:

```bash
micromamba run -n gynandrocalling python \
/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo/scripts/figures/figure3/build_Figure3_morph_supplements.py \
--submit --overwrite
```

Two nano array tasks independently render Optix/Dennis and Yb(v); a dependent
report job packages the outputs. Stable output directory:
`/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/Figure3_for_Arnaud/final_morph_supplements`.

Each locus gets a fixed BA >=0.93 version and a compact version. A threshold-count
table covers 0.90 through 1.00 in 0.01 steps. Compact uses the smallest tested
threshold retaining 1-80 sites. Every passing site, including ties, is retained.
If even BA=1 has more than 80 sites, all such ties remain and the failure to meet
the size target is reported. Pages show consecutive blocks of up to 80 columns;
pagination does not discard sites. Main-figure SNP selection/cohorts are unchanged.
BA=0.93 means 93%, not 0.93%.

The source is the completed reviewed Figure3_for_Arnaud package: membership,
all_SNP_BA and full-interval all_ALT_dosages caches. Source allele order/identities
and specimen counts are validated. All existing allele/filter rules remain.
Optix spans 1,060,000-1,240,000 on NC_140406.1, using the main Dennis BA score.
Yb spans 1,630,000-1,692,000 on NC_140403.1, retaining the main restricted Ybv
homozygote BA scores. This is not a fresh Band or Ray association scan.

Rows pool one reviewed biological specimen once, grouped by population plus
BDR morph for Optix, or population plus dorsal/ventral phenotype for Yb.
Yb is additionally separated by reviewed geography. BA-excluded geographical
exceptions remain in the Yb display when qualified broad calls exist; the
selection cohort does not expand to include them. Global exclusions (including
numata/pardalinus) still apply, and PKs do not enter comparative Yb pools.
Only available data are displayed; missing phenotypes are not inferred.

Full-range Yb calls come from the documented Ybv_BA_broad caches; existing focus
calls are used at those exact alleles where available. Main BA cohort calls are
checked for agreement and then substituted verbatim for those specimens. Source
conflicts are exported. There is no downloading, alignment, or variant calling.

The pale-blue / gold / purple palette encodes mean individual ALT dosage / 2.
This is an allele-frequency gradient: its middle colour does NOT establish that
all specimens are heterozygous. Every row-marker combination exports n_total,
n_called, ALT frequency and the three genotype proportions. Main BA selection
membership and broader display membership are separate exported quantities.
These panels display comparative populations; corrected focal genotypes remain
in the existing main/focus panels and are not extrapolated to uncalled sites.

## Reviewer-focused assessment

These are plausible substantive questions, not a prediction of Jim Mallet's
personal review. Prioritize them over another cosmetic threshold adjustment.

1. **Population structure and ascertainment.** A high BA may distinguish ancestry,
   race or geography rather than a universally causal phenotype allele. Keep the
   Yb eastern exceptions visible, state exactly which samples selected the SNPs,
   and frame performance outside that cohort as generalization, not reselection.
   If one additional validation is undertaken, leave-population-out or
   leave-lineage-out performance would be more informative than random individual
   splits. Neither high in-sample BA nor a larger BA cutoff solves confounding.

2. **Selection versus validation.** A frequency plot made from the SNP-discovery
   individuals is a descriptive illustration, not independent validation. Raising
   BA to make the figure compact is an openly reported display choice, not a
   stronger significance test. Record every passing site and threshold count.
   Hundreds of linked SNPs/dots are not hundreds of independent confirmations.

3. **Phenotype genetics and nomenclature.** Explicitly separate dorsal and ventral
   Yb scores; 1 means the scored pattern is present, not a dominant functional
   allele. Audit unknowns, mixed populations and within-race exceptions. Define
   East/West from documented locality/grouping rules independently of the SNP
   pattern, and distinguish geographic position from phylogeographic ancestry.
   Maroja et al. found cryptic Yb variation in true-breeding H. m. melpomene and
   no complementation between rosina and amaryllis; neither observation warrants
   assuming all yellow-band populations share every linked diagnostic SNP.

4. **Marker association versus mechanism/history.** Avoid calling an associated
   SNP causal or an enhancer without functional evidence. Shared alleles do not
   alone establish the direction or timing of introgression. Existing optix
   modularity/introgression literature and the cortex-region studies provide
   context, not independent proof for each marker in this specimen. Explain what
   the new bilateral specimen contributes beyond the established association maps.

5. **Technical comparability and pooling.** State legacy/new call sources,
   reference harmonization, callability filters and missingness. Check whether
   source and phenotype/geography are confounded. Keep individual genotypes and
   per-marker denominators available behind the compact population display.
   Consensus letters and pooled intermediate AF are not phased haplotypes or
   evidence that every butterfly carries the same genotype.

6. **The specimen-level biological claim.** Retain raw AD and independent depth
   support alongside mixture-corrected calls. Distinguish DNA/read mixing,
   biological mosaicism and assumed founder/ploidy models. The male-on-both-sides
   NDJ butterfly should be described as a bilateral pattern/aneuploid mosaic
   unless sexual mosaicism is independently established. Do not treat two tissues
   from one butterfly as two independent specimens.

Primary literature to address explicitly:

- Wallbank et al. 2016, Evolutionary novelty in a butterfly wing pattern through
  enhancer shuffling: https://doi.org/10.1371/journal.pbio.1002353
- Maroja et al. 2012, Partial complementarity of the mimetic yellow bar phenotype:
  https://doi.org/10.1371/journal.pone.0048627
- Enciso-Romero et al. 2017, Evolution of novel mimicry rings facilitated by
  adaptive introgression: https://doi.org/10.1111/mec.14277

This rendering workflow addresses transparent display and cohort accounting;
it does not claim that the additional cross-population validation or causal
experiments proposed above have been performed.
