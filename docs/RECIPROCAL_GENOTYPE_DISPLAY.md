# Reciprocal-mixture genotype display: verified provenance

## Scope of the existing calls

The supplied Pegasus caller `scripts/134_fit_reciprocal_L_R_mixture_genotypes.py` (version 1.0.0) reads `analysis/original_502/checkpoints/optix_fullVCF_focal_and_WGS_control_AF.long.tsv.gz`. Its expected SHA256 is `2129355063e80b19588d605a7402a316782b38da847f0736a880fa8997a98db9`. Its outputs under `analysis/original_502/reciprocal_L_R_mixture/` are optix-specific. Genome-wide fitting of mixing parameters does not establish genome-wide availability of corrected genotype calls.

The supplied `scripts/92c_build_modeled_LR_frozen_site_table.py` joins frozen comparative SNPs to `optix_reciprocal_L_R_genotype_pairs.posterior95.tsv`. It explicitly requires 117 unique frozen sites and excludes T2/T3-specific output columns. The resulting `checkpoints/focal_calls_at_frozen_sites_biological.modeled_LR_posterior95.wide.tsv` is therefore not a genome-wide corrected-call dataset. Its focal calls use merged TN_Left/TN_Right counts.

## Caller model verified from the supplied excerpt

Each side has three possible diploid ALT dosages: 0, 0.5, 1, corresponding to 0/0, 0/1, 1/1. The caller enumerates nine paired Left/Right states. For true dosages gL and gR, reciprocal contamination gives observed biological fractions:

- Left: (1-cL) gL + cL gR.
- Right: cR gL + (1-cR) gR.

The supplied code defines mapping bias as a log-odds shift, followed by symmetric sequencing error: p = error + (1-2 error) expit(logit(biological_fraction) + bias_log_odds). Fractions are clipped to 1e-6 through 1-1e-6 before this transformation. Read-count likelihoods use a beta-binomial distribution with alpha=p*concentration and beta=(1-p)*concentration; concentration is floored at 2 in the likelihood helper.

Declared thresholds are minimum calling depth 2, minimum fitting depth 10, control minimum GQ 20 and pair posterior cutoff 0.95. The excerpt does not establish all fitting filters, priors, optimizer settings or the final fitted values; these must be recovered from the complete caller and its summary/manifest, not reconstructed from rounded narrative values.

## Verified display and orientation rules

The frozen-site builder joins on chromosome, position, REF and ALT, checks uniqueness, and preserves the frozen site count. Both genotype labels must be valid, the posterior95 flag must be true, and the pair posterior must be >=0.95. Otherwise both displayed genotypes are missing; there is no original-genotype fallback.

When the phenotype-associated allele is ALT, displayed dosage is ALT copy count / 2. When it is REF, dosage is 1 minus ALT dosage. The builder encodes 0 as AA, 0.5 as PA and 1 as PP. These are unphased marker genotypes, not demonstrated causal haplotypes. Its `left_presence_allele_fraction` and `right_presence_allele_fraction` columns retain continuous observed read fractions, even for confident calls: plots of corrected genotype must use the discrete `*_presence_dosage` columns instead.

The builder's `left_GT_raw` and `right_GT_raw` columns actually contain accepted model genotypes, despite their legacy names. They must not be interpreted as the original uncorrected VCF GT.

## Applying the workflow at Yb(d) and Yb(v)

First recover the exact fitted model used by the desired current analysis, including contamination in each direction, error, mapping bias, concentration and paired-state priors. Different genome-wide AF runs may use different fits. Never substitute approximate documentation values or the exploratory symmetric 20% correction.

Reuse existing calls only when reference assembly, coordinates, REF/ALT, sample aggregation and model provenance match. Otherwise run the established paired-state inference locally at the Yb SNPs using the existing BAM/count data. New Hmel2.5 and HmRos alignments require separate coordinate/allele checks; chromosome aliases alone do not authorize coordinate reuse.

Keep SNP selection independent of focal Left/Right. Retain measured control population frequencies, but show focal posterior95 genotype dosage only as 0, 0.5 or 1, with unavailable calls grey. Match the optix palette: 0 #D7EDF5; 0.5 #ECCB63; 1 #662A83; missing #C8C8C8. Match merged-side inference when reproducing the existing optix workflow; do not present separately inferred T2/T3 genotypes as though they were that same existing output.

The initial `yb_interval_genoplot.py` output is an exploratory read-frequency diagnostic, not corrected genotype calls. Its linear symmetric-20% frequency correction cannot serve as the corrected genotype display. Yb intervals in that script are legacy annotations; the Luca primer span is a separate anchor and not a proven exact CRE boundary. Rosina/melpomene SNP differences are race-associated candidates, especially at Yb(v), until phenotype association and independent controls are assessed.
