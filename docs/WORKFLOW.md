# Stepwise analysis workflow

This document records the conceptual sequence of the gynandromorph analysis. It is intended to be more stable and readable than the historical numbered Pegasus scripts.

## 1. Four-library focal design

The focal individual has four sequencing libraries:

- TN Left T2 = SRR39053777 / TN_Left_T2
- TN Right T2 = SRR39053776 / TN_Right_T2
- TN Left T3 = SRR39053775 / TN_Left_T3
- TN Right T3 = SRR39053774 / TN_Right_T3

The same biological individual is therefore observable by side and by tissue replicate.

Useful arrangements:

- Left = Left T2 + Left T3
- Right = Right T2 + Right T3
- T2 merge = Left T2 + Right T2
- T3 merge = Left T3 + Right T3
- all four legs = all four focal libraries

Left/Right preserve the putative biological contrast. T2/T3 cross-side merges deliberately combine the two sides and therefore act as useful mixture diagnostics.

## 2. Site-level allele-frequency representation

The fundamental observation at an informative site is the paired read-count allele frequency:

(AF_Left, AF_Right)

This retains information that would be lost by immediately reducing the data to hard diploid genotype calls.

Important plotting conventions developed during the analysis:

- grouped rather than stacked marginal displays;
- reference / heterozygous / alternate strata where useful;
- empirical control trends shown as smooth reference curves where appropriate;
- model centroids should not obscure the observed density;
- 50% predictive regions should be thin vector rings rather than opaque polygons.

For founder-model cartoons:

- there is no spread at (0,0) or (1,1);
- when one dimension is fixed at zero, sampling spread is one-dimensional;
- two-dimensional spread is shown only when both dimensions can vary.

## 3. Reciprocal mixing

Left and Right libraries can contain reciprocal read contamination or mixing. This is modeled explicitly.

A simple observation model is:

Observed Left  = (1-cL) * true Left  + cL * true Right

Observed Right = cR * true Left + (1-cR) * true Right

Genome-wide fitted reciprocal mixing is approximately 18% in each direction, with a later genome-wide summary around 17.7% right material in Left and 18.3% left material in Right.

These are fitted nuisance parameters, not fixed biological constants, and the code should retain the ability to refit them.

## 4. Founder-state inference

The key autosomal question is whether the observed Left/Right AF distribution is compatible with a constrained three-founder model or requires a four-founder-like signal.

The final analysis separates raw model fit from simulation calibration.

### Three-founder model

One founder chromosome is shared between sides. This constrains which paired Left/Right genotype states are possible.

### Four-founder model

Each side has its own diploid founder pair. This expands the possible paired genotype states.

### Why the calibration matters

The raw likelihood comparison and the calibrated simulation test do not tell exactly the same story.

The raw fit favored the three-founder model, but simulation calibration showed that the strict three-founder generative model did not reproduce the observed statistic well. Later genome-wide and sibling-calibrated analyses placed the observation between the idealized three-founder and unrelated four-founder endpoints.

The current interpretation is therefore not simply "three founders" or "four founders". The important result is that the strict three-founder model is inadequate, while idealized unrelated-four-founder simulations can overshoot the observed signal. Empirical sibling structure explains much of the intermediate signal but not all of it.

See AF_AND_FOUNDER_MODELS.md for details.

## 5. High-confidence reciprocal-mixture genotype display

A separate reciprocal-mixture genotype-calling step is used for focal display.

Rule:

Only reciprocal-mixture Left/Right genotype pairs with posterior >= 0.95 are displayed. There is no fallback to lower-confidence or original calls.

Historical input:

analysis/original_502/reciprocal_L_R_mixture/
    optix_reciprocal_L_R_genotype_pairs.posterior95.tsv.gz

This display layer must remain conceptually separate from the read-count model fitting.

## 6. Chromosome-scale coverage

Relative depth provides an independent diagnostic, especially for sex-chromosome interpretation.

Representative Z-relative coverage:

| library | relative Z coverage |
|---|---:|
| Left T2 | 0.592628 |
| Left T3 | 0.618390 |
| Right T2 | 0.881962 |
| Right T3 | 0.938852 |

The genome-wide founder analysis uses chromosomes 1-20 for founder/mixing inference and treats chromosome 21/Z as a held-out sex-ploidy diagnostic.

## 7. Smudgeplot / k-mer-pair analysis

Smudgeplot provides an alignment-free complement to the AF and coverage analyses.

The major methodological lesson is that both the lower-count cutoff and the assumed 1n coverage can materially change the apparent smudge classes.

See SMUDGEPLOT.md.

## 8. Figure production

Primary figure products should be vector PDF/SVG whenever practical. PNG is useful for quick inspection but is not the archival format.

## 9. Keep evidence layers separate

The repository preserves the distinction among:

1. direct observations: read counts, allele frequencies, depth, k-mer-pair density;
2. fitted nuisance processes: reciprocal mixing;
3. founder-state generative models;
4. high-confidence genotype labels used only for display;
5. biological interpretation.

This prevents fitted or corrected quantities from later being mistaken for independent evidence.
