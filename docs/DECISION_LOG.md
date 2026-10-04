# Analysis decision log

This file records methodological decisions that materially changed interpretation. It is not intended to reproduce every exploratory command.

## 1. Preserve paired Left/Right allele frequencies

Early plotting and validation made it clear that hard genotype labels throw away information needed to understand a left/right mosaic individual.

Decision:

Use paired read-count allele frequencies as the primary observed data for founder-model fitting. Genotype calls are a separate display layer.

## 2. Model reciprocal mixing explicitly

Observed Left and Right allele-frequency clouds are shifted toward one another.

Decision:

Fit separate reciprocal mixing fractions for the two directions rather than assuming zero mixing or forcing a single symmetric value.

Genome-wide calibrated estimates settled near 18% in each direction.

## 3. Separate raw model likelihood from calibrated founder inference

The raw genome-wide likelihood favored the constrained three-founder model.

Simulation calibration then showed that datasets generated under the strict three-founder model rarely reproduce the observed statistic.

Decision:

Do not report the raw likelihood winner as the biological conclusion. Use simulation-calibrated distributions, block resampling, and sibling-aware benchmarks.

## 4. Treat founder structure as a continuum of benchmarks

Later work compared a sequence of increasingly realistic four-founder benchmarks:

1. broad mixing;
2. mixing near the fitted region;
3. exact 20% mixing;
4. empirical full-sibling structure at exact 20%;
5. empirical full siblings with fitted asymmetric mixing.

Decision:

Interpret the observation as intermediate between idealized endpoints. Strict three-founder is inadequate, but an ideal unrelated-four-founder endpoint can overshoot. Sibling structure explains much of the intermediate signal.

## 5. Hold Z out of the founder-model fit

Sex-chromosome coverage differs strongly between sides.

Decision:

Use chromosomes 1-20 for founder-model and reciprocal-mixing inference. Keep chromosome 21/Z as an independent sex-ploidy diagnostic.

## 6. Require posterior >= 0.95 for focal reciprocal-mixture genotype display

The reciprocal-mixture caller provides posterior probabilities for focal Left/Right genotype pairs.

Decision:

Display only posterior >= 0.95 calls. Do not fall back to original GT/GQ or lower-confidence mixture calls.

## 7. Do not select Smudgeplot controls from mapped depth alone

Several apparently high-depth controls produced very poor heterozygous k-mer-pair structure.

Decision:

Judge Smudgeplot controls from FastK/hetmer behavior and the raw k-mer-pair density, not BAM mean depth alone.

## 8. ERR1051012 becomes the empirical diploid control

The older Krista config listed ERR1050907, ERR1050958, ERR1050994, and ERR1051012.

ERR1051012 was the clean useful control.

Decision:

Use ERR1051012 as the empirical diploid benchmark.

## 9. Coverage harmonic can relabel the same cloud

At ERR1051012 L=4, automatic 1n selection could pick a half-coverage harmonic and call the balanced signal AABB.

Forcing 1n about 8.01 on the same .smu data converted the dominant signal to AB without changing the underlying k-mer pairs.

Decision:

Never interpret an apparent ploidy class until the 1n harmonic has been checked.

## 10. One universal lower cutoff is inappropriate

The historical Krista workflow used L=12.

Lower-cutoff diagnostics showed that strongly imbalanced pairs can have a minor arm below 12 while still having substantial total coverage. Candidate 1:3-like components in pseudo-tetraploid and T2/T3 merged data were therefore vulnerable to censorship.

Decision:

Use dataset-appropriate lower cutoffs and explicitly document them.

Current preferred comparison:

- ERR1051012: L6
- pseudo-tetraploid ERR1307111 + ERR1307112: L6
- TN Left: L12
- TN Right: L12
- TN T2 merge: L6
- TN T3 merge: L6
- all four TN legs: L12

## 11. Keep raw-density plotting separate from class calling

Decision:

Maintain a custom vector renderer that plots the .smu geometry directly. Smudgeplot class labels remain diagnostic overlays rather than the sole representation.

## 12. Prefer vector figures and complete diagnostic bundles

Decision:

Use PDF/SVG for final figures. For Smudgeplot, keep both publication-style combined panels and complete native pages containing marginal histograms, 1n coverage, error fraction, and class proportions.

## 13. Preserve upstream lineage without obscuring downstream changes

The project grew from Krista Pipho's gynandromorph_2024 pipeline.

Decision:

Pin the exact upstream commit and document inherited settings separately from later downstream analyses. GynandromorphCode remains a distinct private reproducibility repository rather than pretending the later work was part of the original pipeline.
