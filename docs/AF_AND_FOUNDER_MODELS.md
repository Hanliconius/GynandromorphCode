# Allele-frequency and founder-model analysis

## Goal

Use the joint Left/Right allele-frequency distribution to test constrained founder-state models while fitting reciprocal mixing between sides.

The analysis is based on raw merged Left/Right read counts rather than treating original focal GT/GQ calls as truth.

## Core observed quantity

For each autosomal site:

(AF_Left, AF_Right)

The focal fit uses read counts from the merged Left and Right samples.

## Reciprocal mixing model

Observed Left  = (1-cL) * true Left + cL * true Right

Observed Right = cR * true Left + (1-cR) * true Right

The fitted contamination fractions are estimated separately for each direction.

A later genome-wide calibrated summary gave approximately:

- right material in Left: 17.7%
- left material in Right: 18.3%

Earlier model-conditioned fits were close to this range but differed somewhat between three- and four-founder parameterizations.

## Founder models

### Three-founder model

The two sides share one founder chromosome.

Conceptually:

Left = A + B

Right = A + C

This produces seven allowed paired latent genotype states and excludes opposing homozygous states that require four independent founder chromosomes.

### Four-founder model

The two sides have independent founder pairs.

Conceptually:

Left = A + B

Right = C + D

This permits the full nine paired genotype-state combinations.

## Raw fit versus calibrated test

The raw genome-wide model likelihood favored the three-founder model.

That raw preference is not, by itself, the final biological inference.

Simulation calibration asked whether datasets generated under each founder model reproduce the observed test statistic after realistic depth, sequencing error, bias, overdispersion, and genomic block structure are included.

An early calibrated analysis gave an observed nested LR of about 118.63 and rejected the strict three-founder simulation at an empirical Monte Carlo P of about 0.0164. A high-call sensitivity analysis gave about 0.0244.

This P value is a calibrated tail probability for the test. It is not a posterior probability that four founders are true.

## Sequential distribution ladder

Later work deliberately narrowed the four-founder comparison:

1. broad reciprocal mixing over 0-45%;
2. mixing focused near the observed region;
3. exact 20% reciprocal mixing;
4. empirical full-sibling founder structure at exact 20%;
5. empirical full siblings with fitted asymmetric mixing.

The observed genome-wide statistic lies between the strict three-founder and ideal unrelated-four-founder endpoints.

Broad four-founder simulations can contain the observation, whereas focused exact-20% unrelated four-founder simulations tend to overshoot. Full-sibling structure explains much of the intermediate signal but still does not perfectly reproduce the observation.

This is the current interpretive framing.

## Genome-wide scope

The genome-wide calibrated analysis used a large autosomal fit set and block-based resampling/validation. Chromosomes 1-20 determine founder-model and reciprocal-mixing inference.

Chromosome 21/Z is held out as an independent sex-ploidy diagnostic.

Regional heterogeneity is important; some genomic regions contribute much more strongly than others.

## Important scripts to import from Pegasus

Historical scripts include:

- scripts/158_audit_three_vs_four_founder_inputs.py
- scripts/159_fit_autosomal_three_vs_four_founder.py
- scripts/160_simulate_calibrate_autosomal_founder_models.py

Later sequential-calibration / narrative scripts should also be imported and renamed descriptively.

## Output types worth retaining

The analysis generated or should retain:

- observed model fits;
- site-level evidence;
- 10-kb block evidence;
- optimizer multi-start diagnostics;
- simulation recovery tables;
- block-bootstrap uncertainty;
- checksums / input-audit tables;
- narrative PDFs that explain the model sequence.

## Guardrail

Do not collapse the result to a one-line ploidy call. The founder-model analysis is a calibrated comparison of generative models under reciprocal mixing, not a direct chromosome count.
