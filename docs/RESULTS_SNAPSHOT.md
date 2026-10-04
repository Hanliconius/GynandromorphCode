# Working results snapshot

This is a compact record of the current state of the analysis. It is not a substitute for final statistical reporting.

## Reciprocal mixing and founder models

- genome-wide reciprocal mixing is approximately 18% in each direction;
- later summary estimates are about 17.7% right material in Left and 18.3% left material in Right;
- the raw model likelihood favored the constrained three-founder model;
- simulation calibration rejects the strict three-founder generative model as an adequate explanation of the observed statistic;
- the observation lies between idealized three-founder and unrelated-four-founder endpoints;
- full-sibling structure explains much, but not all, of the intermediate signal;
- regional heterogeneity matters and should remain visible in the final interpretation.

## Relative Z coverage

| library | relative Z coverage |
|---|---:|
| Left T2 | 0.592628 |
| Left T3 | 0.618390 |
| Right T2 | 0.881962 |
| Right T3 | 0.938852 |

Chromosomes 1-20 are used for founder/mixing inference; chromosome 21/Z is retained as a held-out sex-ploidy diagnostic.

## Smudgeplot control calibration

ERR1051012 is the preferred empirical diploid control.

The diagnostic L=4 run showed why coverage harmonics matter: automatic coverage could label the balanced cloud AABB, while manual 1n=8.01 labelled the same underlying cloud AB.

A later cutoff sweep preferred L=6 for the final diploid control.

The final lower-cutoff comparison uses L=6 for the pseudo-tetraploid and cross-side T2/T3 merges, while preserving L=12 for the well-behaved TN Left/Right analyses and the all-four combined analysis.

## Original L=12 TN reference

| dataset | total pairs | error fraction | inferred 1n | major called classes |
|---|---:|---:|---:|---|
| TN Left | 53,496,929 | 0.090 | 20 | AB 44,600,980; AABB 2,356,462 |
| TN Right | 22,222,743 | 0.072 | 20 | AB 16,529,698; AABB 2,825,354 |
| T2 merge | 52,835,318 | 0.551 | 20 | AB 20,312,886; AABB 1,702,456 |
| T3 merge | 56,207,120 | 0.581 | 22 | AB 21,858,988 |

## Interpretation guardrails

- smudge labels depend on assumed/inferred 1n coverage;
- lower-cutoff choice changes which imbalanced pairs enter the analysis;
- manual coverage changes labels but cannot recover k-mer pairs removed earlier by the cutoff;
- the pseudo-tetraploid is a four-copy benchmark made by merging diploid individuals, not an organismal-ploidy control;
- AF, depth, and k-mer-pair analyses are complementary and should not be collapsed into one metric.
