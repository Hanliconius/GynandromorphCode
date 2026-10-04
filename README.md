# GynandromorphCode

Reproducibility repository for genomic analyses of a *Heliconius* gynandromorph, with emphasis on left/right tissue differences, reciprocal read mixing, founder-state models, chromosome-scale coverage, and k-mer-pair (Smudgeplot) diagnostics.

This repository preserves both the final analysis path and the reasoning that led to it. Exploratory branches and failed parameterizations are documented when they changed how the final analysis was interpreted.

## Biological question

The focal individual was sequenced from left and right tissue, with T2 and T3 libraries available for each side. The central questions are:

1. Do the left and right sides differ genetically?
2. Can those differences be explained by reciprocal mixing of reads between sides?
3. Are the autosomal allele-frequency patterns compatible with a constrained three-founder model, or do they require additional four-founder-like structure?
4. Do chromosome-scale depth and k-mer-pair patterns support the same interpretation?
5. Which apparent Smudgeplot classes are biological, and which are artifacts of cutoff or coverage-harmonic choice?

## Focal samples

| role | sample |
|---|---|
| Left T2 | SRR39053777 / TN_Left_T2 |
| Right T2 | SRR39053776 / TN_Right_T2 |
| Left T3 | SRR39053775 / TN_Left_T3 |
| Right T3 | SRR39053774 / TN_Right_T3 |
| diploid Smudgeplot control | ERR1051012 |
| pseudo-tetraploid benchmark | ERR1307111 + ERR1307112 |

## Analysis map

The current workflow is:

1. Prepare and QC the four focal libraries.
2. Measure paired Left/Right allele frequencies from read counts.
3. Estimate reciprocal mixing between the two sides.
4. Compare three-founder and four-founder generative models with simulation calibration rather than relying on raw likelihood alone.
5. Use reciprocal-mixture genotype calls for focal display only when posterior probability is at least 0.95; there is no fallback to lower-confidence calls.
6. Inspect chromosome-scale relative coverage, including the Z chromosome as a held-out diagnostic.
7. Run Smudgeplot arrangements for Left, Right, T2-merged, T3-merged, all four legs, a diploid control, and a pseudo-tetraploid benchmark.
8. Calibrate Smudgeplot lower cutoff and 1n coverage explicitly rather than assuming one cutoff or trusting an automatically chosen harmonic.
9. Generate vector PDF/SVG figures and retain diagnostic reports alongside publication-style summaries.

See:

- [docs/WORKFLOW.md](docs/WORKFLOW.md) for the stepwise analysis;
- [docs/AF_AND_FOUNDER_MODELS.md](docs/AF_AND_FOUNDER_MODELS.md) for the founder-model logic and calibration;
- [docs/SMUDGEPLOT.md](docs/SMUDGEPLOT.md) for the k-mer workflow and cutoff calibration;
- [docs/RESULTS_SNAPSHOT.md](docs/RESULTS_SNAPSHOT.md) for the current working results;
- [docs/UPSTREAM_KRISTA.md](docs/UPSTREAM_KRISTA.md) for the relationship to Krista Pipho's original pipeline;
- [docs/PEGASUS_IMPORT_MANIFEST.md](docs/PEGASUS_IMPORT_MANIFEST.md) for scripts still to migrate from Pegasus.

## Current scientific picture

The final founder-model story is deliberately more nuanced than a raw three-versus-four winner.

The raw genome-wide likelihood favored the constrained three-founder model. Simulation calibration, however, showed that the strict three-founder generative model does not reproduce the observed statistic well. Broad four-founder simulations can contain the observation, while focused unrelated-four-founder simulations can overshoot it. Empirical full-sibling structure explains much of the intermediate signal but not all of it.

Genome-wide reciprocal mixing is approximately 18% in each direction. Chromosomes 1-20 are used for founder/mixing inference, while chromosome 21/Z is retained as an independent sex-ploidy diagnostic.

## Current Smudgeplot comparison

One universal lower cutoff is not appropriate for every dataset. The current preferred comparison is:

| dataset | role | preferred lower cutoff |
|---|---|---:|
| ERR1051012 | empirical diploid control | L=6 |
| ERR1307111 + ERR1307112 | pseudo-tetraploid benchmark | L=6 |
| TN Left | focal side | L=12 |
| TN Right | focal side | L=12 |
| TN T2 merge | left + right T2 | L=6 |
| TN T3 merge | left + right T3 | L=6 |
| all four TN legs | combined focal sample | L=12 |

Two calibration lessons are central:

1. ERR1051012 demonstrated that automatic 1n selection can choose a half-coverage harmonic and relabel the same balanced cloud from AB to AABB.
2. L=12 can preferentially remove the low-count arm of imbalanced k-mer pairs, particularly in merged/pseudo samples.

The later cutoff sweep therefore preferred L=6 for the diploid, pseudo-tetraploid, and cross-side merge comparisons where the lower arm must be retained.

## Upstream lineage

This work extends the earlier public repository:

Krista-Pipho/gynandromorph_2024

The provenance snapshot used here is commit:

65f4313299097fafed0df9e44654ddb98b6dcd61

GynandromorphCode is a separate private repository, not a literal GitHub fork. The exact relationship and inherited pipeline decisions are documented in [docs/UPSTREAM_KRISTA.md](docs/UPSTREAM_KRISTA.md).

## Repository layout

~~~text
.
├── README.md
├── docs/
│   ├── WORKFLOW.md
│   ├── AF_AND_FOUNDER_MODELS.md
│   ├── SMUDGEPLOT.md
│   ├── RESULTS_SNAPSHOT.md
│   ├── UPSTREAM_KRISTA.md
│   └── PEGASUS_IMPORT_MANIFEST.md
├── scripts/
│   ├── README.md
│   └── smudgeplot/
│       ├── README.md
│       └── plot_smudge_smu_viridis.py
└── .gitignore
~~~

## Reproducibility policy

This repository deliberately avoids making reusable code depend on Pegasus-specific absolute paths. Cluster paths, job IDs, FASTQs, BAMs, VCFs, FastK tables, .smu files, and large generated figures are not intended to be versioned.

For scripts imported from the working project:

- keep the original Pegasus path in a provenance header;
- replace numbered-script dependencies with descriptive filenames;
- keep exact scientific parameters documented;
- make reruns idempotent where practical;
- do not silently recompute expensive intermediates;
- prefer vector PDF/SVG output for figures;
- keep direct observations, fitted nuisance parameters, model predictions, and biological interpretation distinct.

## Historical Pegasus root

~~~text
/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
~~~

This path is recorded for provenance only.

## Status

The repository now contains the analysis narrative, current parameter decisions, upstream provenance, result snapshot, migration manifest, and a portable viridis Smudgeplot renderer.

The next code-migration priority is the founder-model fitting/calibration scripts, followed by the reciprocal-mixing AF figure scripts and the final lower-L Smudgeplot rerun scripts.
