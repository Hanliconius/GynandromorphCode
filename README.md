# GynandromorphCode

Reproducibility repository for genomic analyses of a *Heliconius* gynandromorph, with emphasis on left/right tissue differences, reciprocal read mixing, founder-state models, chromosome-scale coverage, and k-mer-pair (Smudgeplot) diagnostics.

This repository is intended to preserve both the **final analysis path** and the **reasoning that led to it**. Exploratory branches and failed parameterizations are documented when they changed how the final analysis was interpreted.

## Biological question

The focal individual was sequenced from left and right tissue, with T2 and T3 libraries available for each side. The central questions are:

1. Do the left and right sides differ genetically?
2. Can those differences be explained by reciprocal mixing of reads between sides?
3. Are the autosomal allele-frequency patterns more consistent with two, three, or four founder chromosomes?
4. Do chromosome-scale depth and k-mer-pair patterns support the same interpretation?
5. Which apparent Smudgeplot classes are biological, and which are artifacts of cutoff or coverage-harmonic choice?

## Focal samples

| role | sample |
|---|---|
| Left T2 | SRR39053777 / TN_Left_T2 |
| Right T2 | SRR39053776 / TN_Right_T2 |
| Left T3 | SRR39053775 / TN_Left_T3 |
| Right T3 | SRR39053774 / TN_Right_T3 |
| current diploid Smudgeplot control | ERR1051012 |
| pseudo-tetraploid benchmark | ERR1307111 + ERR1307112 |

## Analysis map

The current workflow is:

1. **Prepare and QC the four focal libraries.**
2. **Measure left/right allele-frequency differences** using read counts rather than treating hard genotype calls as truth.
3. **Estimate reciprocal mixing** between the two sides.
4. **Compare founder models** (especially three-founder versus four-founder models) using simulations and fitted allele-frequency distributions.
5. **Use reciprocal-mixture genotype calls** for focal display only when posterior probability is at least 0.95; there is no fallback to lower-confidence calls.
6. **Inspect chromosome-scale relative coverage**, including autosomes and sex chromosomes.
7. **Run Smudgeplot arrangements** for Left, Right, T2-merged, T3-merged, all four legs, a diploid control, and a pseudo-tetraploid benchmark.
8. **Calibrate Smudgeplot cutoff and 1n coverage explicitly** rather than assuming one cutoff or trusting an automatically chosen harmonic.
9. **Generate vector figures** (PDF/SVG preferred) and retain diagnostic reports alongside publication-style summaries.

See [docs/WORKFLOW.md](docs/WORKFLOW.md) for the stepwise analysis and [docs/SMUDGEPLOT.md](docs/SMUDGEPLOT.md) for the k-mer workflow.

## Current high-confidence Smudgeplot comparison

The final comparison should not use one universal lower cutoff. The current preferred settings are:

| dataset | role | lower cutoff |
|---|---|---:|
| ERR1051012 | diploid control | L=6 |
| ERR1307111 + ERR1307112 | pseudo-tetraploid benchmark | L=6 |
| TN Left | focal side | L=12 |
| TN Right | focal side | L=12 |
| TN T2 merge | left + right T2 | L=6 |
| TN T3 merge | left + right T3 | L=6 |
| all four TN legs | combined focal sample | L=12 |

The earlier ERR1051012 L=4 / manual 1n=8.01 run was useful diagnostically because it demonstrated that Smudgeplot can choose a half-coverage harmonic and relabel a balanced cloud from AB to AABB. A subsequent cutoff sweep showed that L=6 with automatic coverage is the cleaner final diploid-control setting.

## Repository layout

~~~text
.
├── README.md
├── docs/
│   ├── WORKFLOW.md
│   ├── SMUDGEPLOT.md
│   ├── AF_AND_FOUNDER_MODELS.md
│   ├── RESULTS_SNAPSHOT.md
│   └── PEGASUS_IMPORT_MANIFEST.md
├── scripts/
│   ├── README.md
│   └── smudgeplot/
│       ├── README.md
│       └── plot_smudge_smu_viridis.py
└── .gitignore
~~~

## Reproducibility policy

This repository deliberately avoids embedding Pegasus-specific absolute paths into reusable code. Cluster paths, job IDs, large intermediate files, FASTQs, BAMs, VCFs, FastK tables, and figure products are not intended to be versioned.

For scripts imported from the working Pegasus project:

- replace hard-coded paths with command-line arguments or environment variables;
- replace numbered script dependencies with descriptive filenames;
- keep exact analysis parameters in the documentation;
- make reruns idempotent where practical;
- do not silently re-extract VCFs or overwrite completed intermediates;
- prefer vector PDF/SVG output for figures.

## Working project location on Pegasus

The historical working tree is:

~~~text
/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
~~~

That path is recorded for provenance only; repository scripts should not depend on it.

## Status

This repository currently contains the analysis narrative, parameter decisions, result snapshot, and a portable viridis Smudgeplot renderer. The remaining Pegasus scripts are listed in [docs/PEGASUS_IMPORT_MANIFEST.md](docs/PEGASUS_IMPORT_MANIFEST.md) so they can be copied in and cleaned without losing provenance.
