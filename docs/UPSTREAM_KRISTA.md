# Upstream Krista Pipho workflow

The present repository is a continuation and extension of the earlier public repository:

Krista-Pipho/gynandromorph_2024

Pinned upstream commit used for provenance:

65f4313299097fafed0df9e44654ddb98b6dcd61

The upstream repository describes a Snakemake whole-genome workflow for the H. melpomene gynandromorph, including alignment, QC, coverage, variant calling, phasing, and exploratory ploidy analyses.

## Important upstream files

- README.md: pipeline description and original tutorial
- Snakefile: main Snakemake workflow
- config.yaml: focal TN sample set and older candidate control accessions
- alignment.sh: alignment stage
- qc_alignment.sh: BAM QC/dedup stage
- coverage.sh and scaffold_coverage.sh: coverage analyses
- per_sample_variant_calling.sh: per-sample GVCF generation
- phase_variants.sh: haplotagging / phasing stage
- z_chr_variants_analysis.sh: Z-chromosome exploratory analysis
- smudge_launch.sh: original FastK + Smudgeplot launch script

## Original Smudgeplot settings

At the pinned commit, smudge_launch.sh used:

FastK -v -t4 -k22 -M100 -T60 SAMPLE_DIR/*.fastq -NSAMPLE_DIR/FastK_Table

smudgeplot hetmers -L 12 -t 60 -o SAMPLE_DIR/kmerpairs --verbose SAMPLE_DIR/FastK_Table

smudgeplot all -o SAMPLE_DIR/output SAMPLE_DIR/kmerpairs.smu

The upstream config listed the four TN libraries and, in commented alternative configuration, the control accessions:

ERR1050907
ERR1050958
ERR1050994
ERR1051012

These details are why the later control calibration began with ERR1051012 and why L=12 was treated as the historical baseline rather than an arbitrary new parameter.

## Relationship to GynandromorphCode

GynandromorphCode is not a literal GitHub fork because it was created as a separate private repository and now contains substantial downstream analyses that did not exist upstream.

The useful provenance model is:

Krista upstream workflow
    -> reproduce / audit historical stages
    -> Joe's Pegasus analyses
    -> cleaned descriptive scripts in GynandromorphCode

The upstream commit is therefore pinned here so every inherited decision can be traced.

## Why the full upstream repository is not copied here yet

The upstream repository is public, but it does not expose an explicit license file at the pinned commit. Rather than silently republishing the full source tree inside a different repository, this project currently records the exact upstream repository and commit and imports only newly written or substantially adapted downstream code.

If a literal GitHub fork is desired, it should remain a separate fork of Krista-Pipho/gynandromorph_2024 and be linked from this repository.
