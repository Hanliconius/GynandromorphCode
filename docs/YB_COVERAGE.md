# Yb / cortex normalized coverage replication

## Purpose

Replicate the historical coverage-based examination of the Yb/cortex region
using the exact four focal coverage tracks preserved in Krista Pipho's public
`gynandromorph_2024` Git history, then improve the comparison by normalizing
each library to its own local flanking coverage.

This analysis is independent of the reciprocal-mixture genotype caller.  It is
a direct read-depth / copy-number analysis.

## Historical source

Upstream repository:

`Krista-Pipho/gynandromorph_2024`

Historical snapshot used here:

`99376deebc25f4e0bb7aa9927c7969c693591ef4`

That snapshot contains the four original coverage bigWigs:

- `R_Markdown/TN2L_S1.dedup.bw`
- `R_Markdown/TN2R_S2.dedup.bw`
- `R_Markdown/TN3L_S3.dedup.bw`
- `R_Markdown/TN3R_S4.dedup.bw`

The historical `R_Markdown/Optix coverage plot.Rmd` imported these bigWigs
directly and plotted their coverage scores (on a log scale in the prototype).
It did not normalize each library for overall sequencing depth.

The same historical snapshot also contains
`202407_gynandromorph_snp_analysis/new_color_regions.tsv`, whose chromosome-15
records span `Hmel215003o:1,408,879-1,429,719`.

## Coordinate system and annotated intervals

The historical bigWigs use Hmel2.5 scaffold coordinates.  The first pass plots
`Hmel215003o:1,380,000-1,600,000`, which contains the Yb/cortex region and
provides flanking sequence for normalization.

Annotated intervals retained in the audit/plot are:

- Yb(d): `Hmel215003o:1,403,500-1,405,500`
- historical Krista colour-region variants:
  `Hmel215003o:1,408,879-1,429,719`
- Yb(v): `Hmel215003o:1,429,500-1,430,500`
- ivory Δ78k deletion:
  `Hmel215003o:1,494,902-1,573,177`

The ivory interval is a known copy-number diagnostic: homozygous deletion
coverage should approach zero and heterozygous deletion coverage should
approach one half of the local diploid baseline.

## Normalization

Coverage is binned at 250 bp.  For each library independently, the raw coverage
profile is divided by the median raw coverage across two local flanking windows:

- `1,380,000-1,395,000`
- `1,585,000-1,600,000`

Thus ordinary local diploid coverage is centered near 1.0.  A one-copy deletion
is expected near 0.5 and a two-copy deletion near 0, subject to mapping noise.

The normalization changes only scale; it does not use phenotype or genotype
calls to select positions.

## Outputs

Pegasus output directory:

`analysis/yb_normalized_coverage/results/`

Primary files:

- `Yb_normalized_coverage.pdf`
- `Yb_normalized_coverage.summary.tsv`
- `Yb_normalized_coverage_250bp.long.tsv.gz`
- `Yb_normalized_coverage.audit.json`

## Relationship to focal genotype re-calling

Optix focal genotype figures use the reciprocal-mixture caller and display only
Left/Right genotype pairs with posterior >= 0.95, with no fallback to original
GT/GQ calls.  This Yb coverage analysis does not use those calls: it returns to
the original read-depth evidence preserved in the upstream bigWig tracks.
