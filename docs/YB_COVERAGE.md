# Yb / cortex normalized coverage

## Current analysis reference

The current Yb coverage analysis uses the chromosome-level *Heliconius melpomene rosina* RefSeq assembly:

- accession: `GCF_056956725.1`
- assembly: `ASM5695672v1`

This is the same HmRos reference used by the genome-wide focal/control alignments in
`analysis/genomewide_mixing_input/`.

The historical Hmel2.5 analysis is retained only for provenance and must not be mixed with the HmRos coordinates.

## What was recoverable from Krista Pipho's Git history

The historical repository `Krista-Pipho/gynandromorph_2024` contains Hmel2.5-based work, not HmRos coordinates.

Useful historical evidence includes:

- four original focal coverage bigWigs;
- chromosome-15 colour-region variants spanning
  `Hmel215003o:1,408,879-1,429,719`;
- the broader Hmel2.5 Yb/cortex context.

Published Hmel2.5 intervals are also retained as provenance:

- Yb(d): `Hmel215003o:1,403,500-1,405,500`
- Yb(v): `Hmel215003o:1,429,500-1,430,500`
- ivory Δ78k: `Hmel215003o:1,494,902-1,573,177`

These coordinates are **not** reused directly on HmRos.

## Resolving the HmRos locus

The HmRos workflow downloads/reuses the RefSeq GFF associated with
`GCF_056956725.1`, identifies chromosome 15 from the NCBI assembly report,
and then requires exactly one annotated `cortex` gene on chromosome 15.

The initial HmRos coverage view spans the cortex gene plus 125 kb on each side.
The resolved chromosome/accession and coordinates are written to:

`analysis/yb_hmros_normalized_coverage/metadata/Yb_hmros_coordinates.tsv`

This gives a non-guessed coordinate anchor on the new assembly.

The historical Yb(d), Yb(v), and Krista 2024 colour-region coordinates remain
recorded in that table as provenance only. If exact base-for-base homologous
boundaries are required for the final figure, they should be sequence-lifted
onto HmRos rather than transferred by arithmetic offset.

## Coverage inputs

No historical Hmel2.5 bigWigs are used in the current analysis.

The workflow uses the already-existing deduplicated BAMs aligned to
`GCF_056956725.1`:

- TN_Left_T2
- TN_Left_T3
- TN_Right_T2
- TN_Right_T3
- 18 independently selected WGS controls

The BAMs are under:

`analysis/genomewide_mixing_input/alignment/libraries/`

No reads are downloaded or realigned for this analysis.

## Normalization

The coverage analysis is designed as a panel-of-normals copy-depth analysis,
not as a raw-coverage plot.

### 1. Fixed-bin target depth

Coverage is measured in 250-bp bins across the HmRos cortex neighbourhood.

Only reads/bases passing:

- mapping quality >= 20
- base quality >= 20

contribute to depth.

Zero-depth positions are retained.

### 2. Library-size normalization

Each sample is scaled using its deduplicated mapped-base total from the existing
`samtools stats` output.

For sample *s*:

`library_scaled_depth = raw_bin_depth × median(mapped_bases) / mapped_bases_s`

This removes differences in total usable sequence yield between libraries.

### 3. Local panel-of-normals normalization

For each genomic bin, the median library-scaled depth across the 18 controls is
used as the local expected diploid depth.

For a focal library:

`panel_normalized_depth = focal_library_scaled_depth / median_control_scaled_depth`

Thus:

- ~1.0 = ordinary diploid depth;
- ~0.5 = approximately one-copy depth;
- ~0 = approximately zero-copy depth.

This second normalization corrects much of the local GC/mappability structure
that a single global depth factor cannot remove.

### 4. Unstable-bin masking

A bin is omitted from focal interpretation if the control panel is unreliable.

Current rules require:

- at least 10 controls;
- control median >= 20% of the typical local control depth;
- control relative MAD <= 0.35.

The complete control-panel QC table is retained.

### 5. Replicate handling

T2 and T3 are normalized independently.

Only after normalization are the two libraries averaged to make a Left-side or
Right-side mean. Raw T2/T3 tracks are retained in the plot so replicate
agreement remains visible.

This avoids allowing the deeper library to dominate a side merely because it
contains more reads.

## Outputs

Primary output directory:

`analysis/yb_hmros_normalized_coverage/results/`

Outputs:

- `Yb_hmros_panel_normalized_coverage.pdf`
- `Yb_hmros_panel_normalized_coverage.summary.tsv`
- `Yb_hmros_panel_normalized_250bp.long.tsv.gz`
- `Yb_hmros_side_means_250bp.tsv.gz`
- `Yb_hmros_control_panel_qc.tsv.gz`
- `Yb_hmros_panel_normalized_coverage.audit.json`

## Historical workflow

The earlier `scripts/yb/plot_Yb_normalized_coverage.py` workflow used the
historical Hmel2.5 bigWigs. It is preserved only as a historical replication
path and is not the current analysis.


## Legacy-style HmRos companion view

A four-track HmRos companion to the historical Hmel2.5 Yb plot uses the
existing 250-bp focal depth table and does not realign reads:

- code: `scripts/yb/plot_Yb_HmRos_legacy_replication.py`
- input: `analysis/yb_hmros_normalized_coverage/results/Yb_hmros_panel_normalized_250bp.long.tsv.gz`
- full-chromosome mean-depth calculation: existing focal BAMs, MAPQ>=20 and baseQ>=20
- outputs: `analysis/yb_legacy_replication/HmRos/`
- submit on Pegasus: `bash scripts/yb/submit_Yb_HmRos_legacy_replication.sh`

The tracks preserve the historical order and colors on the native HmRos
`NC_140403.1:1,620,000-1,670,000` axis. Each bin is divided by its sample's
mean depth across the full HmRos chromosome. HmRos uses 250-bp bins from the
existing MAPQ/baseQ-filtered focal alignments; Hmel2.5 uses historical 50-bp
bigWig tracks divided by the Hmel scaffold mean. This makes the baselines more
comparable, but differing assemblies, alignments, and bin sizes still affect
the traces.

Yb(d), Yb(v), and the cortex TSS use transferred HmRos coordinates. The start
of the historical Krista color interval is only bracketed at
`1,638,549-1,640,834`; the plot marks that uncertainty rather than assigning
an exact boundary.
