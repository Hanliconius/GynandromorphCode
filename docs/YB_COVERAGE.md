# Yb coverage workflows

## Matched local-read rebuild on both references (2026-10-12)

The current replication request is to rebuild four focal coverage tracks on
**both Hmel2.5 and HmRos**, using the same local read pairs and processing settings.
The panel-normalized analysis below is a separate earlier analysis, not the
processing recipe for this raw-coverage replication.

Submit from any directory on Pegasus:

```bash
bash /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo/scripts/yb/submit_Yb_local_both_references.sh
```

The submitter reads the four focal records from
`analysis/genomewide_mixing_input/metadata/resolved_download_manifest.tsv`.
It uses the HmRos FASTA from `metadata/reference_fasta.path` and the existing
`analysis/yb_exact_coordinate_transfer/reference/hmel2_5.chromosomes.fa`.
All eight FASTQs and both references were confirmed present on Pegasus.
No download stage is included.

Before submitting any jobs it checks all read/reference paths, native Yb contigs,
plotting packages, BWA, samtools, bamCoverage and Picard. Module setup messages
are retained in `analysis/yb_local_both_references/tool_setup.log`.
If needed, point `YB_PICARD_JAR` to an existing Picard jar.

Stages: two independent reference-index tasks, eight alignment/coverage tasks
(maximum four concurrently), then two plot tasks. All use the cpu partition;
time limits are 6 hours, 2 days and 30 minutes respectively. Job IDs print
immediately; a submission failure cancels jobs already submitted by that invocation.
Run-specific worker/config/plot code is frozen before submission.

Recipe follows upstream commit `65f4313299097fafed0df9e44654ddb98b6dcd61`:
BWA-MEM `-M`, coordinate sorting, Picard MarkDuplicates without removing duplicates,
and bamCoverage with 50-bp bins, no normalization, no MAPQ filter and no duplicate
exclusion. Full genomes are used for alignment. Reads remain independent by leg.
Fresh alignments are generated on both references for comparable processing.

Plots use stacked T2L/T2R/T3L/T3R tracks in blue/red/green/orange, natural-log
positive coverage, and native coordinates:
Hmel215003o:1,395,000-1,440,000 and NC_140403.1:1,620,000-1,670,000.
These are locus zooms; different window lengths do not imply identical basewise
homology. Hmel2.5 feature coordinates and transferred HmRos coordinates are kept
separate. The uncertain HmRos color-interval start remains an envelope.

Outputs are under `analysis/yb_local_both_references/<run>/`, with separate
`Hmel2.5` and `HmRos` directories containing deduplicated BAMs, 50-bp bigWigs,
metrics, vector PDFs, PNGs and audit summaries. The latest run path is recorded
in `analysis/yb_local_both_references/latest_run.path`. Unique run directories
avoid reusing or overwriting historical analyses. Temporary raw BAMs are removed
after successful processing; existing source files are unchanged.

The older `submit_Yb_raw_Krista_replication.sh` includes downloads and covers
only HmRos. Do not use it for this two-reference local rebuild.

## Earlier Yb / cortex normalized coverage

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


## Broad historical Hmel2.5 plot window

The current 45-kb Hmel2.5 Yb zoom is useful for the transferred locus markers,
but it is narrower than the broader historical plotting window discussed for
Krista's figure. A separate four-track rendering uses the existing full-scaffold
historical bigWig bedGraphs across `Hmel215003o:1,380,000-1,600,000`, rebinned to
50 bp, with linear and logarithmic-y PDF views and no additional per-library
normalization:

- code: `scripts/yb/plot_Yb_Hmel2_5_Krista_window.py`
- input: `analysis/yb_Hmel2_5_normalized_coverage/Hmel215003o_bedgraph/`
- submit on Pegasus: `bash scripts/yb/submit_Yb_Hmel2_5_Krista_window.sh`
- outputs: `analysis/yb_legacy_replication/Hmel2.5_Krista_window/`

These are historical bigWig values. They cannot recreate BAM-level filtering
that was applied before the bigWigs were made.
