# Figure 4: NDJ butterfly coverage and paired allele fractions

This workflow uses **NDJ_left and NDJ_right**, not TN_Left/TN_Right. It renders
existing completed NDJ scan outputs on Pegasus. No downloads, alignments, variant
calling, new coverage computation, ploidy fitting, or mixing correction occur.

## Explicit source manifest

Frozen scan:
`/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/NDJ_LR_difference_scan_v1/run_20260920T012415Z`

Inputs are its `config.json`, `COMPLETE.json`, `windows.tsv`, and
`chromosomes/{index:02d}/counts.npz`, with indices resolved from the chromosome
configuration. Completion status and NDJ sample identities must match. Actual
input checksums, filters, ascertainment, window size, per-sample reuse decisions,
normalization and plotted values are exported. These are control-discovered SNPs,
not NDJ-private variant discovery. Depth windows derive from aligned-read depth,
not SNP allele depth; AF derives from the existing paired SNP AD counts.

Recovered source code/provenance:
`cbc9e087942d91776ecfa41ef3dfdd66c773051e` in Hanliconius/GynandromorphCode,
under `recovery/20261014T183358Z/source_1/NDJ_LR_difference_pipeline/` and
`recovery/20261014T183358Z/source_1/scripts/224_plot_NDJ_chr15_1n3n_followup.py`.
The earlier follow-up examined chromosome 15 positions 200,001-2,100,000.
Those exploratory boundaries are retained as dotted guides, not re-estimated
breakpoints. The new zoom extends to 0-3.1 Mb (0-based extent).

## Run

```bash
git -C /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo pull --quiet --ff-only origin main &&
micromamba run -n gynandrocalling python \
  /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo/scripts/figures/figure4/build_Figure4_NDJ_panels.py \
  --submit --overwrite
```

The nano job uses one CPU, 12 GB and a 30-minute limit. Output is a stable
directory, overwritten only with `--overwrite`:
`/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/Figure4_NDJ_panels`.
A failed render has no `complete.json`. Actual Pegasus PDFs still require visual
inspection before publication; local code validation is not figure validation.

## Outputs and presentation

- `01_B_chromosome_coverage.pdf`: chromosome-centred mirrored bars, black outlines,
  blue Left / orange Right, absolute tick labels, diploid guides, Fig2B typography.
  Bars are base-width-weighted means of usable cached windows, not raw whole-BAM
  totals. W was not part of this scan and is explicitly unassessed, not zero.
- `02_genome_delta.pdf`: black points, alternating pale chromosome backgrounds,
  red points for the earlier chromosome-15 interval. The statistic is normalized
  Right minus Left coverage, not an association score or a p-value.
- `03_chr15_delta.pdf`: corresponding chromosome-15 detail and ivory marker.
- `02_genome_relative.pdf`, `03_chr15_relative.pdf`: alternate views showing each
  side's normalized coverage directly, in the same blue/orange as B.
- `04_AF_Figure2_style.pdf`: two pages, whole genome then the earlier 0.2-2.1 Mb
  chromosome-15 region. Figure2H Luka-Bleute density, vector hexagons, log density
  midpoint compression, neutral LOESS marginals. No unsupported original-GT
  classifications, TN data, or TN ploidy/mixing predictive rings are substituted.
- `complete.json`, `reuse_decisions.tsv`, `coverage_windows_plotted.tsv.gz`,
  `chromosome_coverage.tsv`, and AF bin tables document every output.

All four coverage scatter panels are exactly **10 x 2.5 inches (4:1)** including
labels. `--ratio 5` makes them 10 x 2 inches. They are saved without tight cropping.
The fixed axes retain all usable points; extreme values are not silently clipped.

## Normalization and expected bands

For each side, existing relative depths are divided by their median in usable
autosomal windows, excluding chromosome 15 and sex chromosomes. This is a scalar
re-normalization of cached data, recorded in the audit; it differs from the old
follow-up's arithmetic-mean centering. Diploid-background relative coverage is 1;
delta is normalized Right minus normalized Left.

Blue-grey shading is the empirical central 95% of the observed autosomal
background (excluding chromosome 15). It includes technical and biological
variation: it is not a confidence interval or a calibrated null. The relative
view uses the union of the two side-specific envelopes.

Grey model bands assume reciprocal 1n/3n dosage relative to a 2n genomic
background, with independent 0-20% opposite-side DNA/read contribution measured
at that background. For copy numbers pL,pR:

`DL = (1-cL)*pL/2 + cL*pR/2`

`DR = cR*pL/2 + (1-cR)*pR/2`

Thus 1n gives 0.5-0.7, 3n gives 1.3-1.5, and the signed contrast is 0.6-1.0
(or -1.0 to -0.6 in the opposite orientation). These bands assume 2n background
and reciprocal segmental dosage; they are not estimates of mixing or proof of
ploidy. They do not cover every possible independent dosage combination.
Mixing fractions within an aneuploid segment's reads differ from background
mixing fractions; the AF figures therefore do not transplant the old fitted
local AF-centroid mixing parameters into these coverage bands.

Ivory uses the existing mapped-exon annotation in
`scripts/yb/Ybv_gene_features.json`, assembly GCF_056956725.1, NC_140403.1.
The marker spans the mapped exon extent 1,659,008-1,741,475; the annotation is
partial and is not represented as a complete gene model or a breakpoint.

## Targeted 03/04 refinement (2026-10-26)

Use `--submit --overwrite --refresh-03-04` with the same absolute script path.
This preserves the existing run spans and parameters, records checksums for every
01/02 PDF, and refuses to complete if those protected files change. It updates
only 03 coverage PDFs and the 04 AF PDF, in place. The original completion audit
is kept as `previous_complete.json`; the updated `complete.json` is written only
after success.

The cached source depth windows are 100 kb. Two independent nano array tasks
recount the same displayed chromosome-15 interval at 10 kb from the source NDJ
BAMs, using the original samtools depth base/map quality, excluded flags and mate
overlap handling. Dense tables have source/filter signatures and verified cache
reuse; no alignment, download or SNP calling occurs. Normalization remains the
same background normalization used in the first render. The plot retains every
usable 10-kb point and overlays local-linear tricube LOESS with a 200-kb full span
(100-kb bandwidth each side), requiring six local windows. Curves are broken
across gaps exceeding 100 kb. There are no confidence ribbons. Ivory and footer
labels now fit within the fixed canvas.

The AF correction changes rendering only: for the retained hex-bin assignment,
full hex height is **2*sy/3**, not **4*sy/3**. The old latter factor overlapped
neighbouring polygons. Marginal endpoint percentages move into a readable footer,
labels and colourbar are given separate space, and counts/AF orientation are
unchanged. Strong concentration in the AF-zero bins is not removed or recoded.

## Aligned C/C-prime stack and extended context

Use `--submit --overwrite --refresh-03-04 --stack --zoom-end 4100000`.
`05_C_Cprime_stacked.pdf` is one vector page: approved whole-genome delta view
above chr15 relative coverage with 10-kb points and 200-kb LOESS. Plot-area left
and right edges are identical. Each component retains its fixed 4:1 outer panel
shape; the combined page is 10 x 5 inches. Original 01/02 files stay checksum
protected, and this stack refresh leaves the AF PDF untouched.

A valid matching 10-kb cache ending at 3.1 Mb is extended by counting only
3.1-4.1 Mb from each existing BAM. Source/filter signatures and previous table
hashes must match before prefix reuse; the audit records reused prefixes and
newly counted intervals. Normalization and the earlier atypical-region guides
are unchanged. The separately saved 03 plots also show the new 0-4.1 Mb extent.

## Shorter smoothing span

The current renderer uses a **100-kb full LOESS span (50 kb each side)** rather
than the 200-kb span described for the first refinement above. It requires at
least four usable 10-kb windows per local fit and breaks curves across gaps
over 50 kb. Points, depth filters, normalization, mixing bands and the 0-4.1 Mb
zoom are unchanged. Existing validated depth caches remain reusable; shortening
LOESS alone requires no new counting. Use the same stack submission command.
