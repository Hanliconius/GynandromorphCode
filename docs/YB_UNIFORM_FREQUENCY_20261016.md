# Yb SNP frequency genoplots from the shared uniform Rosina recall

This is a new, dated adaptation (2026-10-16). It does not alter the original
Optix caller, historical Yb SNP figures, phenotype files, or coverage analyses.
Implementation: `scripts/yb/yb_uniform_frequency.py` on branch
`optix-joint-recall-20261014`.

## Source and scope

Use `integration_Yb/Yb_window.vcf.gz` from the same shared Rosina run used for
the new Optix analysis. The actual comparative sample set is taken from the
frozen run configuration and must equal the VCF header; the old nominal 43
sample count is not imposed after sample exclusions. Every comparative
genotype comes from this new callset. No legacy SNP score or genotype is read.
There is no new comparative variant calling, download, alignment, or Hmel2.5
SNP analysis. Existing Hmel2.5 coverage results remain separate.

The run must have completed `optix_consistent_v2/corrected_calls.complete.json`
with status `AUDITED_CALLS_AND_MIXTURE_COMPLETE`. The archived model's SHA256
must equal the one recorded in that completion marker.

## Phenotypes and independent site selection

Metadata is joined by exact sample/specimen/accession identifiers. Dorsal and
ventral phenotypes are separate. User confirmations on 2026-10-16 override
metadata explicitly: cythera is dorsal absent / ventral absent; vulcanus is
dorsal absent / ventral present. These sources appear per sample in the new
phenotype manifest. Conflicting unresolved metadata causes an error; unknown
traits are not treated as absence and are excluded from that trait's scoring.
No phenotype implies a homozygous genotype.

At each biallelic nucleotide SNP, comparative diploid genotypes require DP>=2
and GQ>=20. Non-PK, non-focal controls alone determine dominant balanced
accuracy: compare ALT-carrier and REF-carrier orientations and retain the
higher score; omit orientation ties. Require >=80% callable in each phenotype
class and >=2 known individuals in each class. Output both BA>=0.80 and
BA>=0.90 panels. PK individuals are display-only.

Ybd analysis interval is NC_140403.1:1633041-1635646; Ybv is
1658079-1659072. The existing 1600000-1700000 window provides the wider view.
These are retained analysis bounds, not claimed causal deletion endpoints.
The Luca primer anchor span 1633657-1635225 is marked in the dorsal cartoon;
it is not labelled as the exact whole functional CRE or an observed deletion.

## Focal evidence and reciprocal genotype inference

Recount all biallelic Yb sites from the four existing focal BAMs using the
already committed `extend_optix_consistent_postcall.py` count worker unchanged.
A read-only input alias points its locus reader to the Yb VCF. It uses the same
bcftools executable and FASTA, MAPQ30, BQ20, default BAQ and per-BAM depth250
as the new Optix comparative/focal workflow. Overlapping indel records are
skipped and REF mismatches fail. These counts replace the older pysam focal
pileups for this figure only. T2 and T3 allele counts are summed within side.

Apply the unchanged archived nine-state reciprocal beta-binomial likelihood,
equal genotype-pair priors. Transfer **all five fitted parameters from the
completed NEW Optix fit**, not the historical rounded 20% estimate: reciprocal
contamination in each direction, error, ALT bias, concentration. Do not fit
these parameters on phenotype-selected Yb SNPs. This is explicitly a transfer
of the Optix fit to Yb; locus-specific bias/dispersion differences are not
estimated here and posterior confidence is conditional on this transfer.

Both pooled sides need >=2 REF+ALT observations and pair posterior>=0.95.
Retain all confident genotype pairs, including those opposing the biological
expectation. Non-confident calls remain missing; no inverted continuous allele
fraction is displayed as a genotype. The all-site table retains posterior,
genotypes, depths and other-base counts. A subsequent **display callability**
filter requires confident paired calls; it does not use the desired L/R state.

## Display and audit

One phenotype-positive and one phenotype-negative non-PK frequency track per
panel, including newly annotated races; two combined corrected focal tracks;
PK dark/pale and distinct dorsal/ventral phenotype combinations remain
separate. Frequencies are the mean called diploid dosage, not read-depth
weighted; >=80% of each group must be callable or the cell is grey. Focal cells
are discrete 0/0.5/1. The original light-blue / gold / purple palette is retained.
Consecutive SNP columns connect to the locus cartoon. Default maximum120 sites
per panel, evenly selected by genomic rank after filtering; full score tables
are retained. Four threshold/scope combinations for two traits make eight PDF
pages, including wider views and interval zooms.

Each run freezes the adapter, unchanged count worker, original caller, current
metadata, fitted parameter summary, fit completion marker, phenotype manifest
and file hashes. It saves independent SNP scores, corrected calls at every
site, per-panel selected keys, exact plotted matrices, group callability,
membership through the phenotype manifest, and `audit.json`.

## Execution

Invoke the script with `--submit --run /absolute/path/to/shared/run`. Preflight
checks imports, VCF cohort, focal BAM identity/index, model hash, fit bounds and
phenotype class sizes. It submits four independent focal count tasks and one
dependent analysis/PDF job on nano (30min, 2 CPUs, 12G each), without throttling.
The printed output gives job IDs, run directory and exact bundled PDF path.
If a count task reaches the nano limit, use its log to identify it; do not
resubmit the completed comparative recall or switch to stale focal counts.

Development validation uses synthetic tests for missing/quality handling,
REF/ALT orientation, equal-individual frequency and group callability, explicit
cythera/vulcanus phenotypes, PK exclusion from controls, and nine-state posterior
inference with low-depth no-calls. Real analysis and figure generation run only
on Pegasus. No new figure is claimed reviewed before inspecting its output.
