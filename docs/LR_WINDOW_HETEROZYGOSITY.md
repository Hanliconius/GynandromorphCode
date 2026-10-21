# Left/Right regional heterozygosity exploration

This workflow measures the heterozygous fraction among callable discovery SNPs,
not nucleotide diversity per callable base. SNPs are not selected by balanced
accuracy, focal genotype pattern, or phenotype. Existing comparative discovery
still defines ascertainment; differences between discovery cohorts matter.

## Sources and preservation

Project root: `/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling`.

Optix reuses the broad corrected pair table in
`analysis/rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc/optix_consistent_v2/reciprocal_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz`.

Ivory/Yb uses the broad `integration_Yb/Yb_chromosome.vcf.gz` from that same
completed recall, with reference-verified biallelic SNPs at NC_140403.1:1630000-1692000.
It counts the four existing focal BAMs recorded in the Optix postcall run.json,
using the matching Rosina FASTA, MAPQ30/BQ20/default BAQ/depth250 per library.
Counts are summed by side. The archived nine-state reciprocal beta-binomial
likelihood is applied with parameters FIXED to the summary used by the earlier
Yb restoration:
`analysis/original_502/reciprocal_L_R_mixture/optix_reciprocal_L_R_mixture.summary.tsv`.
This extends focal inference to broad Yb SNPs, without fitting parameters again.
Uniform genotype-pair priors are used, as in the archived model.
These are NEW EXPLORATORY Yb calls, not a claim that every old short-marker
call or confidence value is reproduced. Prior calls and figures are never edited.

Optix and Yb retain their respective parameter fits; do not directly compare
their absolute heterozygosity as if the calling procedures were identical.
Within each locus the two sides have exactly the same accepted denominator.
This exploration is conditional on confident calls, so the accepted-site counts
and 0.99 confidence sensitivity must be inspected alongside the primary curves.
Neither linked SNPs nor overlapping windows are independent replicates.

## Windows and quality

* Optix: 1, 5 and 10 kb windows across 1.06-1.24 Mb.
* Ivory/Yb: 1, 2 and 5 kb windows across 1.63-1.692 Mb.
* Both: 50 and 100 discovery-SNP windows.
* Steps are one quarter of window size (integer floor for SNP windows).
* Depth >=10 reference+target ALT reads on EACH pooled side; both diploid GTs
  present; pair posterior >=0.95, with >=0.99 sensitivity on the same pages.
* Yb sites with >5% other-allele reads on either side are excluded.
* At least 20 accepted SNPs per window; otherwise the curve is missing, not zero.
* Physical windows include only complete windows inside the requested region.
* Fixed-SNP boundaries use discovery sites before confidence/depth filtering.

## Launch on Pegasus

```bash
git -C /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo pull --quiet --ff-only origin main &&
micromamba run -n gynandrocalling python \
  /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo/scripts/figures/figure3/plot_LR_window_heterozygosity.py \
  --submit
```

Four unthrottled count tasks run on cpu, followed by a nano report dependency.
Preflight verifies reference alleles, BAM dictionaries, read-group identities,
indices, and documented input paths. Existing output directories are not reused.
No downloads, alignments, comparative re-calling, or parameter optimization occur.

Default output is `analysis/LR_window_heterozygosity_v1` under the project root:

* `Optix_LR_heterozygosity_versions.pdf` (five pages)
* `Ivory_Yb_LR_heterozygosity_versions.pdf` (five pages)
* matching per-window TSVs and PNG previews
* new exploratory Yb genotype pairs, frozen model, parameter/source hashes,
  explicit per-library count commands, job manifest and final completion JSON

Actual-data rendering and density suitability require review after the job runs.
Local numerical checks cover known Ho values, paired denominators, masking,
discovery-window boundaries, REF/ALT invariance and recovery of all nine archived
genotype-pair states. No biological results have been inferred locally.
