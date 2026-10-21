# Integrated Optix/Dennis and ivory/Yb(v) figure and specimen review

Launch on Pegasus in `gynandrocalling`:

```bash
micromamba run -n gynandrocalling python \
/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo/scripts/figures/figure3/build_Optix_Ybv_integrated.py \
--submit --overwrite
```

Stable output: `/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/Optix_Ybv_integrated_figure`.

## Figure variants

* `Optix_Ybv_heterozygosity_annotation_genotypes.pdf`: BA above Ho, followed by gene/interval annotations, grey connectors and the condensed publication SNP panels.
* `Optix_Ybv_heterozygosity_annotation_genotypes_BA_dots.pdf`: the BA track is replaced by red dots at BA >=0.90 on the annotation track. Yb carrier and homozygote classifications occupy separately labelled dot lanes; they are never silently combined into one score.

Fixed genomic spans: Optix NC_140406.1:1060000-1240000; ivory/Yb NC_140403.1:1630000-1692000. Do not change these without author agreement. Ho remains 10-kb windows sliding by 2.5 kb, posterior .95/.99, common callable sites on both sides and >=20 accepted SNPs/window. The common Ho y-axis is tightened to the data with headroom. The native gene/exon GFF and verified partial ivory transfer are preserved; ivory continues beyond the viewport.

Both SNP panels use the Dennis orange/purple frequency palette. Allele meanings remain explicitly labelled: associated versus other allele in Dennis; ALT versus REF in Yb. Actual corrected heterozygotes show two letters in full-width top/bottom half-cells. A pooled frequency of exactly 0.5 is shown as `=`, not as an inferred individual heterozygote. Frozen publication SNPs, pooled values and focal calls remain unchanged. Elevatus replaces the linaresi morph example using actual same-phenotype elevatus sample genotypes; statistics retain the original cohorts.

## Why Yb BA was grey, and the extension

The earlier Yb BA read `ybv_REF_ALT_rebuild_v1` cached genotypes covering only about 3 kb, even though Ho covered 62 kb. Grey meant unscored span, not weak association. That earlier pipeline retained validated chromosome-15 BAMs.

`extend_Ybv_BA_region.py` reads the explicit reference/BAM list in the original rebuild run.json. It submits a BAM/reference/identity preflight, seven unthrottled genomic shards, and a scoring job. Every shard calls the same comparative BAM cohort with the existing Yb method: mpileup MAPQ30/BQ20/depth10000, call -m -A -f GQ. Shards have 200-bp context and nonoverlapping genomic cores. No downloads, alignment, focal counts or mixture fitting occur. The existing quality engine supplies DP/GQ or documented PL-derived quality; same-individual aliases are combined only for concordant calls.

Outputs stay in the stable `Ybv_BA_broad` subfolder. Shard checksums, BAM identities and reference provenance permit reuse on retries. Broad calls must agree with the old 3-kb cache to <=1% discordance among shared qualified calls, otherwise figure rendering is blocked for review. BA uses the unchanged v7 phenotype cohort with >=80% callability per class, separately displaying carrier and homozygote rules. Unscoreable SNPs stay missing, never become BA=0. No blanket shading is inferred from the first/last observed variant now that the full displayed span has been processed.

## Specimen review table

The independent nano review job runs before the wider BA calculation finishes. Outputs in `specimen_review`:

* `specimens_for_review.tsv`: import into Google Sheets with Tab delimiter.
* `specimens_for_review.csv`: equivalent CSV alternative.
* `column_dictionary.tsv`, `metadata_conflicts.tsv`, `README.txt`, `review_manifest.json`.

The exporter unions canonical metadata, all original 502 VCF IDs, the current Optix and Yb cohorts, the existing BAM inventory, the phenotype ledger, cached ENA identity mappings and the user-provided PK specimen table. One row represents each documented specimen, with explicit run aliases grouped. The TN tissue libraries are the one documented focal specimen. Identity conflicts, multiple BioSample IDs and potential duplicate cohort accessions remain flagged for author review. Metadata-only specimens are included and distinguishable from documented sequence resources.

Primary editable fields are names, taxon/population, B/D/R/Yb-dorsal/Yb-ventral/N/Sb/Br codes, sex where documented, BA_Fig3a and BA_Fig3b flags, and review_notes. Fig3a is Dennis; Fig3b is Yb(v). Flags record cohort membership, not callability at every SNP. Original source codes, effective codes, all raw canonical metadata, geography/locality, specimen/run/BioSample identifiers, source paths and 24-/27-marker callability are retained. Biological definitions of N/Sb/Br are left for author confirmation rather than guessed. Codes copied from the supplied inferred-genotype specimen table are identified separately and never substituted for additional measured phenotype data.

Please retain specimen_id when editing. Blank codes mean unknown, and CONFLICT means source records disagree. Returned edits are not applied automatically. The exporter protects either TSV or CSV if it has been manually edited.

## In-place updates and verification

`--overwrite` updates this workflow's generated PDFs in the same directory, blocks active recorded Slurm jobs and refuses a directory owned by another workflow. Source genotype panels and edited review files are preserved. Existing full-span Yb BA is reused when cohort provenance matches. Failed dependencies cancel downstream jobs instead of leaving them pending indefinitely.

Code checks cover identity/phenotype provenance, BA flags, edited-file protection, all seven shard boundaries and dependencies, concordant alias collapse, old-cache concordance, score export, pooled-tie versus actual-heterozygote representation, fixed spans and execution of both composition variants. Scientific computation and actual-data rendering run on Pegasus and require review after completion.
