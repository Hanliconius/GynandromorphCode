# Integrated Optix/Dennis and ivory/Yb(v) figure

Worker: `scripts/figures/figure3/build_Optix_Ybv_integrated.py --submit --overwrite` on Pegasus, environment `gynandrocalling`.

Stable output folder: `/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/Optix_Ybv_integrated_figure`.

Single vector PDF: `Optix_Ybv_heterozygosity_annotation_genotypes.pdf`, with Optix above ivory/Yb(v). Each block has a full-span phenotype BA plot ABOVE the heterozygosity plot, a gene/interval annotation track, grey lines to the condensed publication SNP columns, and the consensus/morph/corrected-genotype rows. Left/Right rows retain 2.7-fold height and unphased top/bottom rectangular half-cells. No nucleotide-specific colour scheme.

**Genomic spans must not be changed without user agreement:**

* Optix NC_140406.1: 1,060,000-1,240,000.
* Ivory/Yb NC_140403.1: 1,630,000-1,692,000.

Heterozygosity uses 10,000-bp windows sliding 2,500 bp, with the existing worker's depth, posterior .95/.99 and >=20-SNP rules. It reuses corrected Optix pairs and the already completed exploratory fixed-parameter Yb pairs from `LR_window_heterozygosity_v1`; no inference or optimization occurs in this figure workflow.

The SNP focus panels reconstruct the established publication versions from `Dennis_figure_set_v1` and `Ybv_figure_set_v1` cached dosage matrices, frozen Dennis frequencies and `ybv_short_ALT_geography_v7` membership. They validate pooled values against those figures. Dennis uses associated-allele purple versus other-allele orange. Yb uses ALT purple versus REF orange. No cohort, allele polarity, marker selection or focal genotype changes are allowed. Morph groups begin with the four largest per phenotype. At the user's request, a displayed linaresi/linearesi group is replaced by an existing same-phenotype elevatus group, preferring the same geography and otherwise the largest available elevatus group. Actual sample IDs and replacement rationale are recorded in complete.json. This changes display examples only, never consensus membership, BA, Ho or focal calls.

The native gene track reads the documented GCF_056956725.1 GFF at `analysis/genomewide_mixing_input/reference/GCF_056956725.1/GCF_056956725.1_ASM5695672v1_genomic.gff.gz`, drawing overlapping native genes and exons with native labels and strand arrows. Yb supplements these with the previously verified cortex and partial ivory exon transfer in `scripts/yb/Ybv_gene_features.json`. Ivory continues beyond the preserved viewport and is labelled accordingly, rather than extending the horizontal axis. The script does not infer an optix label for an unnamed native model. Actual native labels need inspection after the Pegasus rendering.

`--overwrite` refreshes this workflow's output in place, including failed runs. It refuses a directory owned by another workflow or a currently active recorded Slurm job. Source inputs remain untouched. A completion JSON records spans, rows and displayed gene models; run.json records sources and hashes.

Numerical tests cover exact publication reconstruction and pooled frequencies, corrected genotype polarity, gene-to-transcript-to-exon parenting and the fixed spans. Actual-data PDF rendering remains a Pegasus verification step.

## BA and y-axis update

BA reuses the exact scoring helper `stack_LR_heterozygosity_BA.py`, frozen as ba_worker.py. Dennis reuses its existing all-SNP scores. Yb displays the previously specified carrier and homozygote rules separately on the same track, using cached expanded-cohort genotypes and the retained v7 membership. Yb flanks without expanded-cohort BA are hatched. There is no new calling or mixture inference.

Both heterozygosity tracks use one common tighter y-axis from zero, with the upper bound rounded upward in 0.1 units with 10% headroom (minimum 0.3). Values are not clipped to force visual agreement. Genomic x-axis spans remain unchanged. The page is taller to accommodate BA; the relative bold focal-strip heights remain 2.7-fold. Heterozygotes now show two full-width top/bottom rectangles, with centered letters and no diagonal polygons. These halves do not imply phasing.
