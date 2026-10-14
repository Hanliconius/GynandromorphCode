# Optix/Yb 2026 all-sample joint recall — implementation contract

## Why the old VCF cannot be treated as a joint callset

The October 13 post-hoc PK extension used frozen Optix sites / legacy calls for some individuals and new calls for others. The observed ~30.5–31k 0/0 calls with missing GQ are not automatically high-confidence reference genotypes. The saved PK PDF is diagnostic, **not** the basis for the replacement analysis.

## Exact scope

- 20 retained PK dark/pale individuals from scripts/figures/figure3/pk_dark_pale_samples.tsv; exclude P4009 / SRR11243345, retain P2010 / SRR11243362.
- All earlier race controls in the mean-coverage run plus the four supplemental rosina/melpomene controls (ERR2298205, ERR2298206, ERR2298204, ERR2298214), not just samples in the legacy Optix VCF.
- Four focal libraries, T2L/T3L and T2R/T3R, to be pooled on each side for coverage/read-count display, **not** treated as four independent diploid individuals.
- Add 3–5 unique *H. m. cythera* specimens and 3–5 unique *H. m. vulcanus* specimens, identified by ENA BioSample, one selected WGS run per specimen unless a multi-run merge is independently audited.

## Starting archival candidates (NOT yet the final 3–5 each)

Martin et al., Genetics (2016), 203:525–541, table of resequenced butterflies, lists:
- cythera: cyth.CJ2856 (ERS977691), cyth.CJ2857 (ERS977692).
- vulcanus: vul.CJ14632 (ERS977709), vul.CS519 (ERS977710), vul.CS10 (ERS1030540), vul.CS11 (ERS1030541).

These are specimen/BioSample accessions, **not** run accessions and not proof that current ENA FASTQs are accessible. At least one further independent cythera BioSample is required. Additional candidates should preferentially be traced from the pre-existing analysis-142 cohort (10 cythera and 11 vulcanus as previously recorded), then validated against ENA. Do not manufacture run IDs or count technical lanes as individuals.

Luca Livraghi et al. eLife (2021), 10:e68549, Figure 7 includes cythera/vulcanus in topology analyses, but Figure 8 source data 1 is a *different* coverage panel and does not provide the complete cythera/vulcanus WGS sample inventory.

## Reference invariant

The common FASTA is the Rosina assembly used in the corrected PK run: Optix contig NC_140406.1 (old synonym Hmr18), Yb contig NC_140403.1. Check exact REF bases against the old frozen Optix allele table before submitting. A contig-name alias does **not** constitute a liftover. Keep Hmel2.5 alignments/callsets separate.

## End-to-end implementation

1. Freeze a single accession-to-specimen provenance manifest and a single FASTA checksum. Fail on duplicates, non-WGS, unexpected SRA layouts, missing mate/checksum information, ambiguous species/race, or Optix REF mismatch.
2. Inventory existing read files and alignments without changing running jobs. Reuse FASTQs only if ENA run metadata and checksums match. Reuse a BAM **only** with an auditable matching FASTA digest, read group, aligner/duplicate recipe, and completed index; otherwise realign.
3. Download only genuinely missing inputs with Slurm arrays, strict ENA MD5/size checks, partial files and atomic rename. Cap concurrent downloads to avoid saturating shared storage/network.
4. Rebuild all alignments to the *same entire Rosina FASTA*, using frozen read group IDs, BWA-MEM, coordinate sorting, duplicate marking and indexing, with per-library provenance. One Rosina BAM feeds both loci.
5. Use a unified diploid calling workflow for all independent comparative specimens. **Prefer per-sample GATK HaplotypeCaller GVCFs then joint GenotypeGVCFs**, including reference genotype quality, rather than assuming bcftools -f GQ fills every homozygous-reference call (it need not). Freeze versions, MAPQ/BQ rules, GQ/DP fields and filters. Run Optix and Yb as parallel locus shards, then validate cohort/sample concordance.
6. Keep all-site likelihood/callability data as well as the variable-SNP VCF, including exact REF/ALT provenance. Missing and uncertain are different from homozygous reference. No legacy comparative GT values copied into the new callset. Site ascertainment excludes focal TN libraries, and phenotype labels are not used in the caller.
7. At jointly established Optix sites, extract per-library focal allele depths with the same alignment/read quality rules; pool T2+T3 **within** each side for display and fit reciprocal-mixing posteriors independently. Do not force diploid GATK GT calls onto admixture-contaminated tissues.
8. Produce sample-level QC: source accession, BioSample, FASTQ digest, read/mapping metrics, coverage, callable fraction, reference-bias audit, 0/0 vs 0/1 vs 1/1 vs missing breakdown, DP/GQ distributions and optix frozen-site overlap. Fail the final report if required samples or genotyping attributes are absent.
9. Generate new vector-PDF Optix B/D/R / PK and Yb comparison plots. Never overwrite the October 13 diagnostic PDF or reuse its mixed callset.

## Cluster submission

Use Peggy/Pegasus arrays with bounded parallelism: downloads (1 core, e.g. 8–12 concurrent), alignments (16 cores, 64 GB, e.g. at most 4 concurrent), per-sample HaplotypeCaller per locus (e.g. 4 cores, bounded array), joint import/genotyping per locus (independent Slurm jobs), and final audit/plot with afterok dependencies. Size limits from measured inputs and partition quotas, not assumptions. Snapshot manifest, FASTA digest, source, tool versions and Slurm script inside timestamped run directory. Do not release downstream jobs on incomplete upstream tasks. Do not schedule from an unvalidated inventory.

## Existing inputs to audit first

- PK: /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/pk_dark_pale_optix_yb/20261013T220116Z_0b6e0845
- Focal/control mean run: /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_mixing_and_allele_evidence/20261013T194506Z_10e7bd9e (verify actual run parent path)
- Extra controls: run 20261013T201033Z_b28c41fe (verify actual run parent path)
- Current reference, aliases and completed BAMs: corrected PK run + latest yb_local_both_references run
- Original analysis-142 accession / sample table: locate on Pegasus and extract cythera/vulcanus specimens, not site calls

## Current state

This is the reproducibility contract / candidate inventory, **not evidence that the full callset was created**. The next run must first verify full accession inventory and Rosina reference. Until then, do not claim final sample counts, GQ completeness, or submitted Slurm jobs.
