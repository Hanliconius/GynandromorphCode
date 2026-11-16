# NDJ chromosome 15: sequence context and rearrangement signatures

## Question and recovered case context

Does the chromosome 15 coverage transition in the **NDJ_left / NDJ_right segmental aneuploid individual** lie in sequence with properties associated with chromosome breakage or rearrangement? This is separate from the TN four-library founder/Smudgeplot analysis.

The reference is Rosina **GCF_056956725.1**, chromosome 15 **NC_140403.1**. Inputs come from the frozen `NDJ_LR_difference_scan_v1/run_20260920T012415Z/config.json`: `reference`, `bam_left`, `bam_right`, chromosome lengths and sample identities. The correct BAMs are the realigned `NDJ_control_comparison_v1/alignment/NDJ_{left,right}/NDJ_{left,right}.dedup.bam`; older GynL/GynR alignments used a different reference and must not be substituted.

The default target comes from `Figure4_NDJ_panels/chr15_descriptive_segments.json`, field `transition_bracket_bp`. Figure4 estimates a shared descriptive transition using usable 10-kb depth windows. **That bracket is neither a confidence interval nor a sequence-resolved breakpoint.** The older 0.2–2.1-Mb highlighted interval is not used to invent an exact breakpoint. Exact coordinates remain a runtime input from the existing Pegasus result; none were available in the development workspace.

No biological result is claimed by this code delivery. The first run asks whether the uncertain interval has unusual reference sequence, and extracts evidence that may narrow the junction.

## Literature-led feature panel

| Property | Why score it? | Implemented measurement and interpretive limit |
|---|---|---|
| TE/repeat content | In *Leptidea*, fusion boundaries were enriched for LINEs and LTRs; fissions did not show the same pattern [1]. | Optional same-reference RepeatMasker class union coverage; preserve family, orientation and divergence in a separate details table. Repeat abundance alone is not NAHR. |
| Long homology, orientation, duplicated sequence | Recombination between nonallelic homologous segments can generate structural changes; repeat orientation helps constrain possible products [2,3]. | Optional BLAST of overlapping local 5-kb reference queries against the full reference; report nonself partner, identity, aligned length and relative orientation. This is candidate discovery, not an enrichment test or proof that a partner participated. |
| Short direct/inverted/mirror repeats | Potential slipped structures, hairpins/cruciforms and triplex structures; experiments and breakpoint surveys support their relevance [4,5]. | Exact 10-bp arms, spacer 0–100 bp, nearest eligible upstream arm. Mirror subset with >=90% purine or pyrimidine arm and <=8-bp spacer is a triplex proxy. Seed screens are not exhaustive structural predictions. |
| G-rich quadruplex motifs | Experimental yeast work links G4-forming sequences and instability under specific conditions [6]. | Both strands, four runs of 3–7 Gs (or Cs), 1–7-base loops; fraction covered by motif spans. No assertion that the structure formed in the butterfly. |
| Alternating purine/pyrimidine sequence | Z-DNA-forming sequences can induce large deletions experimentally [7]. | Alternating R/Y tracts >=12 bp; broad proxy, deliberately not called a validated Z-DNA prediction. Sequence energetics and cellular conditions matter. |
| Tandem repeats, homopolymers, low complexity | Slippage and secondary structure can impede replication/repair; repeats can also compromise alignment [4,5]. | Exact tandem units 2–6 bp, >=3 copies and >=12 bp total; homopolymers >=8 bp; base-composition entropy; softmask fraction. Entropy is a composition summary, not DUST or a full complexity model. |
| AT-rich runs and composition | Some fragile-site experiments implicate fork stalling at AT-rich sequence, while others implicate cell-specific origin scarcity [8,9]. | AT-only runs >=20 bp, GC, purine/pyrimidine tracts >=20 bp. These cannot establish a fragile site, replication timing, or origin activity. |
| Genes/transcription context | Transcription and replication can affect instability; genomic position in a gene is context, not a transcription measurement [4,9]. | Optional gene, exon, CDS and >=100-kb gene-body union coverage from GFF3. The 100-kb threshold is a descriptive butterfly screen, not a human fragile-site definition. Expression, R-loops, chromatin and replication timing remain unassessed. |
| Telomeric sequence | Terminal healing/fusion are alternative hypotheses; TTAGG repeats are documented in the butterfly study [1]. | TTAGG/CCTAA arrays >=3 copies. A reference array is not evidence of newly added telomere. Confirm healing with a read-supported novel chromosome-end junction. |
| Assembly/mapping confounders | A coverage boundary may reflect gaps, duplicated sequence, reference mismatch or alignment ambiguity. | Ambiguous bases, optional low-mappability BED, BAM MAPQ bins, source/header/M5 checks, clipping and SA evidence. Softmask is not mappability; no mappability annotation means not assessed. |
| Junction microhomology, blunt joining, insertions | Resolved SV junctions carry repair signatures [2]. Short homology can arise through multiple pathways. | A separate scorer accepts *verified oriented* flanks, insertion and observed alternate sequence, validates their reconciliation, and reports compatible overlap up to 50 bp. It cannot infer NHEJ versus MMEJ/MMBIR uniquely. Templated insertion source requires sequence alignment. |

Most mechanistic evidence comes from human or yeast systems. Conserved molecular processes motivate a screen; they do not justify transferring mammalian effect sizes, PRDM9 motifs, localized-centromere assumptions, or cancer-specific rearrangement models to this butterfly. Evolutionary fusion-boundary enrichment also need not describe a somatic event in one individual.

## Execution on Pegasus

Dependencies: Python 3 with `numpy`, `pysam`, `matplotlib`; SLURM. Use the environment that runs Figure4. No downloads, realignment or new variant calling occur. The first job uses one CPU, 12 GB and the `nano` partition for 30 minutes. Runtime on the real chromosome remains unmeasured; partition/time are configurable.

```bash
cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo
git pull --ff-only
python scripts/ndj_breakpoint/run_breakpoint.py --submit
```

The command prints the job ID and output directory under `GynandroCalling/analysis/NDJ_chr15_breakpoint_context_v1/<input-code-digest>/`. Paths can be changed with `--project`, `--source`, `--figure4` and `--output`. If the saved transition fit is missing, regenerate the existing Figure4 refinement first, or explicitly supply an independently justified **0-based half-open** `--interval START0 END0`. No fallback to a remembered approximate coordinate occurs.

The submitter reconciles the known-source registry and any matching resource-ledger rows with the frozen run, checks file/index existence and BAM contig lengths, and writes `reuse_decisions.tsv` before submission. Missing ledger rows do not mean missing data: the selected existing BAMs are validated directly. It snapshots relevant provenance and freezes worker code. On completion it exports a run-specific resource manifest and registers its artifacts in the project ledger under a file lock, preserving other rows and retaining a pre-update snapshot. BAM chromosome M5, when present, is checked against the fetched reference sequence in the worker. If M5 is absent, reference identity relies on the established config plus accession/length evidence, explicitly not a cryptographic sequence identity proof. Supplementary records themselves are excluded to avoid counting a split alignment twice; SA tags on primary alignments are retained.

Completed runs are immutable and reused only with matching inputs/options/code and output checksums. Incomplete runs must be inspected; specify a fresh output directory to retry. Reference/BAM signatures use resolved path, size and mtime, not expensive whole-file checksums; small configs, annotations and outputs use SHA256.

### Repeat and gene annotations

Supply existing annotations **from the exact Rosina assembly**. This avoids silently applying a TE track from Hmel2.5. The converter does not infer assembly from filenames or download a repeat library.

```bash
python scripts/ndj_breakpoint/prepare_annotations.py \
  --reference /actual/path/GCF_056956725.1_ASM5695672v1_genomic.fna \
  --repeatmasker-out /actual/path/rosina.fa.out \
  --gff /actual/path/rosina_genomic.gff.gz \
  --output /actual/path/NDJ_chr15_annotations
python scripts/ndj_breakpoint/run_breakpoint.py --submit \
  --annotations /actual/path/NDJ_chr15_annotations/annotations.json
```

Each annotation option is optional; placeholders above are intentionally not presumed cluster paths. Existing interval tracks such as low-mappability or segmental-duplication regions can be added with `--bed low_mappability=/actual/path/track.bed`. BED has no header and is 0-based half-open. RepeatMasker and GFF input coordinates are converted from 1-based inclusive. Annotation coverage is unioned, so overlapping records are not double-counted. If an annotation does not exist, create it separately with an assembly-appropriate repeat library/mappability workflow; the current run labels it **not assessed**, never zero.

### Long-homology search

After the first job completes, BLAST+ is an optional second job:

```bash
python scripts/ndj_breakpoint/homology_search.py \
  --run-dir /actual/output/directory --submit
```

Default resources: `cpu`, 2 hours, one CPU, 12 GB; configurable. The script uses local `blastn`, a subject FASTA, masked seeding and overlapping 5-kb queries across the bracket plus 100-kb flanks. It reports nonself hits >=200 aligned bp and >=90% identity, retaining all HSPs permitted by the search. Inspect orientation, position relative to the bracket, repeat family and BAM support together. Query overlap produces duplicate candidate hits; the table is not a count of independent duplicated loci. Masked seeding and these thresholds can miss divergent/short/low-complexity homology. Reuse checks include input/code signatures and output hashes.

### Read evidence and junction follow-up

The first job examines both NDJ BAMs over the bracket plus 100 kb on either side, excluding unmapped, duplicate, QC-fail, secondary and supplementary records. It reports low MAPQ separately; candidate support requires MAPQ >=30 and soft clips >=20 bp (or SA tags). Alignments and unique read-name fragments are distinguished. Exact clipping boundaries are 0-based interbase coordinates; raw SA tags retain SAM's **1-based** positions. SA-only anchor starts are explicitly not called breakpoints. Mate positions are 0-based.

`--control-bam NAME=/path/to/compatible.bam` adds an explicitly selected existing control, after index/reference checks. No control identities or downloaded libraries are invented. Control-discovered SNPs in the earlier scan are not equivalent to control BAM support. Compare controls and both sides before interpreting any clipping pileup; differences may reflect read mixing, depth or haplotype divergence.

Use recurrent clip/SA clusters to choose candidate junctions for IGV review and local assembly/read remapping. Confirm both flanks, orientation, independent-fragment support and reference compatibility. Short reads may not resolve long repeats. The pipeline deliberately does **not** turn the largest clip pileup into an automatic final breakpoint.

For a validated junction, `score_junctions.py --input verified.tsv --output signatures.tsv` expects columns `junction_id`, `support_evidence`, `left_flank`, `right_flank`, `inserted`, `junction_sequence`. Flanks must be oriented 5′→3′ in the alternate molecule and bounded so the observed sequence equals `left_flank + inserted + right_flank[k:]`, where k is supported overlapping homology. Empty insertion is a blank field. It rejects inconsistent sequence and reports search-limit censoring. This is a sequence bookkeeping tool, not an assembler or a repair-pathway classifier.

## Backgrounds, outputs and decision rule

The full uncertainty interval is scored with 0-, 5-kb and 50-kb padding on each side. Each comparison uses up to 999 unique uniformly sampled starts on chromosome 15, with the same interval width, >=90% ACGT and exclusion of the widest target interval plus another 100 kb. Controls can overlap each other: they describe genomic sequence, not replicate biological events. A second panel additionally matches GC within 0.03 and ACGT fraction within 0.02. It does not match on the repeat feature being tested. No significance is reported for GC/callability in their matched panel, or for targets below 90% ACGT or panels with <100 controls. The sampling budget is finite and actual counts are recorded; calipers are never silently relaxed.

For each feature: union fraction, observed score, background median, tie-aware percentile, plus-one upper/lower-tail empirical ranks, two-sided `min(1, 2*min(p_upper,p_lower))`, and exploratory BH values across all tested features/scales/backgrounds. Correlated motifs/scales and overlapping controls make these prioritization measures, not calibrated evidence for a causal mechanism. Position relative to chromosome ends, local recombination and other unmeasured covariates may remain confounded. There is no trained “breakability score.”

Local 1-kb profiles and motif BED intervals locate candidates but have **no sitewise p-values**. This avoids choosing the most interesting motif inside an uncertain bracket and reporting an uncorrected point-breakpoint test. Distances are measured to the bracket, not to an invented exact break.

Main review files:

- `REPORT.md`, `module_status.json`, `control_balance.tsv` — interpretation and missing modules.
- `feature_comparisons.tsv`, `feature_distances.tsv`, `feature_intervals.bed`, `local_1kb_features.tsv` — sequence scores and locations.
- `candidate_reads.tsv`, `clip_clusters.tsv`, `read_evidence_1kb.tsv` — reviewable alignment evidence.
- `breakpoint_context.pdf` — vector reference-rank plots and local profiles, generated on Pegasus.
- `run.json`, `reuse_decisions.tsv`, `resource_ledger_snapshot.json`, `complete.json` — input provenance, completion and hashes.
- `homology/nonself_homology.tsv` — optional candidate partner homology, with its own completion record.

An interesting result would be an unusual feature relative to both backgrounds, localized near an independently supported junction, with compatible partner/orientation evidence and no obvious mapping artifact. A common motif elsewhere in a 20-kb bracket is much weaker evidence. A negative scan does not exclude replication/repair mechanisms, epigenetic fragility or a non-reference repeat insertion.

## Primary sources

1. Höök et al. (2023), *High-density linkage maps and chromosome level genome assemblies unveil direction and frequency of extensive structural rearrangements in wood white butterflies (Leptidea spp.)*. [Chromosome Research](https://doi.org/10.1007/s10577-023-09713-z).
2. Kidd et al. (2010), *A human genome structural variation sequencing resource reveals insights into mutational mechanisms*. [Cell](https://doi.org/10.1016/j.cell.2010.10.027).
3. Porubsky et al. (2022), *Recurrent inversion polymorphisms in humans associate with genetic instability and genomic disorders*. [Cell](https://doi.org/10.1016/j.cell.2022.04.017).
4. Bacolla et al. (2004), *Breakpoints of gross deletions coincide with non-B DNA conformations*. [PNAS / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC521098/).
5. Cer et al. (2012), *Searching for non-B DNA-forming motifs using nBMST*. [Methods paper](https://pubmed.ncbi.nlm.nih.gov/22470144/). Our transparent seed screens are a limited implementation, not nBMST itself.
6. Piazza et al. (2012), *Stimulation of gross chromosomal rearrangements by the human CEB1 and CEB25 minisatellites in Saccharomyces cerevisiae depends on G-quadruplexes or Cdc13*. [Primary experiment](https://pubmed.ncbi.nlm.nih.gov/23133402/).
7. Wang, Christensen and Vasquez (2006), *Z-DNA-forming sequences generate large-scale deletions in mammalian cells*. [Primary experiment](https://pubmed.ncbi.nlm.nih.gov/16473937/).
8. Ozeri-Galai et al. (2011), *Failure of origin activation in response to fork stalling leads to chromosomal instability at fragile sites*. [Primary experiment](https://pubmed.ncbi.nlm.nih.gov/21726815/).
9. Letessier et al. (2011), *Cell-type-specific replication initiation programs set fragility of the FRA3B fragile site*. [Nature](https://doi.org/10.1038/nature09745).

Validation: synthetic motif, coordinate-union, ambiguity, matched-background, empirical-tail/BH, verified-junction and synthetic indexed-BAM integration tests. Actual Pegasus NDJ input execution, full-chromosome runtime, BLAST execution and biological plot review are not yet performed.
