# Genome-wide Rosina SNP recall: single shared source for Optix and Yb

**Reference:** Rosina assembly `GCF_056956725.1`. No new SNP calls
on Hmel2.5. This replaces the cancelled locus-only reconstruction.

**Source reads and BAMs:** the existing
`GynandroCalling/analysis/optix_joint_all_samples/20261014T211057Z`
download and BWA/Picard alignment arrays are retained. The 47 eligible
libraries include 43 independent comparative individuals and four focal
tissues. P4003 (`SRR11243351`, original array element 14) is explicitly
excluded. The 43 comparisons include five cythera and five vulcanus.

## Calling method and array design

The executable is
`scripts/figures/figure3/rosina_whole_genome_bcftools.py`.

- Plan **every contig** directly from the frozen Rosina FASTA `.fai`.
- Default 2-Mb genomic-core window arrays, with 200 bp of context on
  either side and strict `POS` cropping back to each non-overlapping
  core; no sharding by individuals and no `%N` Slurm throttles.
- Every window reads the **same list of 43 comparison BAMs** and uses
  `bcftools mpileup -Ou -a FORMAT/DP,FORMAT/AD -q 30 -Q 20`,
  then `bcftools call -m -A -v -f GQ`.
  The historical caller and quality filters are retained. The added
  `-v` requests **variant-only emission** to avoid writing hundreds
  of millions of reference-only genotypes genome-wide; no legacy
  VCF or frozen allele list is imported.
- Emit each window as BCF with a CSI index; merge in reference contig
  order, then create one indexed **genome-wide multi-sample VCF** and
  a derived indexed **genome-wide SNP-only VCF**.
- No GATK, no GVCFs, no Hmel2.5 genotype calls; transposable elements
  are **not masked or removed** from the BAMs. As in prior calls, MAPQ30
  can exclude multimapping TE reads from called SNP evidence.
- Sex chromosomes, W-linked regions and mitochondria are emitted
  with the caller's default diploid assumptions. Those genotypes
  must be interpreted/filtered after sample-sex and ploidy QC, not
  automatically accepted as diploid biological truth.

## Existing BAM reuse: evidence required, not merely BAM presence

The earlier `audit_reusable_rosina_bams.py` mistakenly inspected the
`fasta` key when previous Yb frozen configurations normally use
`references[].source`. The updated audit now checks the correct key,
matching the reference FASTA source, matching full reference dictionary,
sample read group, BAM index, mapped reads on at least ten long
contigs, and duplicate-marking provenance through the previously
generated Picard metrics/complete markers or BAM `@PG`.

A candidate is chosen for reuse only if the updated report has all of:
`status=POTENTIAL_REUSE`,
`source_ref_explicit_match=YES`,
`long_contigs_mapped>=10`, and
`picard_metrics_present=YES` or `has_picard_pg=YES`.
Candidates failing any of these checks use the new Oct-14 alignment.
The 23 earlier candidates are **not yet verified by the corrected
audit**, so the actual number adopted may be lower than 23.

The new script never reruns or cancels original downloads or BAM
alignments. It makes read-only symlinks to approved old or new BAMs.
It prints the **original alignment array indices that would be
redundant** if old BAM reuse is approved; cancellation is a separate
manual decision and must not precede source verification.

## Slurm stages

1. `validate`: all 47 BAMs present, indexed, sample-ID read group,
   identical full Rosina dictionary and genome-wide mapping; checks
   original FAI digest.
2. `call`: parallel 2-Mb windows for every Rosina contig,
   all 43 comparison BAMs per task.
3. `chrom`: independently concat windows per reference contig.
4. `combine`: concat contigs in FASTA order, generate complete and
   SNP-only genome-wide VCFs with CSI indexes.
5. `extract`: Yb and Optix extract stages with NO new SNP calling.
6. `focal`: separately count allele depths in T2L, T3L, T2R, T3R
   and pool by side after counting.
7. `publish`: bcftools stats and a completion JSON; only then
   updates `analysis/rosina_genomewide_joint_recall/latest_complete.path`.

Original alignment tasks are stage dependencies **only for samples
that lack an approved pre-existing BAM**. Failed original element 14
is never included. Batch sizes limit the number of entries per Slurm
array, **not concurrent tasks**.

## Shared output contract for the Yb and Optix analyses

Under `analysis/rosina_genomewide_joint_recall/<run>/`:

- `joint/Rosina.genomewide.43comparative.variants.vcf.gz` + `.csi`
- `joint/Rosina.genomewide.43comparative.snps.vcf.gz` + `.csi`
- `bam_manifest.tsv`: all 47 sample IDs, groups, roles, tissue sides
  and the frozen path to each actual BAM.
- `integration_Yb/Yb_chromosome.vcf.gz` + `.csi`
- `integration_Yb/Yb_window.vcf.gz` + `.csi`
- `integration_Yb/Yb_focal_Left_Right_allele_counts.tsv.gz`
- `integration_Yb/rosina_bam_manifest.tsv`
- `integration_Yb/complete.json` describing Rosina reference,
  SNP source VCF, extracted regions, BAM provenance, and counts.
- Parallel `integration_Optix/` products for the Optix workflow.

Yb analysis scripts can use
`scripts/yb/resolve_shared_rosina_callset.py` to locate the
shared source and manifest; this helper does not call SNPs and does
not change existing Hmel2.5 coverage plots.

## Safe preflight and launch

On Pegasus, from the GynandromorphCode repo:

```bash
root=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
source_run="$root/analysis/optix_joint_all_samples/20261014T211057Z"
micromamba run -n gynandrocalling python \
  scripts/figures/figure3/audit_reusable_rosina_bams.py --run "$source_run"

audit="$source_run/existing_bam_reuse_audit/candidate_audit.tsv"
micromamba run -n gynandrocalling python -m unittest discover \
  -s tests -p 'test_rosina_whole_genome_bcftools.py' -v

micromamba run -n gynandrocalling python \
  scripts/figures/figure3/rosina_whole_genome_bcftools.py \
  --preview --source-run "$source_run" --reuse-audit "$audit"
```

Check that the preview reports the **whole-genome** contig count,
genome size, window count, actual eligible reuse count and remaining
alignment tasks. Only then invoke the same script with `--submit`.
Submitting queues downstream arrays; it does not relaunch downloads
or existing alignments. Recheck dependencies and BAM provenance before
selectively cancelling any redundant original alignment array element.
