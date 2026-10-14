# Optix/Yb 2026-10-14: calling-method correction

**Scientific requirement:** SNPs must be called consistently across the full
comparative cohort using the established project method, rather than combining
an old VCF with newly called individuals.

## Approved prior method

The earlier GynandromorphCode implementation is in
`scripts/figures/figure3/extend_PK_optix_Yb.py`. Its genotype calls are made
using `bcftools mpileup` with `-a FORMAT/DP,FORMAT/AD -q 30 -Q 20`
followed by `bcftools call -m -A -f GQ`. The previously constrained Optix
per-site genotyping also used `-C alleles -T <allele-list>`; that constrained
mode **must not** be mistaken for novel SNP discovery. The new cohort
joint call discovers sites across all retained comparative BAMs without
importing frozen historical alleles, while retaining the same bcftools caller
and MAPQ/base-quality settings.

Preprocessing also follows the established BWA-MEM `-M` full-reference
alignment and Picard MarkDuplicates workflow.

## Incorrect, disabled alternative

`scripts/figures/figure3/optix_joint_recall.py` initially submitted an
unrequested GATK HaplotypeCaller / CombineGVCFs / GenotypeGVCFs pipeline.
This is a **different caller** and was not explicitly approved. The GitHub
submission entrypoint is disabled. The initial setup and BAM alignment jobs
are useful, but the HaplotypeCaller array `138468793` and dependent arrays
`138468794`, `138468795`, and `138468796` must be cancelled on Pegasus.
No GATK calls should be used in any analysis.

## Recovery

Existing raw-read array: `138468791`.
Existing BWA/Picard alignment array: `138468792`.
Original Rosina run:
`/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/optix_joint_all_samples/20261014T211057Z`.

Known-failed and previously problematic P4003 (`SRR11243351`,
comparative sample array index 14) is excluded. The recovered dataset is
43 independent comparative individuals (including 5 cythera and 5 vulcanus)
and 4 focal tissue libraries. The existing successful 47 alignments are
retained. No one should rerun or overwrite a successful alignment to address
the excluded specimen.

Recovery script:
`scripts/figures/figure3/resume_bcftools_optix_joint.py`.
It submits two independent region-wise, 43-sample bcftools joint calls
with dependencies on the **individual 47 selected alignment array tasks**
(`afterok:138468792_0:..._47` excluding task 14). It then uses the
pre-existing focal pileup and QC functions from the frozen, original
aligned-run script. These functions do not perform GATK calling.
New outputs are isolated in
`20261014T211057Z/bcftools_joint_recall/`; existing BAMs are accessed
read-only. Each joint VCF records its bcftools version and exact options.

**Important distinction:** A simultaneous multi-sample bcftools call is
not equivalent to combining independently called VCFs or forcing genotypes
at frozen legacy sites. It keeps the established calling engine and filters
while satisfying the new requirement for uniform cohort-wide discovery.

## What has not been done

This branch change does not SSH to Pegasus, cancel jobs, or submit new jobs.
The user must perform the cancellation and corrected submission commands.
Phenotype-informed SNP selection, reciprocal-mixing posteriors, and revised
figure PDFs remain downstream analysis work, not automatic products of VCF
generation.
