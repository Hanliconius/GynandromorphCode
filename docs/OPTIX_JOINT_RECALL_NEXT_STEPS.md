# Optix/Yb joint recall — next actions after successful 2026-10-14 preflight

## Verified preflight
Run dir: `/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/optix_joint_recall_preflight/20261014T193643Z`

- 38 libraries = 34 comparative specimens + four focal tissue libraries.
- Old 34 comparative = 20 retained PK (excluding P4009) + 10 original controls + four additional controls.
- 5,874 lines containing cythera/vulcanus in older metadata; NOT the number of individuals, sequencing runs or independent accessions.
- The corrected PK run's Rosina FASTA and the old frozen Optix REF alleles were validated; no re-download of the reference is necessary.
- Six historical BioSamples have *not yet* been resolved to paired WGS reads by the successful preflight.

## Resolve specimens & verify current ENA runs

The read-only `scripts/figures/figure3/optix_joint_candidate_review.py` first extracts and ranks distinct accession hints in analysis-142 metadata, then optionally contacts ENA in bounded parallel batches. It validates exactly two mate FASTQs, MD5s, sizes, ENA sample accession, sequencing strategy WGS, and layout PAIRED. It does not download FASTQs or submit jobs.

```bash
repo=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo
root=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
cd "$repo"
git fetch origin optix-joint-recall-20261014
git restore --source=FETCH_HEAD --worktree -- scripts/figures/figure3/optix_joint_candidate_review.py
micromamba run -n gynandrocalling python \
  scripts/figures/figure3/optix_joint_candidate_review.py \
  --preflight "$root/analysis/optix_joint_recall_preflight/20261014T193643Z" \
  --top-per-race 12 --ena-workers 4 --ena
```

Outputs, under a new timestamped `analysis/optix_joint_candidate_review` directory:
- `original142_accession_hints.tsv`, `original142_matching_files.tsv`: hints and per-file counts, not a callset.
- `all_resolved_runs.tsv`, `verified_unique_wgs_candidates.tsv`: ENA-validated runs and one representative sequencing run per BioSample.
- `candidate_races_for_preflight.tsv`: race, BioSample, label and selected run in the input format for the earlier preflight's `--new-races`.
- `ena_unresolved.tsv`: missing or incompatible ENA records.
- `tools_on_login.json`: whether GATK, Java, samtools, bwa, bcftools etc. are visible on the login node.

**Do not start bulk FASTQ downloads or genotyping until the new race candidate count and reference/tool provenance are known.**

## Next implementation step
If 3–5 suitable independent cythera and vulcanus specimens are confirmed, freeze a single `manifest.tsv` with read URLs, expected sizes and checksums, plus the 38 existing libraries, and build bounded arrays for (1) download/MD5, (2) mapping to the identical Rosina reference and duplicate marking, (3) diploid comparative sample GVCFs in Optix/Yb interval shards, (4) separate per-interval joint GATK genotyping, (5) final genotype/GQ/callability audit, (6) focal allele depth and reciprocal mixing posterior, (7) regenerated vector PDFs.

If the WGS archive only has two confirmed cythera individuals, explicitly report that shortfall; do not count additional lanes or 9 capture-only individuals as additional WGS specimens. Decide how to annotate capture data separately if it is included in downstream comparisons.

The old mixed legacy/new SNP table is retained only for comparison and must not enter the new joint callset.
