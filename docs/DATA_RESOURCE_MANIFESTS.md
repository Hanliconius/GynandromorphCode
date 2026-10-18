# Data resource manifests

Before submitting downloads, alignment, calling, or coverage recomputation, reconcile the proposed samples against existing manifests and source configurations. A filename search is supplementary evidence, not proof of absence.

## Source registry

`manifests/known_analysis_sources.tsv` identifies established analysis/configuration sources and the limits of the evidence supporting them. Entries marked `documented_not_revalidated` must be checked on Pegasus. Neither documentation nor a completion file alone establishes that all underlying inputs still exist.

On Pegasus maintain the resolved resource ledger at:

`/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/data/derived/resource_manifest.tsv`

Keep timestamped ledger snapshots and run-specific reuse decisions in each analysis output directory. Commit small source registries, schemas, and compact audit summaries to Git; do not commit sequence data.

## Resource ledger: one row per artifact

Required fields:

- `resource_id`, `biological_sample_id`, `run_accessions`, `sample_aliases`, `library_id`, `role`
- `resource_type`: FASTQ mate, BAM/CRAM, VCF/BCF, allele counts, coverage profile, bigWig, reference
- `path`, `resolved_path`, `exists`, `bytes`, `mtime_utc`, `checksum`, `checksum_type`, `checked_at_utc`
- `source_config`, `source_config_field`, `producing_code_commit`, `completion_marker`, `validation_status`
- `reference_accession`, `contig`, `start_1based`, `end_1based`, `interval_coverage_evidence`
- `sequencing_type`, `sequencing_type_source`, `filters`, `statistic`, `read_group_samples`, `indexed`

Unknown values remain explicit. Resource types are not interchangeable: depth profiles establish coverage availability, not allele counts; an exported short VCF does not establish the span of its source BAM. Capture/targeted libraries are not excluded merely because they are not WGS. Inspect target definitions and observed interval coverage.

## Before new work: one decision per sample and interval

Write `reuse_decisions.tsv` with:

`biological_sample_id, interval, reference, action, selected_resource_id, rejected_resource_ids, evidence, unresolved_questions`

Actions:

1. `reuse_calls`: compatible allele-level calls/counts already exist.
2. `recount_existing_alignment`: compatible BAM/CRAM exists but needed calls do not.
3. `align_existing_reads`: usable read pairs exist; an appropriate alignment does not.
4. `download_missing_reads`: manifests and their referenced paths have been checked and exact missing mates identified.
5. `unresolved`: identity, location, reference, sequencing coverage, or processing compatibility remains uncertain. Do not turn this into an automatic download.

Read manifests recursively: follow `sources`, task/sample BAM and FASTQ fields, symlink targets, and paths outside the project root. Retain the provenance chain. Deduplicate aliases and sequencing runs belonging to the same biological individual; preserve focal libraries and pools as distinct resource roles.

Each submission records reused/new counts, exact mate downloads, interval scope, and resource compatibility evidence. Existing completed work is immutable. Completion updates the ledger with actual outputs and failures. Refresh the ledger before a new run, and revalidate selected files immediately before use.

## 2026-10-18 Yb correction

The broad dorsal VCF contains 502 header samples; 501 had at least one callable dorsal record in the audit. The current eligible non-PK dorsal display includes 473 individuals. This is not the cohort size of the earlier two-interval coverage workflow.

`docs/YB_ALL_SAMPLES_COVERAGE.md` documents profiles spanning both d/v and reusing earlier profiles/BAMs, with no downloads or alignments. Its earlier count was 39 libraries/individuals per reference. `scripts/yb/submit_Yb_luca_median_replication.py` uses the retained version with 38 per reference: 20 PK, 14 comparative controls, and four focal libraries. It traces source configurations and BAM paths explicitly.

`docs/YB_COVERAGE.md` also documents the existing 18-control HmRos BAM panel and input manifest. These are manifest evidence that must be reconciled, not evidence that all 473 dorsal controls have matching local ventral data.

The submitted 448-task ventral expansion was prepared before this reconciliation was complete. Its task count does not establish 448 necessary downloads or alignments. Existing downloads and completed outputs should be retained and accounted for; pending work needs reconciliation against the ledger before continuation.
