# Compute and plotting workflow

For the GynandromorphCode project, use this execution workflow for analyses and figures:

1. Write or edit analysis and plotting code in the development workspace.
2. Commit and push the code to GitHub.
3. Pull the committed code on Pegasus and run analyses, large jobs, and all data-driven plot generation there, using the appropriate environment and Slurm resources.
4. Transfer the requested plots, tables, logs, or other review outputs back from Pegasus so they can be inspected and revised here.

Precomputed outputs such as bigWigs, BAMs, coverage tables, and variant calls are useful inputs, not limits on what can be reproduced. When faithful replication requires it, rerun the upstream steps from raw reads or earlier inputs on Pegasus, reusing the existing alignment and calling array workflows where appropriate. Keep the input source, reference assembly, filters, and processing steps explicit in the code and audit outputs.

Do not run computationally expensive analyses or data-driven plotting in the development workspace. The development workspace is for writing, reviewing, and validating code; Pegasus is the compute and figure-generation environment. Do not describe a plot as generated or reviewed until the actual Pegasus output has been transferred and inspected.

Keep user-facing shell commands concise. Redirect verbose program output to log files and report only useful status, output paths, and short summaries.


# Reuse existing data before new work

Before any new download, alignment, calling, or coverage job, follow `docs/DATA_RESOURCE_MANIFESTS.md` and reconcile the requested samples against `manifests/known_analysis_sources.tsv`, existing run/source manifests, and the Pegasus resource ledger. Trace source paths and aliases, including paths outside the project root. Verify reference, genomic span, statistic/filters, file existence, and completion evidence. A short exported VCF does not establish the coverage span of its source BAM; a filename search that finds nothing does not establish data absence. Capture/targeted libraries require interval checks rather than automatic exclusion.

Write a per-sample/per-interval reuse decision before submission: reuse calls, recount an existing alignment, align existing reads, download specifically missing mates, or unresolved. Unresolved inventory is not authorization to download automatically. Record counts and reasons for newly required work. Keep the resolved Pegasus ledger current after jobs finish, and commit small provenance registries/audit summaries to Git without sequence data. Preserve the largest eligible cohort for each requested comparison; do not reduce a cohort merely to make panel totals match unless explicitly requested.
