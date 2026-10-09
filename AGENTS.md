# Compute and plotting workflow

For the GynandromorphCode project, use this execution workflow for analyses and figures:

1. Write or edit analysis and plotting code in the development workspace.
2. Commit and push the code to GitHub.
3. Pull the committed code on Pegasus and run analyses, large jobs, and all data-driven plot generation there, using the appropriate environment and Slurm resources.
4. Transfer the requested plots, tables, logs, or other review outputs back from Pegasus so they can be inspected and revised here.

Precomputed outputs such as bigWigs, BAMs, coverage tables, and variant calls are useful inputs, not limits on what can be reproduced. When faithful replication requires it, rerun the upstream steps from raw reads or earlier inputs on Pegasus, reusing the existing alignment and calling array workflows where appropriate. Keep the input source, reference assembly, filters, and processing steps explicit in the code and audit outputs.

Do not run computationally expensive analyses or data-driven plotting in the development workspace. The development workspace is for writing, reviewing, and validating code; Pegasus is the compute and figure-generation environment. Do not describe a plot as generated or reviewed until the actual Pegasus output has been transferred and inspected.

Keep user-facing shell commands concise. Redirect verbose program output to log files and report only useful status, output paths, and short summaries.
