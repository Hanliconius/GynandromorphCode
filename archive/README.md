# Archive policy

Archive directories preserve provenance and are not imported by maintained workflows.

## `pegasus_shell_raw/`

Exact shell and Slurm history from Pegasus, including retries, repair jobs, resource experiments, and hard-coded cluster paths. Its manifest is the authoritative inventory.

## `AF_sitecount_provenance_raw/`

Historical Python analyses used to reconstruct AF inclusion, founder models, simulations, and narrative figures. These files are particularly valuable for explaining site counts and plotting decisions, but are not canonical until promoted.

## Rules

- Do not edit archival files to make them runnable.
- Do not call archival scripts from maintained workflows.
- Promote by copying and cleaning one selected implementation into `scripts/<lane>/`.
- Record the archival source path and major changes.
- Keep failed or superseded variants only here.
