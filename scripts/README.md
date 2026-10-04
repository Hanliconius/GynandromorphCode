# Scripts

Scripts are organized by scientific task rather than historical analysis number.

Planned layout:

scripts/
  allele_frequency/
  founder_models/
  reciprocal_mixing/
  coverage/
  optix/
  smudgeplot/

See ../docs/PEGASUS_IMPORT_MANIFEST.md for the historical-to-cleaned mapping.

## Conventions

Cleaned scripts should:

- accept inputs and output roots as arguments;
- avoid hard-coded Pegasus paths;
- avoid silently recomputing expensive intermediates;
- remove partial outputs before a clean rerun when needed;
- use descriptive output names;
- write vector PDF/SVG where practical;
- record critical parameters in stdout/logs;
- preserve biological defaults explicitly.

For Slurm:

- use project-relative paths where possible;
- request explicit resources;
- use one log pair per job or array task;
- give FastK concurrent jobs unique temporary directories;
- prefer direct environment binary paths when concurrent micromamba calls would contend for locks.
