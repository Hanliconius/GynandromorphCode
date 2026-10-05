# Script curation index

This index maps the raw Pegasus history to the maintained narrative workflow. Listed filenames are candidate sources and must be cleaned and validated before promotion.

| Lane | Candidate source material | Canonical intent | Status |
|---|---|---|---|
| Allele frequency | `plot_TN_AF_main_figure.py`, original-GT/LOESS plots, depth-dispersion scripts, AF site-count audits | transparent definition of included sites and L/R AF panels | promote first |
| Founder models | `158_audit*`, `159_fit*`, `160_simulate*`, `202–216` model/narrative series | fitted founder states, simulation checks, and posterior summaries | promote first |
| Reciprocal mixing | genome-wide `173–181` series and posterior-95 pair tables | infer reciprocal mixing while preserving sample identity and uncertainty | promote first |
| Coverage | `plot_relative_coverage_bar_panels.py` and associated summaries | normalized chromosome coverage panels | review; Z excluded |
| Optix | `116*`, `142*`, `143*`, and later comparative-panel scripts | optix comparison selected independently of L/R outcome | review |
| Smudgeplot | preparation/rerun/bundle Slurm scripts plus maintained renderer | reproducible k-mer-pair diagnostic figures | renderer active |
| NDJ | L/R difference scan chain and later `217*`/`224*` analyses | nondisjunction hypotheses and genome-wide L/R diagnostics | exploratory |

## Non-negotiable scientific rules

- Distinguish observed genotype calls, inferred founder states, and simulated values.
- Report exactly why sites enter an AF plot; do not treat implicit VCF inclusion as a biological filter.
- Preserve left/right labels and original sample IDs through every join.
- Display posterior-supported pair calls only at `posterior >= 0.95`; do not silently substitute a lower-confidence model.
- Select comparative optix samples independently of the L/R result being illustrated.
- Keep Z coverage out of the canonical coverage figure until separately validated.
- Record empirical depth handling and missingness filters.

## Version-selection rule

Start from the latest script known to have produced the accepted output. Compare its logic with earlier numbered versions, retain only scientifically meaningful changes, and replace embedded Pegasus paths with arguments or a small configuration file. The archived filename remains recorded in the promoted script header and migration notes.

## Validation expected for every lane

- sample and site counts at each filter;
- reference/coordinate system;
- input hashes or immutable source paths where practical;
- expected tables/figures;
- smoke-test command;
- statement of active, exploratory, superseded, or archival status.
