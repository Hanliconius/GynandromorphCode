# Maintained scripts

This directory is the curated execution layer for the gynandromorph project. Historical Pegasus scripts remain on branch `pegasus-shell-import-2026-10-05` under `archive/pegasus_shell_raw/`; AF/founder-model provenance is retained separately under `archive/AF_sitecount_provenance_raw/`.

## Workflow lanes

| Directory | Scope | Status |
|---|---|---|
| `allele_frequency/` | observed AF summaries, original-GT histograms, site-count audits, depth-aware AF displays | candidates selected; promotion pending |
| `founder_models/` | founder-state fitting, simulation, posterior summaries, narrative figures | candidates selected; promotion pending |
| `reciprocal_mixing/` | reciprocal L/R mixing inference and posterior-supported pair calls | candidates selected; promotion pending |
| `coverage/` | relative coverage summaries and chromosome-level panels | candidates selected; Z held out pending validation |
| `optix/` | optix-region extraction, comparative panels, and unbiased sample selection | candidates selected; promotion pending |
| `smudgeplot/` | k-mer-pair diagnostics and figure rendering | renderer active; upstream jobs pending cleanup |
| `ndj/` | nondisjunction and L/R difference scans | exploratory; promotion pending |

See [the curation index](../docs/CURATION_INDEX.md) for source candidates and scientific rules.

## Promotion standard

A promoted script must have explicit inputs and outputs, configurable project paths, environment/resource notes, a validation command, and a clear relationship to the narrative methods. Numbered retries and cluster repairs remain archival unless they encode a distinct method.

No historical directory is merged wholesale into `scripts/`. The goal is one maintained implementation per step, with source provenance documented.
