# Pegasus import manifest

Two historical sources are intentionally kept outside the maintained script tree.

| Provenance source | Location | Purpose |
|---|---|---|
| shell/Slurm intake | branch `pegasus-shell-import-2026-10-05`, `archive/pegasus_shell_raw/` | exact runnable history, including retries and cluster-specific paths |
| AF/site-count provenance | `archive/AF_sitecount_provenance_raw/` on the import branch | Python analysis history for AF, founder-state, simulation, and narrative figures |

The shell intake contains 176 shell/Slurm files with hashes recorded in its manifest. These files are evidence of what ran, not the maintained interface.

## Current triage

- Allele-frequency and founder-model candidates have been identified.
- Reciprocal-mixing scripts are grouped as one workflow rather than promoted as scattered numbered jobs.
- Coverage, optix, smudgeplot, and NDJ have dedicated lanes.
- Posterior-supported displays retain the agreed `posterior >= 0.95` rule with no fallback substitution.
- The Z-coverage panel remains excluded until its interpretation and normalization are validated.
- Retry suffixes and numbered variants remain provenance until a canonical version is selected.

Promotion details are recorded in [CURATION_INDEX.md](CURATION_INDEX.md), and the maintained directory contract is in [scripts/README.md](../scripts/README.md).
