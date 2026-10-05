# Pegasus import manifest

Historical working root:

`/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling`

Two historical sources are intentionally kept outside the maintained script tree.

| Provenance source | Location | Purpose |
|---|---|---|
| shell/Slurm intake | branch `pegasus-shell-import-2026-10-05`, `archive/pegasus_shell_raw/` | exact runnable history, including retries and cluster-specific paths |
| AF/site-count provenance | `archive/AF_sitecount_provenance_raw/` on the import branch | Python history for AF, founder-state, simulation, and narrative figures |

The shell intake contains 176 shell/Slurm files with hashes recorded in its manifest. These files are evidence of what ran, not the maintained interface.

## General import policy

For each promoted script:

1. record the original Pegasus path in the header;
2. preserve scientific defaults;
3. remove hard-coded absolute paths where practical;
4. convert opaque positional assumptions to named arguments;
5. replace analysis-number filenames with descriptive names;
6. make reruns idempotent;
7. preserve or improve audit logging;
8. verify outputs against historical results before retiring the original.

Do not commit FASTQs, BAMs, VCFs, FastK tables, `.smu` files, or large generated figures.

## Allele-frequency figures

Historical candidates include:

- `plot_TN_AF_main_figure.py`
- `plot_founder_slides_empirical_depth.py`
- `plot_supplementary_haplotype_scenarios_50pct_rows.py`
- `plot_slide_AF_original_GT_loess_left_overlay_ERR2298208.py`
- `plot_slide_AF_leftbars_ERR2298208loess_brokenaxis.py`
- `plot_depth_dependent_AF_dispersion_empirical.py`

Suggested destination: `scripts/allele_frequency/`.

The promoted workflow must also expose the site-count audit explaining how source-VCF membership, genotype calls, missingness, and plot filters produce the final `00-left/00-right` counts.

## Founder-model fitting

Historical candidates include:

- `scripts/158_audit_three_vs_four_founder_inputs.py`
- `scripts/159_fit_autosomal_three_vs_four_founder.py`
- `scripts/160_simulate_calibrate_autosomal_founder_models.py`
- later sequential founder-distribution and narrative scripts, approximately `202–216`

Suggested destination: `scripts/founder_models/`.

## Reciprocal-mixture genotype display

Historical reciprocal-mixture inputs live under `analysis/original_502/reciprocal_L_R_mixture/`, with a genome-wide numbered workflow around `173–181`.

Suggested destination: `scripts/reciprocal_mixing/`.

Critical rule: focal Left/Right displayed genotype pairs require posterior support of at least 0.95, with no fallback.

## Optix

The optix phenotype analyses include numbered scripts around analyses `116`, `142`, and `143`, plus later comparison-panel builders.

Suggested destination: `scripts/optix/`.

Comparison samples must be chosen independently of the illustrated L/R result.

## Coverage

Historical candidate: `plot_relative_coverage_bar_panels.py`.

Suggested destination: `scripts/coverage/plot_relative_coverage.py`.

The Z panel remains excluded until sex/ploidy normalization and interpretation are independently validated.

## Smudgeplot

Historical/current candidates include:

- `scripts/plot_smudge_smu_luka_bleute.py`
- `scripts/plot_smudge_smu_viridis.py`
- `scripts/prepare_smudgeplot_arrangements.sh`
- `scripts/run_smudgeplot_arrangements.sh`
- `scripts/rerun_smudge_control_and_4legs.slurm`
- `scripts/plot_smudge_sixpanel_viridis.py`
- `scripts/plot_smudge_sixpanel_viridis_y120.py`
- `scripts/build_smudge_full_detail_bundle.slurm`
- later lower-L sweep scripts that generated `analysis/smudgeplot_lowerL_20261002/L6_download/`

Suggested destination: `scripts/smudgeplot/`. The existing portable renderer is the active seed for this lane.

## Nondisjunction

Historical candidates include the L/R difference-scan chain and later `217*` and `224*` analyses.

Suggested destination: `scripts/ndj/`.

This lane remains exploratory until the phenotype, null model, genome-wide denominator, and sensitivity analyses are documented.

## Useful historical caches and audits

Examples include:

- `ERR2298208_vs_TNLeft_original_GT_counts.npz`
- `control_benchmarks/best_coverage_control_self_hex45_batlowW.AF.audit.json`
- `TN_merged_control_sites_hex40_batlowW.AF.audit.json`

These normally remain outside git if large; retain checksums and compact manifests instead.

## Current migration priority

1. AF inclusion/site-count audit;
2. founder-model fit and simulation calibration;
3. reciprocal-mixing AF workflow;
4. coverage workflow;
5. Smudgeplot lower-L sweep and final comparison plotting;
6. optix comparison workflow;
7. NDJ exploratory workflow.

Promotion details and cross-cutting scientific rules are recorded in [CURATION_INDEX.md](CURATION_INDEX.md). The maintained directory contract is in [scripts/README.md](../scripts/README.md).
