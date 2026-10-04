# Pegasus import manifest

Historical working root:

/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

The goal is to migrate the scientifically important scripts into descriptive repository paths while keeping exact historical provenance.

## General import policy

For each script:

1. record the original Pegasus path in the header;
2. preserve scientific defaults;
3. remove hard-coded absolute paths where practical;
4. convert opaque positional assumptions to named arguments;
5. replace analysis-number filenames with descriptive names;
6. make reruns idempotent;
7. preserve or improve audit logging;
8. verify outputs against historical results before retiring the original.

Do not commit FASTQs, BAMs, VCFs, FastK tables, .smu files, or large generated figures.

## Allele-frequency figures

Historical scripts to recover include:

plot_TN_AF_main_figure.py
plot_founder_slides_empirical_depth.py
plot_supplementary_haplotype_scenarios_50pct_rows.py
plot_slide_AF_original_GT_loess_left_overlay_ERR2298208.py
plot_slide_AF_leftbars_ERR2298208loess_brokenaxis.py
plot_depth_dependent_AF_dispersion_empirical.py

Suggested destination:

scripts/allele_frequency/

## Founder-model fitting

Historical scripts:

scripts/158_audit_three_vs_four_founder_inputs.py
scripts/159_fit_autosomal_three_vs_four_founder.py
scripts/160_simulate_calibrate_autosomal_founder_models.py

Later sequential founder-distribution scripts and narrative builders should also be imported.

Suggested destination:

scripts/founder_models/

## Reciprocal-mixture genotype display

Historical reciprocal-mixture inputs live under:

analysis/original_502/reciprocal_L_R_mixture/

The optix phenotype analyses include numbered scripts around analyses 116, 142, and 143.

Suggested destinations:

scripts/reciprocal_mixing/
scripts/optix/

Critical rule: focal Left/Right displayed genotype pairs require posterior >= 0.95 with no fallback.

## Coverage

Historical:

plot_relative_coverage_bar_panels.py

Suggested destination:

scripts/coverage/plot_relative_coverage.py

## Smudgeplot

Historical / current scripts to recover:

scripts/plot_smudge_smu_luka_bleute.py
scripts/plot_smudge_smu_viridis.py
scripts/prepare_smudgeplot_arrangements.sh
scripts/run_smudgeplot_arrangements.sh
scripts/rerun_smudge_control_and_4legs.slurm
scripts/plot_smudge_sixpanel_viridis.py
scripts/plot_smudge_sixpanel_viridis_y120.py
scripts/build_smudge_full_detail_bundle.slurm

Also recover the later lower-L sweep scripts that generated:

analysis/smudgeplot_lowerL_20261002/L6_download/

Suggested destination:

scripts/smudgeplot/

## Useful historical caches / audits

Examples:

ERR2298208_vs_TNLeft_original_GT_counts.npz
control_benchmarks/best_coverage_control_self_hex45_batlowW.AF.audit.json
TN_merged_control_sites_hex40_batlowW.AF.audit.json

These should normally remain outside git if large; retain checksums/manifests instead.

## Next import priority

1. founder-model fit + simulation-calibration scripts;
2. reciprocal-mixing AF main-figure scripts;
3. coverage script;
4. Smudgeplot lower-L sweep and final comparison plotting;
5. optix reciprocal-mixture display scripts.
