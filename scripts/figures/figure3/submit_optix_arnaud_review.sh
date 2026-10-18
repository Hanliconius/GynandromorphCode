#!/usr/bin/env bash
set -euo pipefail
repo=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo
run=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc
stage="$run/arnaud_review_code_b69bc47"
out="$run/Optix_Arnaud_review_v1"
if [ -e "$out" ]; then
  printf 'Output already exists; inspect before resubmitting: %s\n' "$out" >&2
  exit 1
fi
mkdir -p "$stage" "$run/logs"
git -C "$repo" archive b69bc474ffa160cdce7451f21f23f178e62a86ef \
 scripts/figures/figure3/optix_arnaud_review.py \
 scripts/figures/figure3/plot_optix_uniform_genoplots.py \
 scripts/figures/figure3/plot_optix_merged_race_frequencies.py \
 scripts/figures/figure3/plot_optix_frequency_genoplots.py \
 scripts/figures/figure3/plot_optix_all_phenotype_controls.py \
 scripts/figures/figure3/pk_dark_pale_samples.tsv | tar -x -C "$stage"
job=$(sbatch --parsable --partition=nano --cpus-per-task=2 --mem=12G --time=00:30:00 \
 --job-name=Optix_review --output="$run/logs/Optix_review.%j.out" --error="$run/logs/Optix_review.%j.err" \
 --wrap="micromamba run -n gynandrocalling python '$stage/scripts/figures/figure3/optix_arnaud_review.py' --run '$run' --pk-catalog '$stage/scripts/figures/figure3/pk_dark_pale_samples.tsv' --output '$out'")
printf 'JOBID\t%s\nMAIN_PDF\t%s\nREVIEW_DIR\t%s\n' "$job" "$out/merged/Dennis_main_candidate_BA90.pdf" "$out"
