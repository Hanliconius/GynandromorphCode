#!/usr/bin/env bash
set -euo pipefail
repo=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo
run=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc
stage="$run/arnaud_nonPK_code_079e661"
out="$run/Optix_Arnaud_review_nonPK_v2"
if [ -e "$out" ]; then
  printf 'Output already exists; inspect before resubmitting: %s\n' "$out" >&2
  exit 1
fi
mkdir -p "$stage" "$run/logs"
git -C "$repo" archive 079e661382271c2dcf42e8aaebe838271db953ce \
 scripts/figures/figure3/optix_arnaud_review_nonPK.py \
 scripts/figures/figure3/plot_optix_uniform_genoplots.py \
 scripts/figures/figure3/plot_optix_merged_race_frequencies.py \
 scripts/figures/figure3/plot_optix_frequency_genoplots.py \
 scripts/figures/figure3/plot_optix_publication_nonPK.py \
 scripts/figures/figure3/pk_dark_pale_samples.tsv | tar -x -C "$stage"
job=$(sbatch --parsable --partition=nano --cpus-per-task=2 --mem=12G --time=00:30:00 \
 --job-name=Optix_nonPK --output="$run/logs/Optix_nonPK.%j.out" --error="$run/logs/Optix_nonPK.%j.err" \
 --wrap="micromamba run -n gynandrocalling python '$stage/scripts/figures/figure3/optix_arnaud_review_nonPK.py' --run '$run' --pk-catalog '$stage/scripts/figures/figure3/pk_dark_pale_samples.tsv' --output '$out'")
printf 'JOBID\t%s\nMAIN_PDF\t%s\nREVIEW_DIR\t%s\n' "$job" "$out/merged/Dennis_main_candidate_BA90.pdf" "$out"
