#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

BAND_BA="$PROJECT_ROOT/analysis/B_western_postman_vs_western_Amazon_DR_reanalysis/checkpoints/B_western_BA_and_adjusted_LR.tsv"
DR_BA="$PROJECT_ROOT/analysis/original_502/Wallbank_Dennis_Ray_BA/Wallbank_Dennis_Ray_dominant_BA.tsv"
LR="$PROJECT_ROOT/analysis/original_502/reciprocal_L_R_mixture/optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"

echo "===== FIGURE 3 BA>0.99 + BOTH RECALLED CALLABLE PREFLIGHT ====="

for f in "$BAND_BA" "$DR_BA" "$LR"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing/empty input: $f" >&2
        exit 1
    fi
    ls -lh "$f"
done

"$PYTHON" -m py_compile "$HERE/Fig3_optix_BA099_union_both_callable.py"
bash -n "$HERE/run_Fig3_optix_BA099_both_callable.nano.slurm"

echo "PREFLIGHT=PASS"

job=$(
    sbatch --parsable         --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,FIGURE3_CODE_DIR=$HERE"         "$HERE/run_Fig3_optix_BA099_both_callable.nano.slurm"
)
job=${job%%;*}

echo
echo "JOBID=$job"
echo "OUT=$PROJECT_ROOT/logs/Fig3_BA099_callable.$job.out"
echo "ERR=$PROJECT_ROOT/logs/Fig3_BA099_callable.$job.err"
echo "PDF=$PROJECT_ROOT/analysis/figure3/optix_BA099_callable_pattern_dosage/optix_BAgt099_union_both_recalled_callable_five_loci.pdf"
echo "SITES=$PROJECT_ROOT/analysis/figure3/optix_BA099_callable_pattern_dosage/BAgt099_union_both_recalled_callable_sites.tsv"
echo
squeue -j "$job"
