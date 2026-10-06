#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

INFORMATIVE="$PROJECT_ROOT/analysis/original_502/optix_deltaAF_haplotype_analysis/optix_B_D_R_deltaAF_informative_SNPs.tsv"
LR="$PROJECT_ROOT/analysis/original_502/reciprocal_L_R_mixture/optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"

echo "===== FIGURE 3 ORIGINAL INFORMATIVE-SITE PREFLIGHT ====="
for f in "$INFORMATIVE" "$LR"; do
    if [[ ! -s "$f" ]]; then
        echo "ERROR: missing/empty input: $f" >&2
        exit 1
    fi
done

"$PYTHON" -m py_compile "$HERE/Fig3_optix_first_informative.py"
bash -n "$HERE/run_Fig3_optix_first_informative.nano.slurm"
echo "PREFLIGHT=PASS"

job=$(
    sbatch --parsable         --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,FIGURE3_CODE_DIR=$HERE"         "$HERE/run_Fig3_optix_first_informative.nano.slurm"
)
job=${job%%;*}

echo "JOBID=$job"
echo "OUT=$PROJECT_ROOT/logs/Fig3_first_informative.$job.out"
echo "ERR=$PROJECT_ROOT/logs/Fig3_first_informative.$job.err"
echo "PDF=$PROJECT_ROOT/analysis/figure3/optix_first_informative_sites/optix_first_informative_B_D_R_collapsed_plus_recalled_LR.pdf"
squeue -j "$job"
