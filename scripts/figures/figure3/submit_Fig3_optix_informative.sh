#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

INFORMATIVE="$PROJECT_ROOT/analysis/original_502/optix_deltaAF_haplotype_analysis/optix_B_D_R_deltaAF_informative_SNPs.tsv"

mkdir -p "$PROJECT_ROOT/logs"

echo "===== FIGURE 3 INFORMATIVE-ONLY PREFLIGHT ====="

if [[ ! -s "$INFORMATIVE" ]]; then
    echo "ERROR: original stage-117a informative-site table is missing:" >&2
    echo "  $INFORMATIVE" >&2
    exit 1
fi

echo "INFORMATIVE_SITE_TABLE=$INFORMATIVE"
echo "Selection rule: |delta AF| >= 0.80; both comparative groups >= 0.80 callable; DP >= 2; GQ >= 20"
echo "TN_Left/TN_Right were excluded from site selection and polarization."

"$PYTHON" -m py_compile "$HERE/Fig3_optix_informative_only.py"
bash -n "$HERE/run_Fig3_optix_informative.nano.slurm"

echo "PREFLIGHT=PASS"

job=$(
    sbatch         --parsable         --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,FIGURE3_CODE_DIR=$HERE"         "$HERE/run_Fig3_optix_informative.nano.slurm"
)
job=${job%%;*}

echo
echo "FIGURE3_OPTIX_INFORMATIVE_JOB=$job"
echo "LOG_OUT=$PROJECT_ROOT/logs/Fig3_optix_informative.$job.out"
echo "LOG_ERR=$PROJECT_ROOT/logs/Fig3_optix_informative.$job.err"
echo "OUTPUT_DIR=$PROJECT_ROOT/analysis/figure3/optix_informative_pattern_dosage"
echo
squeue -j "$job"
