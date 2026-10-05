#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}

mkdir -p "$PROJECT_ROOT/logs" "$PROJECT_ROOT/analysis/figure2/panels"

PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

echo "===== FIGURE 2 BOSS-UPDATE PREFLIGHT ====="

bash -n "$HERE/run_Figure2_boss_updates.nano.slurm"

"$PYTHON" -m py_compile     "$HERE/common.py"     "$HERE/Fig2G_simulated_AF.py"     "$HERE/Fig2H_observed_AF.py"

for path in     "$PROJECT_ROOT/analysis/smudgeplot_standard_rerun_20261002/ERR1051012_diploid/L6/kmerpairs.smu"     "$PROJECT_ROOT/analysis/smudgeplot_lowerL_20261002/pseudoTetraploid_ERR1307111_ERR1307112/L6/kmerpairs.smu"     "$PROJECT_ROOT/analysis/smudgeplot_arrangements/TN_Left/kmerpairs.smu"     "$PROJECT_ROOT/analysis/smudgeplot_arrangements/TN_Right/kmerpairs.smu"     "$PROJECT_ROOT/analysis/smudgeplot_lowerL_20261002/TN_T2_merge/L6/kmerpairs.smu"     "$PROJECT_ROOT/analysis/smudgeplot_lowerL_20261002/TN_T3_merge/L6/kmerpairs.smu"
do
    [[ -s "$path" ]] || {
        echo "ERROR: missing/empty required input: $path" >&2
        exit 1
    }
done

echo "PREFLIGHT=PASS"

job=$(
    sbatch         --parsable         --output="$PROJECT_ROOT/logs/Fig2_boss.%A_%a.out"         --error="$PROJECT_ROOT/logs/Fig2_boss.%A_%a.err"         --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,FIGURE2_CODE_DIR=$HERE"         "$HERE/run_Figure2_boss_updates.nano.slurm"
)

job=${job%%;*}

echo "FIGURE2_BOSS_JOB=$job"
echo
echo "TASKS"
echo "0  Fig2C1 diploid control"
echo "1  Fig2C2 pseudo-tetraploid"
echo "2  Fig2D1 Left"
echo "3  Fig2D2 Right"
echo "4  Fig2D3 T2 merge"
echo "5  Fig2D4 T3 merge"
echo "6  Fig2G simulated AF"
echo "7  Fig2H observed AF"
echo
echo "LOGS=$PROJECT_ROOT/logs/Fig2_boss.${job}_*.{out,err}"
echo "OUTPUTS=$PROJECT_ROOT/analysis/figure2/panels"
echo
squeue -j "$job"
