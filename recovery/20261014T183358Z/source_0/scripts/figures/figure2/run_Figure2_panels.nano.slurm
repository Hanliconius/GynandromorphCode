#!/usr/bin/env bash
#SBATCH --job-name=Fig2panel
#SBATCH --partition=nano
#SBATCH --array=0-5%6
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00

set -euo pipefail

PROJECT_ROOT="${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}"
CODE_DIR="${FIGURE2_CODE_DIR:?FIGURE2_CODE_DIR was not exported by the submission helper}"
PYTHON="${PYTHON:-$HOME/.local/share/mamba/envs/gynandrocalling/bin/python}"

if [[ ! -x "$PYTHON" ]]; then
    echo "ERROR: Python executable not found: $PYTHON" >&2
    exit 1
fi

if [[ ! -d "$CODE_DIR" ]]; then
    echo "ERROR: Figure 2 code directory not found: $CODE_DIR" >&2
    exit 1
fi

panels=(
    "Fig2B_chromosome_coverage.py"
    "Fig2C_smudge_controls.py"
    "Fig2D_smudge_gynandromorph.py"
    "Fig2E_AF_marginal.py"
    "Fig2G_simulated_AF.py"
    "Fig2H_observed_AF.py"
)

task="${SLURM_ARRAY_TASK_ID}"

if (( task < 0 || task >= ${#panels[@]} )); then
    echo "ERROR: invalid SLURM_ARRAY_TASK_ID=$task" >&2
    exit 1
fi

script="${panels[$task]}"
script_path="$CODE_DIR/$script"

if [[ ! -s "$script_path" ]]; then
    echo "ERROR: panel script not found: $script_path" >&2
    exit 1
fi

echo "FIGURE2_PANEL_START"
echo "task=$task"
echo "script=$script"
echo "host=$(hostname)"
echo "FIGURE2_CODE_DIR=$CODE_DIR"
echo "GYNANDRO_ROOT=$PROJECT_ROOT"

"$PYTHON" "$script_path" \
    --project-root "$PROJECT_ROOT"

echo "FIGURE2_PANEL_COMPLETE"
echo "task=$task"
echo "script=$script"
