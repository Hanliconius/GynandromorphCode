#!/bin/bash
#SBATCH --job-name=smudge6fig
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/smudge6fig_%j.out
#SBATCH --error=logs/smudge6fig_%j.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

cd "$ROOT"

test -s scripts/plot_smudge_sixpanel_viridis.py

echo "HOST       = $(hostname)"
echo "PYTHON     = $PYTHON"
echo "START      = $(date)"
echo

"$PYTHON" \
    scripts/plot_smudge_sixpanel_viridis.py

echo
echo "===== OUTPUTS ====="

ls -lh \
    analysis/smudgeplot_4legs_combined/gynandromorph_controls_sixpanel_viridis.pdf \
    analysis/smudgeplot_4legs_combined/gynandromorph_controls_sixpanel_viridis.svg \
    analysis/smudgeplot_4legs_combined/gynandromorph_controls_sixpanel_viridis.png

echo
echo "FINISH     = $(date)"
