#!/bin/bash
#SBATCH --job-name=smudge6y120
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/smudge6y120_%j.out
#SBATCH --error=logs/smudge6y120_%j.err

set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

cd "$ROOT"

"$PYTHON" scripts/plot_smudge_sixpanel_viridis_y120.py

ls -lh \
  analysis/smudgeplot_4legs_combined/gynandromorph_controls_sixpanel_viridis_y120.*
