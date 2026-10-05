#!/usr/bin/env bash
#SBATCH --job-name=Fig3optix
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Fig3_optix.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Fig3_optix.%j.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${FIGURE3_CODE_DIR:?FIGURE3_CODE_DIR must be exported by the submit helper}
PYTHON=${PYTHON:-$HOME/.local/share/mamba/envs/gynandrocalling/bin/python}

mkdir -p "$PROJECT_ROOT/logs"

echo "FIGURE3_OPTIX_START"
echo "job=${SLURM_JOB_ID}"
echo "host=$(hostname)"
echo "PROJECT_ROOT=$PROJECT_ROOT"
echo "CODE_DIR=$CODE_DIR"

"$PYTHON" "$CODE_DIR/Fig3_optix_collapsed_pattern_dosage.py"

echo "FIGURE3_OPTIX_COMPLETE"
