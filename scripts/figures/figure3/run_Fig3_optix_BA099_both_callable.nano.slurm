#!/usr/bin/env bash
#SBATCH --job-name=Fig3BA99
#SBATCH --partition=nano
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=00:20:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Fig3_BA099_callable.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/logs/Fig3_BA099_callable.%j.err

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
CODE_DIR=${FIGURE3_CODE_DIR:?FIGURE3_CODE_DIR must be exported}
PYTHON=${PYTHON:-$HOME/.local/share/mamba/envs/gynandrocalling/bin/python}

mkdir -p "$PROJECT_ROOT/logs"

echo "FIG3_BA099_CALLABLE_START"
echo "job=${SLURM_JOB_ID}"
echo "host=$(hostname)"
echo "PROJECT_ROOT=$PROJECT_ROOT"

"$PYTHON" "$CODE_DIR/Fig3_optix_BA099_both_callable.py"

echo "FIG3_BA099_CALLABLE_COMPLETE"
