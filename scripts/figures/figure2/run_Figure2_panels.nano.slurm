#!/usr/bin/env bash
#SBATCH --job-name=Fig2panel
#SBATCH --partition=nano
#SBATCH --array=0-7%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00

set -euo pipefail

PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}

# Slurm copies the submitted batch script into /var/spool/slurm before
# execution, so BASH_SOURCE[0] points at the spool copy on compute nodes.
# The submission helper therefore exports the real repository code path.
CODE_DIR=${FIGURE2_CODE_DIR:-${SLURM_SUBMIT_DIR}/scripts/figures/figure2}

PYTHON=${PYTHON:-$HOME/.local/share/mamba/envs/gynandrocalling/bin/python}

if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: expected Python not executable: ${PYTHON}" >&2
    exit 1
fi

if [[ ! -d "${CODE_DIR}" ]]; then
    echo "ERROR: Figure 2 code directory not found: ${CODE_DIR}" >&2
    exit 1
fi

panels=(
    Fig2A_side_specific_sequencing.py
    Fig2B_chromosome_coverage.py
    Fig2C_smudge_controls.py
    Fig2D_smudge_gynandromorph.py
    Fig2E_AF_marginal.py
    Fig2F_founder_scenarios.py
    Fig2G_simulated_AF.py
    Fig2H_observed_AF.py
)

script=${panels[${SLURM_ARRAY_TASK_ID}]}
script_path="${CODE_DIR}/${script}"

if [[ ! -s "${script_path}" ]]; then
    echo "ERROR: Figure 2 panel script not found: ${script_path}" >&2
    exit 1
fi

echo "FIGURE2_PANEL_START panel=${SLURM_ARRAY_TASK_ID} script=${script} host=$(hostname)"
echo "FIGURE2_CODE_DIR=${CODE_DIR}"
echo "GYNANDRO_ROOT=${PROJECT_ROOT}"

"${PYTHON}" "${script_path}" --project-root "${PROJECT_ROOT}"

echo "FIGURE2_PANEL_COMPLETE panel=${SLURM_ARRAY_TASK_ID} script=${script}"
