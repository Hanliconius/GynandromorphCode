#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON=${PYTHON:-$HOME/.local/share/mamba/envs/gynandrocalling/bin/python}

if [[ ! -x "${PYTHON}" ]]; then
    PYTHON=python
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

for script in "${panels[@]}"
do
    echo "===== ${script} ====="
    "${PYTHON}" "${HERE}/${script}" --project-root "${PROJECT_ROOT}"
done

echo
echo "===== FIGURE 2 PANELS COMPLETE ====="
ls -lh "${PROJECT_ROOT}/analysis/figure2/panels"/Fig2?.pdf
