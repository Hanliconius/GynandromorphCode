#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PYTHON=${PYTHON:-$HOME/.local/share/mamba/envs/gynandrocalling/bin/python}

if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: expected Python not executable: ${PYTHON}" >&2
    exit 1
fi

echo "===== PYTHON SYNTAX ====="
"${PYTHON}" -m py_compile     "${HERE}/common.py"     "${HERE}"/Fig2?_*.py

echo "===== SHELL SYNTAX ====="
bash -n "${HERE}/run_Figure2_panels.sh"
bash -n "${HERE}/run_Figure2_panels.nano.slurm"
bash -n "${HERE}/submit_Figure2_panels.sh"

echo "FIGURE2_CODE_PREFLIGHT=PASS"
