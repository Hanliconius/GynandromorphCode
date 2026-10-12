#!/usr/bin/env bash
set -euo pipefail
CODE=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo/scripts/yb
ROOT=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
mkdir -p "$ROOT/analysis/yb_local_both_references"
LOG="$ROOT/analysis/yb_local_both_references/tool_setup.log"
{
    if type module >/dev/null 2>&1; then
        module load BWA/0.7.17 || true
        module load samtools || true
        module load Java/1.8.0_60 || true
        module load deepTools || module load deeptools || true
        module load picard || module load Picard || true
    fi
} > "$LOG" 2>&1
printf 'YB_CHECK_START\n'
MODE=${1:---submit}
[[ "$MODE" == "--submit" || "$MODE" == "--check-tools" ]] || { echo "Use --submit or --check-tools" >&2; exit 2; }
exec micromamba run -n gynandrocalling python "$CODE/rebuild_Yb_local_both_references.py" "$MODE"
