#!/usr/bin/env bash

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PROJECT_ROOT="${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}"

mkdir -p \
    "$PROJECT_ROOT/logs" \
    "$PROJECT_ROOT/analysis/figure2/panels"

# Panels A and F have been dropped from Figure 2.
rm -f \
    "$PROJECT_ROOT/analysis/figure2/panels/Fig2A.pdf" \
    "$PROJECT_ROOT/analysis/figure2/panels/Fig2F.pdf"

submission="$(
    sbatch \
        --parsable \
        --array=0-5%6 \
        --output="$PROJECT_ROOT/logs/Fig2_panel.%A_%a.out" \
        --error="$PROJECT_ROOT/logs/Fig2_panel.%A_%a.err" \
        --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,FIGURE2_CODE_DIR=$HERE" \
        "$HERE/run_Figure2_panels.nano.slurm"
)"

jobid="${submission%%;*}"

echo "FIGURE2_ARRAY_JOB=$jobid"
echo "FIGURE2_CODE_DIR=$HERE"
echo "GYNANDRO_ROOT=$PROJECT_ROOT"
echo
echo "TASK MAP:"
echo "  0 = B chromosome coverage"
echo "  1 = C Smudgeplot controls"
echo "  2 = D Smudgeplot gynandromorph"
echo "  3 = E AF marginal"
echo "  4 = G simulated AF"
echo "  5 = H observed AF"
echo
echo "Check with:"
echo "  squeue -j $jobid"
echo "  sacct -j $jobid -X --format=JobID,JobName,Partition,State,Elapsed,ExitCode"
