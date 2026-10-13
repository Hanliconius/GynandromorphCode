#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
OUT="$PROJECT/analysis/figure3/optix_simple_evidence/$(date -u +%Y%m%dT%H%M%SZ)_$$"
mkdir -p "$OUT/logs"
cp "$HERE/plot_optix_simple_evidence.py" "$HERE/Fig3_optix_full_range_frequency.py" "$HERE/Fig3_optix_BA099_both_callable.py" "$OUT/"
printf -v CMD '%q ' micromamba run -n gynandrocalling python "$OUT/plot_optix_simple_evidence.py" --output-dir "$OUT"
JOB=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=24G --time=02:00:00 --job-name=optixSimple --chdir="$OUT" --output="$OUT/logs/plot.%j.out" --error="$OUT/logs/plot.%j.err" --wrap="$CMD")
JOB=${JOB%%;*}
printf 'JOBID\t%s\nPDF\t%s\nERR\t%s\n' "$JOB" "$OUT/optix_simple_evidence.pdf" "$OUT/logs/plot.$JOB.err"
