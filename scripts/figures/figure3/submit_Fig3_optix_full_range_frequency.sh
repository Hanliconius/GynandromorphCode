#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"
VCF="$PROJECT/external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf"
[[ -s "$VCF" ]] || { echo "ERROR: local original optix VCF missing: $VCF" >&2; exit 2; }
[[ -x "$PYTHON" ]] || { echo "ERROR: Python missing: $PYTHON" >&2; exit 2; }
"$PYTHON" -m py_compile "$HERE/Fig3_optix_full_range_frequency.py" "$HERE/Fig3_optix_BA099_both_callable.py"
OUT="$PROJECT/analysis/figure3/optix_full_range_frequency/$(date -u +%Y%m%dT%H%M%SZ)_$$"
mkdir -p "$OUT/logs"
cp "$HERE/Fig3_optix_full_range_frequency.py" "$HERE/Fig3_optix_BA099_both_callable.py" "$OUT/"
printf -v COMMAND '%q ' "$PYTHON" "$OUT/Fig3_optix_full_range_frequency.py" --vcf "$VCF" --output-dir "$OUT"
JOB=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=16G --time=01:00:00 \
    --job-name=optixFreq --chdir="$OUT" --output="$OUT/logs/track.%j.out" \
    --error="$OUT/logs/track.%j.err" --wrap="$COMMAND")
JOB=${JOB%%;*}
printf 'JOBID\t%s\nPDF\t%s\nOUT\t%s\nERR\t%s\n' \
    "$JOB" "$OUT/optix_full_range_BDR_frequency_adjusted_LR.pdf" \
    "$OUT/logs/track.$JOB.out" "$OUT/logs/track.$JOB.err"
