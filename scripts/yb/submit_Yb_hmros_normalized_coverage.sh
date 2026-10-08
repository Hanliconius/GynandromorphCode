#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

echo "===== HmRos Yb NORMALIZED COVERAGE PREFLIGHT ====="

for f in     "$HERE/resolve_Yb_hmros_coordinates.py"     "$HERE/plot_Yb_hmros_panel_normalized_coverage.py"
do
    "$PYTHON" -m py_compile "$f"
done

bash -n "$HERE/run_Yb_hmros_normalized_coverage.cpu.slurm"

REF_PATH="$PROJECT_ROOT/analysis/genomewide_mixing_input/metadata/reference_fasta.path"
MANIFEST="$PROJECT_ROOT/analysis/genomewide_mixing_input/metadata/resolved_download_manifest.tsv"

for f in "$REF_PATH" "$MANIFEST"
do
    [[ -s "$f" ]] || {
        echo "ERROR: missing required input: $f" >&2
        exit 1
    }
done

REFERENCE=$(head -n 1 "$REF_PATH")

if [[ "$REFERENCE" != *"/GCF_056956725.1/"* ]]
then
    echo "ERROR: active reference is not the new HmRos GCF_056956725.1 assembly: $REFERENCE" >&2
    exit 1
fi

echo "REFERENCE=$REFERENCE"
echo "NORMALIZATION=library mapped-bases scaling + per-bin median of 18 controls"
echo "DEPTH_FILTERS=MAPQ>=20, baseQ>=20"
echo "T2/T3=normalized separately before side averaging"
echo "PREFLIGHT=PASS"
echo

job=$(
    sbatch         --parsable         --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,YB_CODE_DIR=$HERE"         "$HERE/run_Yb_hmros_normalized_coverage.cpu.slurm"
)
job=${job%%;*}

echo "JOBID=$job"
echo "OUT=$PROJECT_ROOT/logs/Yb_hmros_normalized.$job.out"
echo "ERR=$PROJECT_ROOT/logs/Yb_hmros_normalized.$job.err"
echo "COORDS=$PROJECT_ROOT/analysis/yb_hmros_normalized_coverage/metadata/Yb_hmros_coordinates.tsv"
echo "PDF=$PROJECT_ROOT/analysis/yb_hmros_normalized_coverage/results/Yb_hmros_panel_normalized_coverage.pdf"
echo "SUMMARY=$PROJECT_ROOT/analysis/yb_hmros_normalized_coverage/results/Yb_hmros_panel_normalized_coverage.summary.tsv"
echo
squeue -j "$job"
