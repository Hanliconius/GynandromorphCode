#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"

echo "===== Yb EXACT COORDINATE TRANSFER PREFLIGHT ====="

"$PYTHON" -m py_compile \
  "$HERE/analyze_Yb_Hmel2_5_troughs.py" \
  "$HERE/map_Hmel2_5_to_HmRos.py"

bash -n "$HERE/run_Yb_exact_coordinate_transfer.cpu.slurm"

for f in \
  "$PROJECT_ROOT/analysis/yb_Hmel2_5_normalized_coverage/results/Yb_Hmel2.5_Figure8_normalized_50bp.long.tsv.gz" \
  "$PROJECT_ROOT/analysis/yb_Hmel2_5_normalized_coverage/results/Yb_Hmel2.5_Figure8_normalized_50bp.side_means.tsv.gz" \
  "$PROJECT_ROOT/analysis/genomewide_mixing_input/metadata/reference_fasta.path"
do
  [[ -s "$f" ]] || {
    echo "ERROR: missing input: $f" >&2
    exit 1
  }
done

echo "PREFLIGHT=PASS"
echo "OLD_REFERENCE=Hmel2.5 / H. m. melpomene"
echo "NEW_REFERENCE=GCF_056956725.1 HmRos / H. m. rosina"
echo "FIRST_STEP=locate exact 50-bp coverage troughs"
echo "SECOND_STEP=sequence-based coordinate transfer if Hmel2.5 FASTA is already present"
echo

job=$(
  sbatch \
    --parsable \
    --partition=cpu \
    --cpus-per-task=1 \
    --mem=8G \
    --time=00:30:00 \
    --export="ALL,GYNANDRO_ROOT=$PROJECT_ROOT,YB_CODE_DIR=$HERE" \
    "$HERE/run_Yb_exact_coordinate_transfer.cpu.slurm"
)
job=${job%%;*}

echo "JOBID=$job"
echo "OUT=$PROJECT_ROOT/logs/Yb_exact_liftover.$job.out"
echo "ERR=$PROJECT_ROOT/logs/Yb_exact_liftover.$job.err"
echo "RESULTS=$PROJECT_ROOT/analysis/yb_exact_coordinate_transfer/results"
echo
squeue -j "$job"
