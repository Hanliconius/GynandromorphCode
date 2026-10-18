#!/usr/bin/env bash
set -euo pipefail
repo=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo
root=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
stage="$root/analysis/yb_literature_code_b01f846"
out="$root/analysis/yb_literature_checks_v2"
if [ -e "$out" ]; then printf 'Output exists: %s\n' "$out" >&2; exit 1; fi
mkdir -p "$stage/logs"
git -C "$repo" show b01f846077627af1b26c6ccdee96962f6ec399d3:scripts/yb/plot_Yb_literature_checks.py > "$stage/plot_Yb_literature_checks.py"
job=$(sbatch --parsable --partition=nano --cpus-per-task=1 --mem=6G --time=00:25:00 --job-name=Yb_check \
 --output="$stage/logs/Yb_check.%j.out" --error="$stage/logs/Yb_check.%j.err" \
 --wrap="micromamba run -n gynandrocalling python '$stage/plot_Yb_literature_checks.py' --output '$out'")
printf 'JOBID\t%s\nPDF\t%s\nCYDNO_INVENTORY\t%s\nLOG\t%s\n' "$job" "$out/Yb_d_v_nonPK_compact_Rosina.pdf" "$out/cydno_pachinus_inventory.tsv" "$stage/logs/Yb_check.${job%%;*}.out"
