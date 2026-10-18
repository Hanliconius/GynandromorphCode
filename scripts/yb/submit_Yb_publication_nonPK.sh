#!/usr/bin/env bash
set -euo pipefail
repo=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo
root=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
stage="$root/analysis/yb_publication_code_df8f2b5"
out="$root/analysis/yb_publication_nonPK_v1"
if [ -e "$out" ]; then printf 'Output exists: %s\n' "$out" >&2; exit 1; fi
mkdir -p "$stage" "$stage/logs"
git -C "$repo" show df8f2b5b378ba40429a4b49bc24a62dc3a05651e:scripts/yb/plot_Yb_publication_nonPK.py > "$stage/plot_Yb_publication_nonPK.py"
job=$(sbatch --parsable --partition=nano --cpus-per-task=1 --mem=6G --time=00:15:00 --job-name=Yb_panel \
 --output="$stage/logs/Yb_panel.%j.out" --error="$stage/logs/Yb_panel.%j.err" \
 --wrap="micromamba run -n gynandrocalling python '$stage/plot_Yb_publication_nonPK.py' --output '$out'")
printf 'JOBID\t%s\nPDF\t%s\nLOG\t%s\n' "$job" "$out/Yb_d_v_nonPK_compact_Rosina.pdf" "$stage/logs/Yb_panel.${job%%;*}.out"
