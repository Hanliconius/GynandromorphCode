#!/usr/bin/env bash
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RUN=${1:?Provide the absolute existing Yb run directory}
[[ "$RUN" == /* && -s "$RUN/run.json" ]] || { echo "ERROR: invalid existing run: $RUN" >&2; exit 2; }
PYTHON="$HOME/.local/share/mamba/envs/gynandrocalling/bin/python"
missing=0
for ref in Hmel2.5 HmRos; do
    for sample in TN2L_S1 TN2R_S2 TN3L_S3 TN3R_S4; do
        for suffix in .dedup.50bp.bw .complete.json; do
            file="$RUN/$ref/coverage/$sample$suffix"
            if [[ ! -s "$file" ]]; then
                printf 'MISSING\t%s\n' "$file" >&2
                missing=1
            fi
        done
    done
done
[[ "$missing" == 0 ]] || { echo "No plot jobs submitted; check remaining alignment status." >&2; exit 2; }
"$PYTHON" -W error::SyntaxWarning -m py_compile "$HERE/plot_Yb_local_both_references.py"
RETRY="$RUN/plot_retry_$(date -u +%Y%m%dT%H%M%SZ)_$$"
mkdir -p "$RETRY"
cp "$HERE/plot_Yb_local_both_references.py" "$RETRY/"
{
    printf '#!/usr/bin/env bash\nset -euo pipefail\n'
    printf 'RUN=%q\nPYTHON=%q\nPLOT=%q\n' "$RUN" "$PYTHON" "$RETRY/plot_Yb_local_both_references.py"
    cat <<'WORKER'
refs=(Hmel2.5 HmRos)
ref=${refs[$SLURM_ARRAY_TASK_ID]}
exec "$PYTHON" "$PLOT" --reference "$ref" --input-dir "$RUN/$ref/coverage" --output-dir "$RUN/$ref/plots"
WORKER
} > "$RETRY/worker.sh"
job=$(sbatch --parsable --partition=cpu --cpus-per-task=1 --mem=8G --time=00:30:00 \
    --array=0-1%2 --job-name=YbReplot --output="$RETRY/plot.%A_%a.out" \
    --error="$RETRY/plot.%A_%a.err" "$RETRY/worker.sh")
printf 'PLOT_ARRAY\t%s\nLOG_DIR\t%s\n' "${job%%;*}" "$RETRY"
