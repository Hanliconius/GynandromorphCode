#!/usr/bin/env bash
set -euo pipefail
repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
root=${GYNANDRO_ROOT:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling}
mkdir -p "$root/logs" "$root/analysis/arnaud_mixing_supplement"
jobfile=$(mktemp)
trap 'rm -f "$jobfile"' EXIT
cat > "$jobfile" <<'JOB'
#!/usr/bin/env bash
set -euo pipefail
micromamba run -n gynandrocalling python "$1/scripts/figures/supplementary/build_arnaud_mixing_figures.py" --project-root "$2" "${@:3}"
JOB
sbatch --parsable --partition=nano --time=00:30:00 --cpus-per-task=2 --mem=12G --job-name=arnaud_AF --output="$root/logs/arnaud_AF_%j.log" "$jobfile" "$repo" "$root" "$@"
