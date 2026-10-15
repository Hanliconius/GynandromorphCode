#!/usr/bin/env bash
# Read-only validation of the existing audit plus an isolated Ivory checkout.
set -euo pipefail
audit="${1:-/gpfs/automountdir/gpfs/scratch/martinlab/joe/publication_audit_20261014_143331}"
r_env="${2:-syntangle_test}"
[[ -d "$audit" ]] || { echo "Audit directory missing: $audit" >&2; exit 1; }
out="$(mktemp -d "$audit/followup_XXXXXXXX")"
echo "FOLLOWUP_DIR=$out"
retry="$audit/recheck_20261014_192703"
{
  for project in SynTangle genome-observatory-live; do
    printf '%s ' "$project"
    git -C "$audit/$project/validation_checkout" rev-parse HEAD
  done
  if [[ -f "$retry/installed_versions.txt" ]]; then
    cat "$retry/installed_versions.txt"
  else
    echo "MISSING: installed_versions.txt"
  fi
} > "$out/python_evidence.txt"
echo "PYTHON_EVIDENCE=$out/python_evidence.txt"
grep -E '^(SynTangle |genome-observatory-live |matplotlib==|numpy==|scipy==|cvxpy==|scs==|highspy==|pycountry==)' "$out/python_evidence.txt" || true
git clone --quiet --single-branch --branch publication-readiness-20261014 \
  git@github.com:Hanliconius/ivory-mir193_paper.git "$out/ivory" \
  > "$out/clone.log" 2>&1
git -C "$out/ivory" rev-parse HEAD > "$out/ivory_commit.txt"
set +e
micromamba run -n "$r_env" Rscript "$out/ivory/tools/validate_maintained_code.R" "$out/ivory" \
  > "$out/ivory_validation.log" 2>&1
status=$?
set -e
echo "IVORY_EXIT_CODE=$status"
echo "IVORY_LOG=$out/ivory_validation.log"
if grep -q '^MAINTAINED_R_FILES' "$out/ivory_validation.log"; then
  grep -E '^(MAINTAINED_R_FILES|PARSE_FAIL|PACKAGE|FULL_DATA_REPRODUCTION)' "$out/ivory_validation.log"
else
  tail -n 8 "$out/ivory_validation.log"
fi
exit "$status"
