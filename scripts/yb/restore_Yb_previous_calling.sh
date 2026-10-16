#!/usr/bin/env bash
# Restore the pre-uniform Yb workflow; no calling/selection changes.
set -euo pipefail
repo=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo
root=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
pin=738cdaa95c04debceba8a6d949efb856b2cd0b34
counts="$root/analysis/yb_control_overlays/20261013T174941Z_c116f9a5/interval_SNP_genoplot/20261013T232705Z"
fit="$root/analysis/original_502/reciprocal_L_R_mixture/optix_reciprocal_L_R_mixture.summary.tsv"
focal="$root/analysis/yb_local_both_references/20261012T194933Z_7fa9de96"
stamp=$(date -u +%Y%m%dT%H%M%SZ)
restore="$root/analysis/yb_previous_calling_restored/$stamp"
mkdir -p "$restore/code" "$restore/logs"
git -C "$repo" cat-file -e "$pin^{commit}"
git -C "$repo" archive "$pin" | tar -x -C "$restore/code"
micromamba run -n gynandrocalling python - "$restore" "$counts" "$fit" "$focal" "$pin" <<'PY'
import hashlib,json,sys
from pathlib import Path
import numpy,pandas,pysam,scipy,matplotlib
out,counts,fit,focal=map(Path,sys.argv[1:5])
files=[counts/f'{ref}_Yb{i}_SNP_counts.tsv.gz'
       for ref in ('Hmel2.5','HmRos') for i in ('d','v')]
files += [fit,out/'code/scripts/yb/yb_corrected_genotypes.py',
          out/'code/scripts/reciprocal_mixing/134_fit_reciprocal_L_R_mixture_genotypes.py']
for p in files:
    if not p.is_file():raise SystemExit('Missing original input: '+str(p))
manifest=dict(code_commit=sys.argv[5],counts_run=str(counts),fit_summary=str(fit),
    focal_run=str(focal),mode='Exact previous caller and stored counts; original five PKs and extra race strips',
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
    new_uniform_calls_used=False,new_Optix_fit_used=False,
    SNP_selection_changed=False,original_outputs_preserved=True)
(out/'restoration.json').write_text(json.dumps(manifest,indent=2)+'\n')
PY
cat > "$restore/worker.sh" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec micromamba run -n gynandrocalling python \
  "$restore/code/scripts/yb/yb_corrected_genotypes.py" \
  --counts-run "$counts" --fit-summary "$fit" --focal-run "$focal" \
  --add-ready-pks --add-race-strips
EOF
job=$(sbatch --parsable --partition=nano --time=00:30:00 --mem=8G \
  --job-name=YbRestore --output="$restore/logs/restore.%j.out" \
  --error="$restore/logs/restore.%j.err" "$restore/worker.sh")
printf 'RESTORE_JOB\t%s\nLOG\t%s/logs/restore.%s.out\nAUDIT\t%s/restoration.json\n' \
  "$job" "$restore" "${job%%;*}" "$restore"
