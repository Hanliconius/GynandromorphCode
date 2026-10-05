#!/usr/bin/env bash

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
gyn_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python
gyn_rscript=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/Rscript
analysis_script=scripts/142_B_western_postman_vs_western_Amazon_DR_adjusted_LR.py
plotter=scripts/142_plot_B_western_postman_vs_western_Amazon_DR_adjusted_LR.R
slurm_script=scripts/142_B_western_postman_vs_western_Amazon_DR.nano.slurm

cd "${project}" || {
    printf 'Cannot enter project: %s\n' "${project}"
}

mkdir -p logs

compile_status=127
preflight_status=127
plotter_status=127
static_status=127
submit_status=127

if [[ -s "${analysis_script}" ]]
then
    "${gyn_python}" -m py_compile "${analysis_script}" \
        > logs/142_clean_compile.log 2>&1
    compile_status=$?
else
    printf 'Missing analysis script: %s\n' "${analysis_script}" \
        > logs/142_clean_compile.log
fi

cat logs/142_clean_compile.log

if [[ "${compile_status}" -eq 0 ]]
then
    "${gyn_python}" "${analysis_script}" --prepare-only \
        > logs/142_clean_preflight.log 2>&1
    preflight_status=$?
else
    printf 'Preflight not run because compilation failed.\n' \
        > logs/142_clean_preflight.log
fi

cat logs/142_clean_preflight.log

if [[ "${preflight_status}" -eq 0 && -s "${plotter}" ]]
then
    "${gyn_rscript}" \
        -e "invisible(parse(file='${plotter}')); cat('R syntax check passed.\n')" \
        > logs/142_clean_plotter_parse.log 2>&1
    plotter_status=$?
else
    printf 'Plotter parse not run: preflight failed or plotter is missing.\n' \
        > logs/142_clean_plotter_parse.log
fi

cat logs/142_clean_plotter_parse.log

if [[ "${preflight_status}" -eq 0 ]]
then
    "${gyn_python}" - <<'PYEOF' \
        > logs/142_clean_static_independence_audit.log 2>&1
from pathlib import Path

path = Path("scripts/142_B_western_postman_vs_western_Amazon_DR_adjusted_LR.py")
text = path.read_text()

score_at = text.index("scanner.main()")
focal_read_at = text.index("mixture = pd.read_csv")

if focal_read_at <= score_at:
    raise RuntimeError("Focal mixture calls are read before comparative scoring")

required = (
    'if sample in original.FOCAL_SAMPLES:',
    'sample == LOCKED_EXCLUSION',
    '"reciprocal_mixture_posterior95_only"',
    'where(confidence, "./.")',
)

missing = [token for token in required if token not in text]
if missing:
    raise RuntimeError("Static safeguards missing: " + repr(missing))

print("comparative_scoring_before_focal_read=true")
print("focal_samples_explicitly_excluded=true")
print("locked_ERR3653308_exclusion=true")
print("reciprocal_mixture_only=true")
print("posterior_threshold=0.95")
print("raw_GT_fallback=false")
PYEOF
    static_status=$?
else
    printf 'Static audit not run because preflight failed.\n' \
        > logs/142_clean_static_independence_audit.log
fi

cat logs/142_clean_static_independence_audit.log

printf '\nFINAL VALIDATION STATUS\n'
printf 'python_compile\t%s\n' "${compile_status}"
printf 'analysis_preflight\t%s\n' "${preflight_status}"
printf 'plotter_parse\t%s\n' "${plotter_status}"
printf 'independence_audit\t%s\n' "${static_status}"

if [[ "${compile_status}" -eq 0 \
   && "${preflight_status}" -eq 0 \
   && "${plotter_status}" -eq 0 \
   && "${static_status}" -eq 0 \
   && -s "${slurm_script}" ]]
then
    job_id=$(sbatch --parsable "${slurm_script}")
    submit_status=$?
    if [[ "${submit_status}" -eq 0 ]]
    then
        printf '\nSubmitted job: %s\n' "${job_id}"
        printf 'Output: logs/142_B_western_postman_vs_western_Amazon_DR.nano.%s.out\n' "${job_id}"
        printf 'Error:  logs/142_B_western_postman_vs_western_Amazon_DR.nano.%s.err\n' "${job_id}"
    else
        printf '\nsbatch failed with status %s; no job was submitted.\n' "${submit_status}"
    fi
else
    printf '\nNo job submitted because one or more validations failed.\n'
fi

printf 'submission\t%s\n' "${submit_status}"
