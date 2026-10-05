#!/bin/bash
#SBATCH --job-name=founder_3v4
#SBATCH --partition=cpu
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --output=logs/159_three_vs_four_founder.%j.out
#SBATCH --error=logs/159_three_vs_four_founder.%j.err

set -u

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

working_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

date --iso-8601=seconds

"${working_python}" \
    scripts/159_fit_autosomal_three_vs_four_founder.py \
    2>&1 \
    | tee "logs/159_three_vs_four_founder.${SLURM_JOB_ID}.combined.log"

run_status=${PIPESTATUS[0]}

printf '\nrun_status=%s\n' "${run_status}"

output_status=0

for output in \
    analysis/autosomal_three_vs_four_founder/observed_model_fits.tsv \
    analysis/autosomal_three_vs_four_founder/observed_optimizer_multistart.tsv \
    analysis/autosomal_three_vs_four_founder/observed_site_model_evidence.tsv.gz \
    analysis/autosomal_three_vs_four_founder/observed_10kb_block_model_evidence.tsv \
    analysis/autosomal_three_vs_four_founder/observed_output_checksums.sha256
do
    if [[ -s "${output}" ]]
    then
        printf 'PRESENT\t%s\t%s bytes\n' \
            "${output}" \
            "$(wc -c < "${output}")"
    else
        printf 'MISSING_OR_EMPTY\t%s\n' "${output}"
        output_status=1
    fi
done

if [[
    "${run_status}" -eq 0 &&
    "${output_status}" -eq 0
]]
then
    printf '\nOBSERVED MODEL FITS\n'

    cat \
        analysis/autosomal_three_vs_four_founder/observed_model_fits.tsv

    printf '\nTEN-KB BLOCK EVIDENCE\n'

    cat \
        analysis/autosomal_three_vs_four_founder/observed_10kb_block_model_evidence.tsv
fi

printf '\nFINAL STATUS\n'
printf 'run_status=%s\n' "${run_status}"
printf 'output_status=%s\n' "${output_status}"

date --iso-8601=seconds

test \
    "${run_status}" -eq 0 \
    -a "${output_status}" -eq 0
