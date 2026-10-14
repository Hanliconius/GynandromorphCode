#!/usr/bin/env bash

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
cd "${project}"

root="${project}/analysis/pixy_pi_17controls_allsites"
chunks="${root}/metadata/chunks.tsv"

bash scripts/230_prepare_pixy_pi_17controls.sh

chunk_max=$(
    awk 'END {print NR - 2}' \
        "${chunks}"
)

chunk_n=$((chunk_max + 1))

call_submission=$(
    sbatch \
        --parsable \
        --array="0-${chunk_max}%48" \
        scripts/231_call_pixy_allsites_chunks.nano.slurm
)
call_job=${call_submission%%;*}

pixy_submission=$(
    sbatch \
        --parsable \
        --dependency="afterok:${call_job}" \
        --array="0-${chunk_max}%80" \
        scripts/232_run_pixy_pi_chunks.nano.slurm
)
pixy_job=${pixy_submission%%;*}

collect_submission=$(
    sbatch \
        --parsable \
        --dependency="afterok:${pixy_job}" \
        scripts/233_collect_pixy_pi.nano.slurm
)
collect_job=${collect_submission%%;*}

echo
echo "===== SUBMITTED PIXY PI PIPELINE ====="
echo "chunk_tasks=${chunk_n}"
echo "all_sites_array=${call_job}       nano, max 48 simultaneous"
echo "pixy_array=${pixy_job}            nano, max 80 simultaneous; afterok:${call_job}"
echo "collect=${collect_job}             nano; afterok:${pixy_job}"
echo
echo "squeue -j ${call_job},${pixy_job},${collect_job}"
echo
echo "After completion:"
echo "cat ${root}/summary/pi_genomewide_by_population.tsv"
