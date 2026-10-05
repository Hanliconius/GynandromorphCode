#!/usr/bin/env bash
#SBATCH --job-name=rayFixed
#SBATCH --partition=nano
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/113_replicate_Ray_with_current_contigua.nano.%j.out
#SBATCH --error=logs/113_replicate_Ray_with_current_contigua.nano.%j.err

set -u

project="/scratch/martinlab/joe/GynandroCalling"
python_bin="/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/gynandrocalling/bin/python"

cd "${project}"

printf 'Host: %s\n' "$(hostname)"
printf 'Started: %s\n' "$(date --iso-8601=seconds)"

"${python_bin}" \
    scripts/113_replicate_Ray_with_current_contigua.py

printf 'Finished: %s\n' "$(date --iso-8601=seconds)"
