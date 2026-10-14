#!/bin/bash
#SBATCH --job-name=melpSNPcount
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/melpomene_polymorphic_sites.%j.out
#SBATCH --error=logs/melpomene_polymorphic_sites.%j.err

set -euo pipefail

bcftools=/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/bcftools

bcf=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/genomewide_mixing_input/variant_calls/joint_24/genomewide/all_24_samples.joint_at_control_sites.chromosomes_1_to_21.bcf

samples="ERR3790122,ERR3790123,ERR2298257,ERR2298260,ERR3790118,ERR2298201,ERR2298216,ERR1307113,ERR1307112,ERR1307111,ERR3790124,ERR2298208,ERR2298221"

echo "PURE_MELPOMENE_SAMPLES=13"

count=$(
    "$bcftools" view \
        -s "$samples" \
        -c 1:minor \
        -H \
        "$bcf" \
    | wc -l
)

printf 'pure_melpomene_polymorphic_sites\t%s\n' "$count"
