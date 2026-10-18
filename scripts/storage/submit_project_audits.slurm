#!/bin/bash
#SBATCH --job-name=Gyn_audit
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/storage_audit_code/audit.%j.out
#SBATCH --error=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/storage_audit_code/audit.%j.err
set -euo pipefail
# SLURM copies the script, so receive the frozen source directory explicitly.
code_dir=${1:?Pass absolute frozen scripts directory}
root=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
repo=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandromorphCode_repo
stamp=$(date -u +%Y%m%dT%H%M%SZ)_${SLURM_JOB_ID}
seq="$root/analysis/sequencing_audits/$stamp"
out="$root/analysis/storage_audits/$stamp"
micromamba run -n gynandrocalling python "$code_dir/yb/audit_Yb_expansion_sequencing.py" \
 --run "$root/analysis/yb_ventral_population_expansion_v1" --output "$seq"
micromamba run -n gynandrocalling python "$code_dir/storage/audit_project_storage.py" \
 --root "$root" --output "$out" \
 --protect-run "$root/analysis/yb_ventral_population_expansion_v1" \
 --sequencing-report "$seq" --push-repo "$repo"
