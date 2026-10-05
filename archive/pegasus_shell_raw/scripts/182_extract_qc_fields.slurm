#!/usr/bin/env bash
#SBATCH --job-name=gyn182ex
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --output=logs/182_extract_qc_fields.%A_%a.out
#SBATCH --error=logs/182_extract_qc_fields.%A_%a.err

set -euo pipefail
start=$SECONDS
trap 'rc=$?; printf "STAGE_END status=%s elapsed_seconds=%s\n" "$rc" "$((SECONDS-start))"' EXIT

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
bcftools=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/bcftools
cd "$project"

printf 'STAGE_START task=%s host=%s time=%s\n' \
    "$SLURM_ARRAY_TASK_ID" "$(hostname)" "$(date -Is)"

manifest=analysis/genomewide_mixing_input/metadata/chromosome_calling_manifest.tsv

row=$(awk -F '\t' -v task="$SLURM_ARRAY_TASK_ID" '
NR==1 {
    for(i=1;i<=NF;i++) h[$i]=i
    if(!h["array_id"] || !h["assigned_molecule"] || !h["contig"]) exit 2
    next
}
$(h["array_id"])==task {
    print $(h["assigned_molecule"]),$(h["contig"])
    n++
}
END {if(n!=1) exit 3}
' "$manifest")

read -r molecule contig <<< "$row"

[[ "$molecule" =~ ^([1-9]|1[0-9]|2[01])$ ]]
[[ "$contig" =~ ^NC_[0-9]+\.[0-9]+$ ]]

source_dir="analysis/genomewide_mixing_input/variant_calls/joint_24/by_chromosome/chr${molecule}_${contig}"
bcf="$source_dir/all_24_samples.joint_at_control_sites.bcf"

[[ -s "$bcf" && -s "$bcf.csi" && -s "$source_dir/calling.complete.tsv" ]]

outdir=analysis/genomewide_mixing_corrected_v1/inputs
mkdir -p "$outdir"
output="$outdir/chr${molecule}_${contig}.qc.tsv.gz"

if [[ -e "$output" ]]; then
    printf 'STOP: output already exists; inspect before rerunning: %s\n' \
        "$output" >&2
    exit 1
fi

work=$(mktemp -d "$outdir/.extract_${molecule}.XXXXXX")

"$bcftools" query -l "$bcf" > "$work/samples.txt"
[[ $(wc -l < "$work/samples.txt") -eq 24 ]]

"$bcftools" view -h "$bcf" > "$work/source.header.vcf"

for field in GT AD DP PL SP; do
    grep -q "^##FORMAT=<ID=${field}," "$work/source.header.vcf"
done

expected=$("$bcftools" index -n "$bcf")

{
    printf '#CHROM\tPOS\tREF\tALT\tQUAL\tFILTER'

    while IFS= read -r sample; do
        printf '\t%s' "$sample"
    done < "$work/samples.txt"

    printf '\n'

    "$bcftools" query \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT:%AD:%DP:%PL:%SP]\n' \
        "$bcf"
} | gzip -1c > "$work/counts.tsv.gz"

gzip -t "$work/counts.tsv.gz"

actual=$(gzip -dc "$work/counts.tsv.gz" | awk -F '\t' '
NF!=30 {bad=1}
END {
    if(bad) exit 1
    print NR-1
}')

[[ "$actual" -eq "$expected" && "$actual" -gt 0 ]]

printf 'field\tvalue\nsource\t%s\nrows\t%s\nsample_format\tGT:AD:DP:PL:SP\nfilter_applied\tFalse\n' \
    "$bcf" "$actual" > "$work/audit.tsv"

mv "$work/source.header.vcf" "$output.source.header.vcf"
mv "$work/samples.txt" "$output.samples.txt"
mv "$work/counts.tsv.gz" "$output"
mv "$work/audit.tsv" "$output.complete.tsv"
rmdir "$work"

printf 'COMPLETE chromosome=%s rows=%s bytes=%s output=%s\n' \
    "$molecule" "$actual" "$(wc -c < "$output")" "$output"
