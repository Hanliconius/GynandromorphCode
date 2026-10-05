#!/usr/bin/env bash
#SBATCH --job-name=gyn183qc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=06:00:00
#SBATCH --output=logs/183_validate_corrected_inputs.%j.out
#SBATCH --error=logs/183_validate_corrected_inputs.%j.err

set -euo pipefail

start=$SECONDS
trap 'rc=$?; printf "STAGE_END status=%s elapsed_seconds=%s\n" "$rc" "$((SECONDS-start))"' EXIT

cd /gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling

root=analysis/genomewide_mixing_corrected_v1
mkdir -p "$root/audit"

report="$root/audit/input_validation.${SLURM_JOB_ID}.tsv"

printf 'chromosome\tstatus\trows\tbytes\tfile\n' > "$report"

passed=0
failed=0

for chromosome in {1..21}
do
    accession_number=$((140388 + chromosome))
    contig="NC_${accession_number}.1"
    input="$root/inputs/chr${chromosome}_${contig}.qc.tsv.gz"

    status=PASS
    rows=NA
    bytes=0

    if [[ ! -s "$input" ]]
    then
        status=MISSING_OR_EMPTY
    elif [[ ! -s "$input.complete.tsv" ]]
    then
        status=MISSING_COMPLETION_AUDIT
    elif [[ ! -s "$input.samples.txt" ]]
    then
        status=MISSING_SAMPLE_LIST
    elif [[ ! -s "$input.source.header.vcf" ]]
    then
        status=MISSING_SOURCE_HEADER
    elif ! gzip -t "$input"
    then
        status=GZIP_FAILED
    else
        rows=$(awk -F '\t' '$1=="rows" {print $2}' \
            "$input.complete.tsv")

        sample_count=$(wc -l < "$input.samples.txt")
        unique_count=$(sort -u "$input.samples.txt" | wc -l)

        if [[ ! "$rows" =~ ^[1-9][0-9]*$ ]]
        then
            status=INVALID_ROW_AUDIT
        elif [[ "$sample_count" -ne 24 || "$unique_count" -ne 24 ]]
        then
            status=INVALID_SAMPLE_LIST
        else
            for focal in \
                TN_Left_T2 TN_Left_T3 \
                TN_Right_T2 TN_Right_T3 \
                TN_Left_merged TN_Right_merged
            do
                if ! grep -Fxq "$focal" "$input.samples.txt"
                then
                    status=MISSING_FOCAL_SAMPLE
                fi
            done
        fi
    fi

    if [[ -s "$input" ]]
    then
        bytes=$(wc -c < "$input")
    fi

    if [[ "$status" == PASS ]]
    then
        passed=$((passed + 1))
    else
        failed=$((failed + 1))
    fi

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$chromosome" "$status" "$rows" "$bytes" "$input" \
        | tee -a "$report"
done

printf '\nVALIDATION SUMMARY\n'
printf 'passed_chromosomes=%s\n' "$passed"
printf 'failed_chromosomes=%s\n' "$failed"
printf 'report=%s\n' "$report"
printf 'scope=extraction integrity; not scientific model validation\n'

if [[ "$passed" -eq 21 && "$failed" -eq 0 ]]
then
    printf 'ENDPOINT=CORRECTED_INPUTS_READY\n'
    printf 'NEXT=implement and test restored model; do not reuse old fitting scripts\n'
else
    printf 'ENDPOINT=INPUT_REPAIR_REQUIRED\n'
fi

[[ "$passed" -eq 21 && "$failed" -eq 0 ]]
