#!/bin/bash
#SBATCH --job-name=NDJalign
#SBATCH --partition=cpu
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=logs/217_NDJalign.%A_%a.out
#SBATCH --error=logs/217_NDJalign.%A_%a.err
set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
alternate=/gpfs/automountdir/gpfs/scratch/martinlab/joe/Gynandrmorph
tools=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
python_bin=${project}/pairwise-af/.venv-af/bin/python
cd "${project}"
reference=$(< analysis/genomewide_mixing_input/metadata/reference_fasta.path)
output_root=${project}/analysis/NDJ_control_comparison_v1/alignment

for executable in "${tools}/bwa" "${tools}/samtools" "${python_bin}"
do
    test -x "${executable}"
done
"${python_bin}" -c 'import pysam'
for input in "${reference}" "${reference}.fai" "${reference}.amb" "${reference}.ann" "${reference}.bwt" "${reference}.pac" "${reference}.sa" "${alternate}/GynL.R1.fq.gz" "${alternate}/GynL.R2.fq.gz" "${alternate}/GynR.R1.fq.gz" "${alternate}/GynR.R2.fq.gz"
do
    if [[ ! -s "${input}" ]]; then
        printf 'MISSING_OR_EMPTY\t%s\n' "${input}" >&2
        false
    fi
done

if [[ "${1:-}" == --preflight ]]; then
    printf 'PREFLIGHT_PASS\treference=%s\n' "${reference}"
else
    case "${SLURM_ARRAY_TASK_ID:?Run through sbatch}" in
        0) old=GynL; sample=NDJ_left ;;
        1) old=GynR; sample=NDJ_right ;;
        *) printf 'Invalid array task\n' >&2; false ;;
    esac
    mkdir -p "${output_root}"
    final=${output_root}/${sample}
    if [[ -e "${final}" ]]; then
        printf 'Output already exists; refusing overwrite: %s\n' "${final}" >&2
        false
    fi
    lock=${output_root}/${sample}.running
    mkdir "${lock}"
    work=""
    cleanup() {
        status=$?
        trap - EXIT
        if [[ -n "${work}" && -d "${work}" ]]; then
            rm -rf -- "${work}"
        fi
        rmdir "${lock}"
        printf 'STAGE_END\tsample=%s\tstatus=%s\n' "${sample}" "${status}"
        exit "${status}"
    }
    trap cleanup EXIT
    work=$(mktemp -d "${output_root}/.${sample}.${SLURM_JOB_ID}.XXXXXX")
    printf 'STAGE_START\tsample=%s\treference=%s\n' "${sample}" "${reference}"

    # BWA requires literal backslash-t sequences here, NOT embedded tabs.
    rg="@RG\\tID:${sample}\\tSM:${sample}\\tLB:${sample}\\tPL:ILLUMINA\\tPU:${sample}"
    "${tools}/bwa" mem -M -t 12 -R "${rg}" "${reference}" \
        "${alternate}/${old}.R1.fq.gz" "${alternate}/${old}.R2.fq.gz" \
        | "${tools}/samtools" sort -n -@ 2 -m 2G \
            -T "${work}/namesort" -o "${work}/names.bam" -
    "${tools}/samtools" fixmate -m -@ 3 "${work}/names.bam" "${work}/fixmate.bam"
    rm -f -- "${work}/names.bam"
    "${tools}/samtools" sort -@ 7 -m 3G -T "${work}/coordsort" \
        -o "${work}/coordinates.bam" "${work}/fixmate.bam"
    rm -f -- "${work}/fixmate.bam"
    "${tools}/samtools" markdup -r -s -@ 7 -f "${work}/markdup.txt" \
        "${work}/coordinates.bam" "${work}/${sample}.dedup.bam"
    rm -f -- "${work}/coordinates.bam"
    "${tools}/samtools" index -@ 7 "${work}/${sample}.dedup.bam"
    "${tools}/samtools" quickcheck -v "${work}/${sample}.dedup.bam"
    "${tools}/samtools" flagstat -@ 7 "${work}/${sample}.dedup.bam" > "${work}/flagstat.txt"

    "${python_bin}" - "${work}/${sample}.dedup.bam" "${reference}" "${sample}" "${work}/alignment.complete.json" <<'PY'
import json
from pathlib import Path
import sys
import pysam
bam, reference, sample, audit = sys.argv[1:]
with open(reference + '.fai') as handle:
    expected = {f[0]: int(f[1]) for f in (line.split('\t') for line in handle)}
with pysam.AlignmentFile(bam, 'rb') as handle:
    if dict(zip(handle.references, handle.lengths)) != expected:
        raise RuntimeError('Reference dictionary mismatch')
    groups = handle.header.to_dict().get('RG', [])
    if not groups or any(g.get('SM') != sample or g.get('ID') != sample for g in groups):
        raise RuntimeError('Read-group validation failed')
    handle.check_index()
    if handle.mapped == 0:
        raise RuntimeError('No mapped reads')
    result = {'sample': sample, 'reference': reference, 'mapped_alignments': handle.mapped,
              'bam_filename': Path(bam).name, 'duplicate_policy': 'samtools markdup -r'}
Path(audit).write_text(json.dumps(result, indent=2) + '\n')
print('BAM_VALIDATION_PASS\t' + sample, flush=True)
PY
    mv -T "${work}" "${final}"
    work=""
    printf 'BAM\t%s/%s.dedup.bam\n' "${final}" "${sample}"
fi
