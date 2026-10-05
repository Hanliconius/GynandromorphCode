#!/usr/bin/env bash
set -euo pipefail

ROOT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
cd "${ROOT}"

OUTROOT="${ROOT}/analysis/smudgeplot_arrangements"
LAUNCHER="${ROOT}/external/gynandromorph_2024/smudge_launch.sh"

mkdir -p "${OUTROOT}"

# ------------------------------------------------------------
# Smudgeplot environment
# ------------------------------------------------------------

SMUDGE_ENV="smudgeplot"

# Sanity check before doing anything expensive
micromamba run -n "${SMUDGE_ENV}" bash -c '
    echo "FastK:"
    command -v FastK
    FastK 2>&1 | head -n 2 || true

    echo
    echo "smudgeplot:"
    command -v smudgeplot
    smudgeplot --version 2>&1 || true
'

# ------------------------------------------------------------
# IMPORTANT:
# We have not yet filled in the real FASTQ paths here.
# Do not run real jobs until these are replaced.
# ------------------------------------------------------------

TN_LEFT_T2_R1="/full/path/to/TN_Left_T2_R1.fastq.gz"
TN_LEFT_T2_R2="/full/path/to/TN_Left_T2_R2.fastq.gz"

TN_LEFT_T3_R1="/full/path/to/TN_Left_T3_R1.fastq.gz"
TN_LEFT_T3_R2="/full/path/to/TN_Left_T3_R2.fastq.gz"

TN_RIGHT_T2_R1="/full/path/to/TN_Right_T2_R1.fastq.gz"
TN_RIGHT_T2_R2="/full/path/to/TN_Right_T2_R2.fastq.gz"

TN_RIGHT_T3_R1="/full/path/to/TN_Right_T3_R1.fastq.gz"
TN_RIGHT_T3_R2="/full/path/to/TN_Right_T3_R2.fastq.gz"

CANDIDATE1_NAME="ERR2298208"
CANDIDATE1_R1="/full/path/to/ERR2298208_R1.fastq.gz"
CANDIDATE1_R2="/full/path/to/ERR2298208_R2.fastq.gz"

CANDIDATE2_NAME="SECOND_CONTROL"
CANDIDATE2_R1="/full/path/to/SECOND_CONTROL_R1.fastq.gz"
CANDIDATE2_R2="/full/path/to/SECOND_CONTROL_R2.fastq.gz"

MANIFEST="${OUTROOT}/smudgeplot_arrangements.manifest.tsv"

# ------------------------------------------------------------
# Write REAL tabs with printf
# ------------------------------------------------------------

{
    printf "label\tr1_csv\tr2_csv\n"

    printf "TN_Left_merge\t%s,%s\t%s,%s\n" \
        "${TN_LEFT_T2_R1}" "${TN_LEFT_T3_R1}" \
        "${TN_LEFT_T2_R2}" "${TN_LEFT_T3_R2}"

    printf "TN_Right_merge\t%s,%s\t%s,%s\n" \
        "${TN_RIGHT_T2_R1}" "${TN_RIGHT_T3_R1}" \
        "${TN_RIGHT_T2_R2}" "${TN_RIGHT_T3_R2}"

    printf "TN_T2_merge\t%s,%s\t%s,%s\n" \
        "${TN_LEFT_T2_R1}" "${TN_RIGHT_T2_R1}" \
        "${TN_LEFT_T2_R2}" "${TN_RIGHT_T2_R2}"

    printf "TN_T3_merge\t%s,%s\t%s,%s\n" \
        "${TN_LEFT_T3_R1}" "${TN_RIGHT_T3_R1}" \
        "${TN_LEFT_T3_R2}" "${TN_RIGHT_T3_R2}"

    printf "%s_single\t%s\t%s\n" \
        "${CANDIDATE1_NAME}" \
        "${CANDIDATE1_R1}" \
        "${CANDIDATE1_R2}"

    printf "pseudoTetraploid_%s_%s\t%s,%s\t%s,%s\n" \
        "${CANDIDATE1_NAME}" "${CANDIDATE2_NAME}" \
        "${CANDIDATE1_R1}" "${CANDIDATE2_R1}" \
        "${CANDIDATE1_R2}" "${CANDIDATE2_R2}"

} > "${MANIFEST}"

echo
echo "MANIFEST:"
column -t -s $'\t' "${MANIFEST}"
echo

# ------------------------------------------------------------
# HARD STOP while placeholder paths remain
# ------------------------------------------------------------

if grep -q '/full/path/to/' "${MANIFEST}"; then
    echo "STOP: manifest still contains placeholder FASTQ paths."
    echo "We need to replace those with the real existing files first."
    exit 1
fi

# ------------------------------------------------------------
# Runner
#
# The original launcher itself runs:
#
# FastK -v -t4 -k22 -M100 -T60 ...
# smudgeplot hetmers ...
# smudgeplot all ...
#
# ------------------------------------------------------------

run_one() {
    local label="$1"
    local r1_csv="$2"
    local r2_csv="$3"

    local outdir="${OUTROOT}/${label}"
    mkdir -p "${outdir}"

    echo "RUN     ${label}"
    echo "R1      ${r1_csv}"
    echo "R2      ${r2_csv}"
    echo "OUTDIR  ${outdir}"

    # We will wire the exact FASTQ handling into this once the
    # real paths are filled in.
    micromamba run -n "${SMUDGE_ENV}" \
        bash "${LAUNCHER}" \
        "${label}" \
        "${outdir}" \
        "${r1_csv}" \
        "${r2_csv}"
}

tail -n +2 "${MANIFEST}" |
while IFS=$'\t' read -r label r1_csv r2_csv; do
    run_one "${label}" "${r1_csv}" "${r2_csv}"
done

