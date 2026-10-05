#!/bin/bash
#SBATCH --job-name=smudgeDipCtrl
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/smudgeDipCtrl.%j.out
#SBATCH --error=logs/smudgeDipCtrl.%j.err

set -euo pipefail

PROJECT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

A="ERR2298208"
B="ERR1307111"

# Depth-controlled design:
#   single diploid = 20 M pairs of A
#   two diploids   = 10 M pairs of A + 10 M pairs of B
SINGLE_PAIRS=20000000
PAIR_EACH=10000000

WORK="$PROJECT/analysis/smudgeplot_diploid_benchmarks"

mkdir -p \
    "$WORK/raw" \
    "$WORK/subsampled" \
    "$WORK/single_ERR2298208" \
    "$WORK/two_diploids_ERR2298208_ERR1307111"


# ============================================================
# Resolve tools
# ============================================================

resolve_tool () {
    local name="$1"
    local hit=""

    hit="$(command -v "$name" 2>/dev/null || true)"

    if [[ -z "$hit" ]]; then
        hit="$(
            find "$HOME/.local/share/mamba/envs" \
                -type f \
                -path "*/bin/$name" \
                -perm -u+x \
                -print 2>/dev/null \
                | head -1
        )"
    fi

    if [[ -z "$hit" ]]; then
        echo "ERROR: could not find $name" >&2
        exit 1
    fi

    printf '%s\n' "$hit"
}

FASTK="$(resolve_tool FastK)"
SMUDGEPLOT="$(resolve_tool smudgeplot)"
SEQTK="$(resolve_tool seqtk)"
FASTERQ="$(resolve_tool fasterq-dump)"

PYTHON="$PROJECT/pairwise-af/.venv-af/bin/python"

if [[ ! -x "$PYTHON" ]]; then
    echo "ERROR: missing plotting Python: $PYTHON" >&2
    exit 1
fi

echo "FASTK=$FASTK"
echo "SMUDGEPLOT=$SMUDGEPLOT"
echo "SEQTK=$SEQTK"
echo "FASTERQ=$FASTERQ"
echo "PYTHON=$PYTHON"


# ============================================================
# Resolve existing FASTQs, or regenerate from the already-local
# SRA file if final FASTQs are not available.
# ============================================================

resolve_existing_fastq () {
    local acc="$1"
    local mate="$2"

    mapfile -t hits < <(
        find \
            "$PROJECT/analysis/genomewide_mixing_input" \
            -type f \
            \( \
                -name "${acc}_${mate}.fastq" \
                -o -name "${acc}_${mate}.fastq.gz" \
                -o -name "${acc}_${mate}.fq" \
                -o -name "${acc}_${mate}.fq.gz" \
            \) \
            ! -path '*/tmp/*' \
            2>/dev/null \
            | sort -u
    )

    if [[ "${#hits[@]}" -eq 1 ]]; then
        printf '%s\n' "${hits[0]}"
        return 0
    fi

    if [[ "${#hits[@]}" -gt 1 ]]; then
        echo \
            "ERROR: multiple non-temporary FASTQs for ${acc} mate ${mate}" \
            >&2

        printf '%s\n' "${hits[@]}" >&2
        exit 1
    fi

    return 1
}


prepare_accession () {
    local acc="$1"

    local existing1=""
    local existing2=""

    existing1="$(
        resolve_existing_fastq "$acc" 1 || true
    )"

    existing2="$(
        resolve_existing_fastq "$acc" 2 || true
    )"

    if [[ -n "$existing1" && -n "$existing2" ]]; then
        echo "REUSING_FASTQ $acc" >&2

        printf '%s\n' \
            "$existing1" \
            "$existing2"

        return
    fi

    local rawdir="$WORK/raw/$acc"
    mkdir -p "$rawdir"

    local r1="$rawdir/${acc}_1.fastq"
    local r2="$rawdir/${acc}_2.fastq"

    if [[ -s "$r1" && -s "$r2" ]]; then
        echo "REUSING_EXTRACTED_FASTQ $acc" >&2

        printf '%s\n' \
            "$r1" \
            "$r2"

        return
    fi

    mapfile -t sras < <(
        find \
            "$PROJECT/analysis/genomewide_mixing_input/sra/$acc" \
            -type f \
            -name "${acc}.sra" \
            2>/dev/null \
            | sort -u
    )

    if [[ "${#sras[@]}" -ne 1 ]]; then
        echo \
            "ERROR: expected exactly one local SRA for $acc; found ${#sras[@]}" \
            >&2

        printf '%s\n' "${sras[@]}" >&2
        exit 1
    fi

    echo "EXTRACTING $acc FROM ${sras[0]}" >&2

    rm -f \
        "$rawdir/${acc}_1.fastq.partial" \
        "$rawdir/${acc}_2.fastq.partial"

    tmpdir="$rawdir/fasterq_tmp"
    rm -rf "$tmpdir"
    mkdir -p "$tmpdir"

    "$FASTERQ" \
        --split-files \
        --threads 16 \
        --temp "$tmpdir" \
        --outdir "$rawdir" \
        "${sras[0]}"

    rm -rf "$tmpdir"

    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "ERROR: fasterq-dump did not make paired FASTQs for $acc" >&2
        exit 1
    fi

    printf '%s\n' \
        "$r1" \
        "$r2"
}


mapfile -t A_FASTQ < <(
    prepare_accession "$A"
)

mapfile -t B_FASTQ < <(
    prepare_accession "$B"
)

if [[ "${#A_FASTQ[@]}" -ne 2 ]]; then
    echo "ERROR: failed to resolve paired FASTQs for $A" >&2
    exit 1
fi

if [[ "${#B_FASTQ[@]}" -ne 2 ]]; then
    echo "ERROR: failed to resolve paired FASTQs for $B" >&2
    exit 1
fi

A_R1="${A_FASTQ[0]}"
A_R2="${A_FASTQ[1]}"

B_R1="${B_FASTQ[0]}"
B_R2="${B_FASTQ[1]}"

echo
echo "INPUTS"
echo "A_R1=$A_R1"
echo "A_R2=$A_R2"
echo "B_R1=$B_R1"
echo "B_R2=$B_R2"


# ============================================================
# Deterministic paired subsampling.
#
# seqtk uses the same seed separately on R1 and R2, preserving
# paired record selection because the mate files have the same
# record ordering.
# ============================================================

sample_fastq () {
    local input="$1"
    local output="$2"
    local number="$3"
    local seed="$4"

    if [[ -s "$output" ]]; then
        echo "REUSING_SUBSAMPLE $output"
        return
    fi

    local partial="${output}.partial"

    rm -f "$partial"

    "$SEQTK" sample \
        -s"$seed" \
        "$input" \
        "$number" \
        > "$partial"

    if [[ ! -s "$partial" ]]; then
        echo "ERROR: empty seqtk output: $partial" >&2
        exit 1
    fi

    mv "$partial" "$output"
}


# Single diploid:
# 20M paired fragments from ERR2298208.

SINGLE_R1="$WORK/subsampled/${A}.single20m_1.fastq"
SINGLE_R2="$WORK/subsampled/${A}.single20m_2.fastq"

sample_fastq \
    "$A_R1" \
    "$SINGLE_R1" \
    "$SINGLE_PAIRS" \
    101

sample_fastq \
    "$A_R2" \
    "$SINGLE_R2" \
    "$SINGLE_PAIRS" \
    101


# Two-diploid pool:
# 10M paired fragments from each individual.

PAIR_A_R1="$WORK/subsampled/${A}.pair10m_1.fastq"
PAIR_A_R2="$WORK/subsampled/${A}.pair10m_2.fastq"

PAIR_B_R1="$WORK/subsampled/${B}.pair10m_1.fastq"
PAIR_B_R2="$WORK/subsampled/${B}.pair10m_2.fastq"

sample_fastq \
    "$A_R1" \
    "$PAIR_A_R1" \
    "$PAIR_EACH" \
    202

sample_fastq \
    "$A_R2" \
    "$PAIR_A_R2" \
    "$PAIR_EACH" \
    202

sample_fastq \
    "$B_R1" \
    "$PAIR_B_R1" \
    "$PAIR_EACH" \
    303

sample_fastq \
    "$B_R2" \
    "$PAIR_B_R2" \
    "$PAIR_EACH" \
    303


# ============================================================
# Verify sampled FASTQ sizes exactly.
# ============================================================

check_fastq () {
    local file="$1"
    local expected="$2"

    lines="$(
        wc -l < "$file"
    )"

    expected_lines="$(
        awk \
            -v n="$expected" \
            'BEGIN { printf "%.0f", n * 4 }'
    )"

    if [[ "$lines" != "$expected_lines" ]]; then
        echo \
            "ERROR: $file has $lines lines; expected $expected_lines" \
            >&2
        exit 1
    fi
}

check_fastq "$SINGLE_R1" "$SINGLE_PAIRS"
check_fastq "$SINGLE_R2" "$SINGLE_PAIRS"

check_fastq "$PAIR_A_R1" "$PAIR_EACH"
check_fastq "$PAIR_A_R2" "$PAIR_EACH"

check_fastq "$PAIR_B_R1" "$PAIR_EACH"
check_fastq "$PAIR_B_R2" "$PAIR_EACH"

echo "SUBSAMPLING_PASS"


# ============================================================
# Same FastK / smudgeplot workflow as the old analysis.
# ============================================================

run_smudge () {
    local out="$1"
    shift

    mkdir -p "$out"

    if [[ ! -s "$out/kmerpairs.smu" ]]; then

        # Remove possible partial previous FastK products.
        rm -f \
            "$out"/FastK_Table.* \
            "$out"/kmerpairs.* \
            2>/dev/null || true

        echo
        echo "FASTK $out"

        "$FASTK" \
            -v \
            -t4 \
            -k22 \
            -M100 \
            -T60 \
            "$@" \
            -N"$out/FastK_Table"

        echo
        echo "SMUDGEPLOT HETMERS $out"

        "$SMUDGEPLOT" hetmers \
            -L 12 \
            -t 60 \
            -o "$out/kmerpairs" \
            --verbose \
            "$out/FastK_Table"

    else
        echo "REUSING_SMU $out/kmerpairs.smu"
    fi

    if [[ ! -s "$out/output.smudge_report.tsv" ]]; then

        echo
        echo "SMUDGEPLOT ALL $out"

        "$SMUDGEPLOT" all \
            -o "$out/output" \
            "$out/kmerpairs.smu"
    fi
}


run_smudge \
    "$WORK/single_ERR2298208" \
    "$SINGLE_R1" \
    "$SINGLE_R2"

run_smudge \
    "$WORK/two_diploids_ERR2298208_ERR1307111" \
    "$PAIR_A_R1" \
    "$PAIR_A_R2" \
    "$PAIR_B_R1" \
    "$PAIR_B_R2"


# ============================================================
# Vector luka_bleute rendering
# ============================================================

"$PYTHON" \
    "$PROJECT/scripts/plot_smudge_smu_luka_bleute.py" \
    --input \
    "$WORK/single_ERR2298208/kmerpairs.smu" \
    --output-prefix \
    "$WORK/single_ERR2298208_luka_bleute" \
    --title \
    "One diploid"

"$PYTHON" \
    "$PROJECT/scripts/plot_smudge_smu_luka_bleute.py" \
    --input \
    "$WORK/two_diploids_ERR2298208_ERR1307111/kmerpairs.smu" \
    --output-prefix \
    "$WORK/two_diploids_ERR2298208_ERR1307111_luka_bleute" \
    --title \
    "Two distinct diploids"


# ============================================================
# Final audit
# ============================================================

echo
echo "FINAL OUTPUTS"

find "$WORK" \
    -maxdepth 2 \
    -type f \
    \( \
        -name '*.smu' \
        -o -name '*.smudge_report.tsv' \
        -o -name '*luka_bleute.pdf' \
        -o -name '*luka_bleute.svg' \
    \) \
    -printf '%p\t%s bytes\n' \
    | sort

echo
echo "COMPLETE"
