#!/bin/bash
#SBATCH --job-name=smudgeStd
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-6%2
#SBATCH --output=logs/smudgeStd.%A_%a.out
#SBATCH --error=logs/smudgeStd.%A_%a.err

set -euo pipefail

PROJECT="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
BASE="$PROJECT/analysis/smudgeplot_standardised"

mkdir -p "$BASE"

export PATH="$HOME/.local/share/mamba/envs/smudgeplot/bin:$PATH"

PYTHON="$PROJECT/pairwise-af/.venv-af/bin/python"
PLOTTER="$PROJECT/scripts/plot_smudge_smu_viridis.py"

FASTK="$(command -v FastK || true)"
SMUDGE="$(command -v smudgeplot || true)"

resolve_tool () {
    local name="$1"
    local hit

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
        echo "ERROR: cannot find $name" >&2
        exit 1
    fi

    printf '%s\n' "$hit"
}

SEQTK="$(resolve_tool seqtk)"
FASTERQ="$(resolve_tool fasterq-dump)"

if [[ -z "$FASTK" ]]; then
    echo "ERROR: FastK not found" >&2
    exit 1
fi

if [[ -z "$SMUDGE" ]]; then
    echo "ERROR: smudgeplot not found" >&2
    exit 1
fi

if [[ ! -x "$PYTHON" ]]; then
    echo "ERROR: Python missing: $PYTHON" >&2
    exit 1
fi

if [[ ! -f "$PLOTTER" ]]; then
    echo "ERROR: plotter missing: $PLOTTER" >&2
    exit 1
fi


# ============================================================
# Gynandromorph source FASTQs
# ============================================================

FASTQ_ROOT="$PROJECT/analysis/genomewide_mixing_input/fastq"

T2L_R1="$FASTQ_ROOT/TN_Left_T2/TN_Left_T2_1.fastq.gz"
T2L_R2="$FASTQ_ROOT/TN_Left_T2/TN_Left_T2_2.fastq.gz"

T3L_R1="$FASTQ_ROOT/TN_Left_T3/TN_Left_T3_1.fastq.gz"
T3L_R2="$FASTQ_ROOT/TN_Left_T3/TN_Left_T3_2.fastq.gz"

T2R_R1="$FASTQ_ROOT/TN_Right_T2/TN_Right_T2_1.fastq.gz"
T2R_R2="$FASTQ_ROOT/TN_Right_T2/TN_Right_T2_2.fastq.gz"

T3R_R1="$FASTQ_ROOT/TN_Right_T3/TN_Right_T3_1.fastq.gz"
T3R_R2="$FASTQ_ROOT/TN_Right_T3/TN_Right_T3_2.fastq.gz"


# ============================================================
# Resolve an ERR accession.
#
# Prefer existing final FASTQ. If unavailable, extract from the
# local SRA into THIS TASK'S directory, avoiding array races.
# ============================================================

resolve_accession () {
    local acc="$1"
    local task_source="$2"

    local r1=""
    local r2=""

    mapfile -t r1hits < <(
        find "$PROJECT/analysis/genomewide_mixing_input" \
            -type f \
            \( \
                -name "${acc}_1.fastq" \
                -o -name "${acc}_1.fastq.gz" \
            \) \
            ! -path '*/tmp/*' \
            2>/dev/null \
            | sort -u
    )

    mapfile -t r2hits < <(
        find "$PROJECT/analysis/genomewide_mixing_input" \
            -type f \
            \( \
                -name "${acc}_2.fastq" \
                -o -name "${acc}_2.fastq.gz" \
            \) \
            ! -path '*/tmp/*' \
            2>/dev/null \
            | sort -u
    )

    if [[ "${#r1hits[@]}" -eq 1 && "${#r2hits[@]}" -eq 1 ]]; then
        printf '%s\n%s\n' \
            "${r1hits[0]}" \
            "${r2hits[0]}"
        return
    fi

    mkdir -p "$task_source/$acc"

    r1="$task_source/$acc/${acc}_1.fastq"
    r2="$task_source/$acc/${acc}_2.fastq"

    if [[ ! -s "$r1" || ! -s "$r2" ]]; then

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
                "ERROR: cannot uniquely resolve local SRA for $acc" \
                >&2
            exit 1
        fi

        tmp="$task_source/$acc/fasterq_tmp"

        rm -rf "$tmp"
        mkdir -p "$tmp"

        "$FASTERQ" \
            --split-files \
            --threads 16 \
            --temp "$tmp" \
            --outdir "$task_source/$acc" \
            "${sras[0]}"

        rm -rf "$tmp"
    fi

    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "ERROR: failed to obtain FASTQs for $acc" >&2
        exit 1
    fi

    printf '%s\n%s\n' "$r1" "$r2"
}


# ============================================================
# Panel specification
#
# Every task receives EXACTLY 20M paired fragments total.
# ============================================================

case "$SLURM_ARRAY_TASK_ID" in

    0)
        NAME="01_one_diploid"
        TITLE="One diploid"
        SOURCES=("ERR2298208")
        PAIRS=(20000000)
        ;;

    1)
        NAME="02_two_diploids"
        TITLE="Two distinct diploids"
        SOURCES=("ERR2298208" "ERR1307111")
        PAIRS=(10000000 10000000)
        ;;

    2)
        NAME="03_four_leg_mix"
        TITLE="Four-leg mix"
        SOURCES=("T2L" "T3L" "T2R" "T3R")
        PAIRS=(5000000 5000000 5000000 5000000)
        ;;

    3)
        NAME="04_ipsi_left"
        TITLE="Ipsi-left"
        SOURCES=("T2L" "T3L")
        PAIRS=(10000000 10000000)
        ;;

    4)
        NAME="05_ipsi_right"
        TITLE="Ipsi-right"
        SOURCES=("T2R" "T3R")
        PAIRS=(10000000 10000000)
        ;;

    5)
        NAME="06_contra_T2"
        TITLE="Contra-T2"
        SOURCES=("T2L" "T2R")
        PAIRS=(10000000 10000000)
        ;;

    6)
        NAME="07_contra_T3"
        TITLE="Contra-T3"
        SOURCES=("T3L" "T3R")
        PAIRS=(10000000 10000000)
        ;;

    *)
        echo "ERROR: unknown task $SLURM_ARRAY_TASK_ID" >&2
        exit 1
        ;;
esac


TASK="$BASE/$NAME"
SOURCE_TMP="$TASK/source_tmp"
SUB="$TASK/subsampled"

mkdir -p \
    "$TASK" \
    "$SOURCE_TMP" \
    "$SUB"

echo "TASK=$SLURM_ARRAY_TASK_ID"
echo "NAME=$NAME"
echo "TITLE=$TITLE"
echo "SOURCES=${SOURCES[*]}"
echo "PAIRS=${PAIRS[*]}"


# ============================================================
# Resolve named source to R1/R2
# ============================================================

get_pair () {

    local source="$1"

    case "$source" in

        T2L)
            printf '%s\n%s\n' \
                "$T2L_R1" "$T2L_R2"
            ;;

        T3L)
            printf '%s\n%s\n' \
                "$T3L_R1" "$T3L_R2"
            ;;

        T2R)
            printf '%s\n%s\n' \
                "$T2R_R1" "$T2R_R2"
            ;;

        T3R)
            printf '%s\n%s\n' \
                "$T3R_R1" "$T3R_R2"
            ;;

        ERR*)
            resolve_accession \
                "$source" \
                "$SOURCE_TMP"
            ;;

        *)
            echo "ERROR: unknown source $source" >&2
            exit 1
            ;;
    esac
}


# ============================================================
# Deterministic paired subsampling
# ============================================================

FASTK_INPUTS=()

MANIFEST="$TASK/input_manifest.tsv"

printf \
    'panel\tsource\trequested_pairs\tseed\tR1\tR2\n' \
    > "$MANIFEST"

for index in "${!SOURCES[@]}"; do

    source="${SOURCES[$index]}"
    npairs="${PAIRS[$index]}"

    seed="$((1000 + SLURM_ARRAY_TASK_ID * 100 + index))"

    mapfile -t pair < <(
        get_pair "$source"
    )

    if [[ "${#pair[@]}" -ne 2 ]]; then
        echo "ERROR: failed to resolve $source" >&2
        exit 1
    fi

    r1="${pair[0]}"
    r2="${pair[1]}"

    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "ERROR: missing input for $source" >&2
        echo "R1=$r1" >&2
        echo "R2=$r2" >&2
        exit 1
    fi

    out1="$SUB/${index}_${source}_1.fastq"
    out2="$SUB/${index}_${source}_2.fastq"

    echo "SUBSAMPLE $source : $npairs pairs"

    "$SEQTK" sample \
        -s"$seed" \
        "$r1" \
        "$npairs" \
        > "${out1}.partial"

    "$SEQTK" sample \
        -s"$seed" \
        "$r2" \
        "$npairs" \
        > "${out2}.partial"

    mv "${out1}.partial" "$out1"
    mv "${out2}.partial" "$out2"

    # Exact read-count validation.
    expected_lines="$((npairs * 4))"

    lines1="$(wc -l < "$out1")"
    lines2="$(wc -l < "$out2")"

    if [[ "$lines1" -ne "$expected_lines" ]]; then
        echo \
            "ERROR: $source R1 has $lines1 lines; expected $expected_lines" \
            >&2
        exit 1
    fi

    if [[ "$lines2" -ne "$expected_lines" ]]; then
        echo \
            "ERROR: $source R2 has $lines2 lines; expected $expected_lines" \
            >&2
        exit 1
    fi

    printf \
        '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$NAME" \
        "$source" \
        "$npairs" \
        "$seed" \
        "$r1" \
        "$r2" \
        >> "$MANIFEST"

    FASTK_INPUTS+=(
        "$out1"
        "$out2"
    )
done


# ============================================================
# Confirm total standardized depth
# ============================================================

total_pairs=0

for n in "${PAIRS[@]}"; do
    total_pairs="$((total_pairs + n))"
done

if [[ "$total_pairs" -ne 20000000 ]]; then
    echo \
        "ERROR: panel total is $total_pairs, not 20000000" \
        >&2
    exit 1
fi

echo "TOTAL_PAIRED_FRAGMENTS=$total_pairs"


# ============================================================
# FastK
#
# Same biological settings as recovered old analysis:
#   k = 22
#   -t4
#   -M100
#   -T60
# ============================================================

if [[ ! -s "$TASK/kmerpairs.smu" ]]; then

    echo
    echo "=== FASTK ==="

    "$FASTK" \
        -v \
        -t4 \
        -k22 \
        -M100 \
        -T60 \
        "${FASTK_INPUTS[@]}" \
        -N"$TASK/FastK_Table"

    echo
    echo "=== SMUDGEPLOT HETMERS ==="

    "$SMUDGE" hetmers \
        -L 12 \
        -t 60 \
        -o "$TASK/kmerpairs" \
        --verbose \
        "$TASK/FastK_Table"

else
    echo "REUSING $TASK/kmerpairs.smu"
fi


if [[ ! -s "$TASK/kmerpairs.smu" ]]; then
    echo "ERROR: missing kmerpairs.smu" >&2
    exit 1
fi


# ============================================================
# Standard smudgeplot outputs
# ============================================================

echo
echo "=== SMUDGEPLOT ALL ==="

"$SMUDGE" all \
    -o "$TASK/output" \
    "$TASK/kmerpairs.smu"


# ============================================================
# Vector VIRIDIS rendering
# ============================================================

echo
echo "=== VECTOR VIRIDIS PLOT ==="

"$PYTHON" \
    "$PLOTTER" \
    --input "$TASK/kmerpairs.smu" \
    --output-prefix "$TASK/${NAME}_viridis" \
    --title "$TITLE"


# ============================================================
# Remove huge temporary FASTQs after successful analysis
# ============================================================

rm -rf \
    "$SUB" \
    "$SOURCE_TMP"


# ============================================================
# Audit
# ============================================================

echo
echo "=== OUTPUTS ==="

find "$TASK" \
    -maxdepth 1 \
    -type f \
    \( \
        -name '*.smu' \
        -o -name '*.smudge_report.tsv' \
        -o -name '*.pdf' \
        -o -name '*.svg' \
        -o -name '*.png' \
        -o -name 'input_manifest.tsv' \
    \) \
    -printf '%f\t%s bytes\n' \
    | sort

echo
echo "COMPLETE $NAME"
