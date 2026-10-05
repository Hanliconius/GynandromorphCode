#!/usr/bin/env bash
#SBATCH --job-name=gyn166rc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/166c_prepare_correct_rosina_reference.%j.out
#SBATCH --error=logs/166c_prepare_correct_rosina_reference.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
root="${project}/analysis/genomewide_mixing_input"
reference_root="${root}/reference"
metadata_dir="${root}/metadata"
reference_dir="${reference_root}/GCF_056956725.1"

tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
bwa_bin="${tool_bin}/bwa"
samtools_bin="${tool_bin}/samtools"

curl_candidates=(
    /usr/bin/curl
    "${tool_bin}/curl"
)

accession=GCF_056956725.1
assembly=GCF_056956725.1_ASM5695672v1
base_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/056/956/725/${assembly}"

fasta_gz="${reference_dir}/${assembly}_genomic.fna.gz"
fasta="${reference_dir}/${assembly}_genomic.fna"
report="${reference_dir}/${assembly}_assembly_report.txt"
dict_file="${fasta%.fna}.dict"
complete_file="${reference_dir}/reference.complete"

cd "${project}"

mkdir -p "${reference_dir}" "${metadata_dir}"

download_file() {
    local url=$1
    local destination=$2
    local partial="${destination}.part"
    local downloaded=0
    local curl_bin
    local attempt

    if [[ -s "${destination}" ]]
    then
        printf 'Reusing existing download: %s\n' "${destination}"
        return
    fi

    for curl_bin in "${curl_candidates[@]}"
    do
        if [[ ! -x "${curl_bin}" ]]
        then
            continue
        fi

        printf 'Trying downloader: %s\n' "${curl_bin}"

        for attempt in 1 2 3 4 5 6
        do
            printf 'Download attempt %s: %s\n' \
                "${attempt}" \
                "${url}"

            if "${curl_bin}" \
                --fail \
                --location \
                --http1.1 \
                --connect-timeout 60 \
                --continue-at - \
                --output "${partial}" \
                "${url}"
            then
                downloaded=1
                break 2
            fi
        done
    done

    if [[ "${downloaded}" -ne 1 || ! -s "${partial}" ]]
    then
        printf 'ERROR: download failed or was empty: %s\n' \
            "${url}" >&2
        false
    fi

    mv -- "${partial}" "${destination}"
}

download_file \
    "${base_url}/${assembly}_genomic.fna.gz" \
    "${fasta_gz}"

download_file \
    "${base_url}/${assembly}_assembly_report.txt" \
    "${report}"

gzip -t "${fasta_gz}"

if ! grep -q '^# Assembly level: Chromosome' "${report}"
then
    printf 'ERROR: NCBI report does not identify this as a chromosome-level assembly.\n' >&2
    grep '^# Assembly level:' "${report}" >&2
    false
fi

if [[ ! -s "${fasta}" ]]
then
    fasta_tmp="${fasta}.tmp.${SLURM_JOB_ID}"

    gzip -dc "${fasta_gz}" > "${fasta_tmp}"

    mv -- "${fasta_tmp}" "${fasta}"
fi

if [[ ! -s "${fasta}.fai" ]]
then
    "${samtools_bin}" faidx "${fasta}"
fi

sequence_count=$(wc -l < "${fasta}.fai")

total_bases=$(
    awk '
        {
            total += $2
        }

        END {
            printf "%.0f", total
        }
    ' "${fasta}.fai"
)

large_sequences=$(
    awk '
        $2 >= 5000000 {
            count++
        }

        END {
            print count + 0
        }
    ' "${fasta}.fai"
)

if [[ "${large_sequences}" -lt 20 ]]
then
    printf 'ERROR: only %s sequences are at least 5 Mb; expected at least 20.\n' \
        "${large_sequences}" >&2
    false
fi

if [[ ! -s "${fasta}.bwt" ]]
then
    "${bwa_bin}" index "${fasta}"
fi

if [[ ! -s "${dict_file}" ]]
then
    "${samtools_bin}" dict \
        -o "${dict_file}" \
        "${fasta}"
fi

for required in \
    "${fasta}" \
    "${fasta}.fai" \
    "${fasta}.amb" \
    "${fasta}.ann" \
    "${fasta}.bwt" \
    "${fasta}.pac" \
    "${fasta}.sa" \
    "${dict_file}" \
    "${report}"
do
    if [[ ! -s "${required}" ]]
    then
        printf 'ERROR: missing or empty reference product: %s\n' \
            "${required}" >&2
        false
    fi
done

if [[ \
    -s "${metadata_dir}/reference_fasta.path" \
    && ! -e "${metadata_dir}/reference_fasta.fragmented_GCA_900068215.1.path" \
]]
then
    cp -p \
        "${metadata_dir}/reference_fasta.path" \
        "${metadata_dir}/reference_fasta.fragmented_GCA_900068215.1.path"
fi

if [[ \
    -s "${metadata_dir}/reference_manifest.tsv" \
    && ! -e "${metadata_dir}/reference_manifest.fragmented_GCA_900068215.1.tsv" \
]]
then
    cp -p \
        "${metadata_dir}/reference_manifest.tsv" \
        "${metadata_dir}/reference_manifest.fragmented_GCA_900068215.1.tsv"
fi

path_tmp="${metadata_dir}/reference_fasta.path.tmp.${SLURM_JOB_ID}"

printf '%s\n' "${fasta}" > "${path_tmp}"

mv -- \
    "${path_tmp}" \
    "${metadata_dir}/reference_fasta.path"

manifest_tmp="${metadata_dir}/reference_manifest.tsv.tmp.${SLURM_JOB_ID}"

printf '%s\n' \
    $'accession\tassembly_name\tassembly_level\tsource_url\tfasta' \
    > "${manifest_tmp}"

printf '%s\t%s\tChromosome\t%s\t%s\n' \
    "${accession}" \
    ASM5695672v1 \
    "${base_url}" \
    "${fasta}" \
    >> "${manifest_tmp}"

mv -- \
    "${manifest_tmp}" \
    "${metadata_dir}/reference_manifest.tsv"

sha256sum \
    "${fasta_gz}" \
    "${fasta}" \
    "${report}" \
    > "${reference_dir}/reference_files.sha256"

printf 'assembly_accession\t%s\nassembly_name\tASM5695672v1\nassembly_level\tChromosome\nfasta\t%s\nsequences\t%s\nbases\t%s\nsequences_ge_5Mb\t%s\ncompleted_at\t%s\n' \
    "${accession}" \
    "${fasta}" \
    "${sequence_count}" \
    "${total_bases}" \
    "${large_sequences}" \
    "$(date --iso-8601=seconds)" \
    > "${complete_file}"

printf 'CORRECT ROSINA REFERENCE COMPLETE\n'

cat "${complete_file}"

printf '\nTWENTY-TWO LONGEST REFERENCE SEQUENCES\n'

sort -k2,2nr "${fasta}.fai" \
    | sed -n '1,22p'
