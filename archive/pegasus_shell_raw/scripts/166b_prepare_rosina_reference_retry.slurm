#!/usr/bin/env bash
#SBATCH --job-name=gyn166rf
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/166b_prepare_rosina_reference_retry.%j.out
#SBATCH --error=logs/166b_prepare_rosina_reference_retry.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
root="${project}/analysis/genomewide_mixing_input"
reference_dir="${root}/reference/GCA_900068215.1"
metadata_dir="${root}/metadata"

tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
bwa_bin="${tool_bin}/bwa"
samtools_bin="${tool_bin}/samtools"

curl_candidates=(
    /usr/bin/curl
    "${tool_bin}/curl"
)

assembly=GCA_900068215.1_Heliconius_melpomene_rosina.546
base_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/068/215/${assembly}"

fasta_gz="${reference_dir}/${assembly}_genomic.fna.gz"
fasta="${reference_dir}/${assembly}_genomic.fna"
report="${reference_dir}/${assembly}_assembly_report.txt"
complete_file="${reference_dir}/reference.complete"

cd "${project}"

mkdir -p "${reference_dir}" "${metadata_dir}"

if [[ \
    -s "${complete_file}" \
    && -s "${fasta}" \
    && -s "${fasta}.fai" \
    && -s "${fasta}.bwt" \
    && -s "${fasta%.fna}.dict" \
]]
then
    printf 'Reference already complete: %s\n' "${fasta}"
    exit 0
fi

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
        "${curl_bin}" --version

        for attempt in 1 2 3 4 5 6
        do
            printf 'Download attempt %s with %s\n' \
                "${attempt}" \
                "${curl_bin}"

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

    if [[ "${downloaded}" -ne 1 ]]
    then
        printf 'ERROR: every download attempt failed: %s\n' \
            "${url}" >&2

        false
    fi

    if [[ ! -s "${partial}" ]]
    then
        printf 'ERROR: empty download: %s\n' "${url}" >&2
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

if [[ ! -s "${fasta}" ]]
then
    gzip -dc "${fasta_gz}" \
        > "${fasta}.tmp.${SLURM_JOB_ID}"

    mv -- \
        "${fasta}.tmp.${SLURM_JOB_ID}" \
        "${fasta}"
fi

if [[ ! -s "${fasta}.fai" ]]
then
    "${samtools_bin}" faidx "${fasta}"
fi

if [[ ! -s "${fasta}.bwt" ]]
then
    "${bwa_bin}" index "${fasta}"
fi

dict_file="${fasta%.fna}.dict"

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

printf '%s\n' "${fasta}" \
    > "${metadata_dir}/reference_fasta.path.tmp.${SLURM_JOB_ID}"

mv -- \
    "${metadata_dir}/reference_fasta.path.tmp.${SLURM_JOB_ID}" \
    "${metadata_dir}/reference_fasta.path"

sha256sum \
    "${fasta_gz}" \
    "${fasta}" \
    "${report}" \
    > "${reference_dir}/reference_files.sha256"

printf 'assembly_accession\tGCA_900068215.1\nassembly_name\tHeliconius_melpomene_rosina.546\nfasta\t%s\nsequences\t%s\nbases\t%s\ncompleted_at\t%s\n' \
    "${fasta}" \
    "$(wc -l < "${fasta}.fai")" \
    "$(awk '{ total += $2 } END { printf "%.0f", total }' "${fasta}.fai")" \
    "$(date --iso-8601=seconds)" \
    > "${complete_file}"

printf 'REFERENCE COMPLETE\n'
cat "${complete_file}"
