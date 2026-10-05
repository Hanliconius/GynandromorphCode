#!/usr/bin/env bash
#SBATCH --job-name=gyn166ref
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/166_prepare_rosina_reference.%j.out
#SBATCH --error=logs/166_prepare_rosina_reference.%j.err

set -euo pipefail

tool_bin=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin
export PATH="${tool_bin}:/usr/bin:/bin:${PATH}"

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
root="${project}/analysis/genomewide_mixing_input"
reference_dir="${root}/reference/GCA_900068215.1"
metadata_dir="${root}/metadata"

assembly_base=GCA_900068215.1_Heliconius_melpomene_rosina.546

base_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/068/215/${assembly_base}"

fasta_gz="${reference_dir}/${assembly_base}_genomic.fna.gz"
fasta="${reference_dir}/${assembly_base}_genomic.fna"
report="${reference_dir}/${assembly_base}_assembly_report.txt"
complete="${reference_dir}/reference.complete"

cd "${project}"

mkdir -p \
    "${reference_dir}" \
    "${metadata_dir}"

for tool in curl gzip bwa samtools sha256sum
do
    command -v "${tool}" >/dev/null 2>&1 || {
        printf 'Required command is unavailable: %s\n' \
            "${tool}" >&2
        false
    }
done

already_complete=0

if [[ \
    -s "${complete}" \
    && -s "${fasta}" \
    && -s "${fasta}.fai" \
    && -s "${fasta}.bwt" \
]]
then
    already_complete=1
fi

if [[ "${already_complete}" -eq 1 ]]
then
    printf 'REFERENCE ALREADY COMPLETE\n'
else
    curl -fL \
        --retry 5 \
        "${base_url}/${assembly_base}_genomic.fna.gz" \
        -o "${fasta_gz}.part"

    mv \
        "${fasta_gz}.part" \
        "${fasta_gz}"

    curl -fL \
        --retry 5 \
        "${base_url}/${assembly_base}_assembly_report.txt" \
        -o "${report}.part"

    mv \
        "${report}.part" \
        "${report}"

    gzip -t "${fasta_gz}"

    gzip -dc "${fasta_gz}" \
        > "${fasta}.part"

    mv \
        "${fasta}.part" \
        "${fasta}"

    samtools faidx "${fasta}"

    bwa index "${fasta}"

    samtools dict \
        -o "${fasta%.fna}.dict" \
        "${fasta}"

    printf '%s\n' "${fasta}" \
        > "${metadata_dir}/reference_fasta.path"

    sha256sum \
        "${fasta}" \
        "${fasta}.fai" \
        "${report}" \
        > "${metadata_dir}/reference_checksums.sha256"

    date --iso-8601=seconds \
        > "${complete}"
fi

printf 'REFERENCE READY\n'
printf 'fasta=%s\n' "${fasta}"

printf 'sequences=%s\n' \
    "$(wc -l < "${fasta}.fai")"

printf 'bases=%s\n' \
    "$(awk '
        {
            total += $2
        }

        END {
            printf "%.0f", total
        }
    ' "${fasta}.fai")"

cat "${metadata_dir}/reference_checksums.sha256"
