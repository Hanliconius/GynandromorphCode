#!/usr/bin/env bash
#SBATCH --job-name=gyn165map
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=logs/165_collect_resolved_downloads.%j.out
#SBATCH --error=logs/165_collect_resolved_downloads.%j.err

set -euo pipefail

project=/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling
root="${project}/analysis/genomewide_mixing_input"
mapping_dir="${root}/metadata/resolved_runs"
output="${root}/metadata/resolved_download_manifest.tsv"

working_python=/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/pixy/bin/python

cd "${project}"

"${working_python}" - \
    "${mapping_dir}" \
    "${output}" <<'PYEOF'
import sys
from pathlib import Path

import pandas as pd


mapping_dir = Path(sys.argv[1])
output = Path(sys.argv[2])

files = sorted(
    mapping_dir.glob("*.tsv")
)

if len(files) != 22:
    raise RuntimeError(
        "Expected 22 resolved mappings; "
        f"found {len(files)}"
    )

frames = [
    pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    for path in files
]

data = pd.concat(
    frames,
    ignore_index=True,
)

if len(data) != 22:
    raise RuntimeError(
        f"Expected 22 mapping rows; found {len(data)}"
    )

for column in (
    "array_id",
    "source_accession",
    "resolved_run",
    "library_id",
):
    if data[column].duplicated().any():
        raise RuntimeError(
            f"Duplicate values in {column}"
        )

data["_array_id"] = pd.to_numeric(
    data["array_id"],
    errors="raise",
)

data = (
    data.sort_values("_array_id")
    .drop(columns="_array_id")
)

if data["array_id"].tolist() != [
    str(value)
    for value in range(22)
]:
    raise RuntimeError(
        "Resolved array IDs are not exactly 0 through 21"
    )

expected_focal = {
    "TN_Left_T2": (
        "Left",
        "T2",
        "TN_Left_merged",
    ),
    "TN_Right_T2": (
        "Right",
        "T2",
        "TN_Right_merged",
    ),
    "TN_Left_T3": (
        "Left",
        "T3",
        "TN_Left_merged",
    ),
    "TN_Right_T3": (
        "Right",
        "T3",
        "TN_Right_merged",
    ),
}

focal = data.loc[
    data["role"].eq("focal")
]

if set(focal["library_id"]) != set(expected_focal):
    raise RuntimeError(
        "Resolved focal library set is incorrect"
    )

for row in focal.itertuples(index=False):
    expected = expected_focal[
        row.library_id
    ]

    observed = (
        row.side,
        row.leg,
        row.merge_group,
    )

    if observed != expected:
        raise RuntimeError(
            "Incorrect focal mapping for "
            f"{row.library_id}: {observed}"
        )

for row in data.itertuples(index=False):
    fastq1 = Path(row.fastq1)

    if not fastq1.is_file() or fastq1.stat().st_size == 0:
        raise RuntimeError(
            f"Missing FASTQ1: {fastq1}"
        )

    if row.layout == "PAIRED":
        fastq2 = Path(row.fastq2)

        if (
            not fastq2.is_file()
            or fastq2.stat().st_size == 0
        ):
            raise RuntimeError(
                f"Missing FASTQ2: {fastq2}"
            )

    elif row.layout != "SINGLE":
        raise RuntimeError(
            f"Unsupported layout: {row.layout}"
        )

data.to_csv(
    output,
    sep="\t",
    index=False,
)

print(
    "RESOLVED DOWNLOAD MANIFEST COMPLETE"
)

print(
    data.to_string(
        index=False
    )
)
PYEOF

sha256sum "${output}" \
    > "${output}.sha256"

date --iso-8601=seconds \
    > "${root}/metadata/downloads.complete"

printf 'output=%s\n' "${output}"

printf 'rows=%s\n' \
    "$(awk '
        NR > 1 {
            count++
        }

        END {
            print count + 0
        }
    ' "${output}")"
