#!/usr/bin/env python3
"""Download and merge paired ENA FASTQs for one BioSample array task."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import subprocess
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen


FIELDS = "run_accession,sample_accession,library_layout,fastq_ftp"


def fetch_runs(biosample: str) -> list[dict[str, str]]:
    query = urlencode({
        "accession": biosample,
        "result": "read_run",
        "fields": FIELDS,
        "format": "tsv",
    })
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?{query}"
    with urlopen(url, timeout=90) as response:
        text = response.read().decode("utf-8")
    rows = [row for row in csv.DictReader(io.StringIO(text), delimiter="\t")
            if row.get("run_accession")]
    if not rows:
        raise RuntimeError(f"ENA returned no read runs for {biosample}")
    if any(row.get("sample_accession") != biosample for row in rows):
        raise RuntimeError(f"ENA sample accession mismatch for {biosample}")
    if any(row.get("library_layout") != "PAIRED" for row in rows):
        raise RuntimeError(f"Expected paired-end runs for {biosample}")
    return rows


def fastq_pair(row: dict[str, str]) -> tuple[str, str]:
    files = [part for part in row.get("fastq_ftp", "").replace(",", ";").split(";") if part]
    read1 = [part for part in files if "_1.fastq.gz" in part or "_1.fq.gz" in part]
    read2 = [part for part in files if "_2.fastq.gz" in part or "_2.fq.gz" in part]
    if len(read1) != 1 or len(read2) != 1:
        raise RuntimeError(f"Expected one paired FASTQ per run; got {files}")
    return read1[0], read2[0]


def as_https(path: str) -> str:
    return path if path.startswith("https://") else f"https://{path.removeprefix('ftp://')}"


def download(url_path: str, target: Path) -> None:
    partial = target.with_name(target.name + ".part")
    subprocess.run([
        "curl", "--fail", "--location", "--retry", "4", "--retry-delay", "5",
        "--output", str(partial), as_https(url_path),
    ], check=True)
    with gzip.open(partial, "rb") as handle:
        handle.read(1)
    partial.replace(target)


def concatenate(parts: list[Path], target: Path) -> None:
    partial = target.with_name(target.name + ".part")
    with partial.open("wb") as out:
        for part in parts:
            with part.open("rb") as src:
                while chunk := src.read(8 * 1024 * 1024):
                    out.write(chunk)
    with gzip.open(partial, "rb") as handle:
        handle.read(1)
    partial.replace(target)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--task-id", required=True, type=int)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    with args.manifest.open(newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    if not 0 <= args.task_id < len(samples):
        raise RuntimeError(f"Array task {args.task_id} is outside the sample manifest")
    sample = samples[args.task_id]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    runs = fetch_runs(sample["biosample"])
    r1_parts: list[Path] = []
    r2_parts: list[Path] = []
    manifest = args.output_dir / f"{sample['sample']}.ENA_runs.tsv"
    with manifest.open("w") as handle:
        handle.write("run_accession\tbiosample\tr1_url\tr2_url\n")
        for row in runs:
            r1_url, r2_url = fastq_pair(row)
            handle.write(f"{row['run_accession']}\t{sample['biosample']}\t{r1_url}\t{r2_url}\n")
            r1 = args.output_dir / f"{sample['sample']}.{row['run_accession']}_1.fastq.gz"
            r2 = args.output_dir / f"{sample['sample']}.{row['run_accession']}_2.fastq.gz"
            for url, target in ((r1_url, r1), (r2_url, r2)):
                if not target.is_file() or target.stat().st_size == 0:
                    download(url, target)
            r1_parts.append(r1)
            r2_parts.append(r2)

    out_r1 = args.output_dir / f"{sample['sample']}_R1.fastq.gz"
    out_r2 = args.output_dir / f"{sample['sample']}_R2.fastq.gz"
    concatenate(r1_parts, out_r1)
    concatenate(r2_parts, out_r2)
    print(f"YB_FASTQ_READY\t{sample['sample']}\t{len(runs)}\t{out_r1}\t{out_r2}")


if __name__ == "__main__":
    main()
