#!/usr/bin/env python3
"""Resolve cortex/Yb neighbourhood on the H. m. rosina RefSeq assembly.

The final coverage analysis uses GCF_056956725.1 (ASM5695672v1), not Hmel2.5.
This script uses the NCBI assembly report to identify chromosome 15 and the
RefSeq GFF to identify the annotated cortex gene. It deliberately does not
pretend that historical Hmel2.5 Yb(d)/Yb(v) coordinates are directly portable.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
import re

import pandas as pd


def parse_attributes(text: str) -> dict[str, str]:
    result = {}
    for item in text.rstrip(";").split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        result[key] = value
    return result


def chromosome15_ids(report: Path) -> set[str]:
    ids: set[str] = set()

    with report.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            (
                sequence_name,
                sequence_role,
                assigned_molecule,
                molecule_type,
                genbank,
                relationship,
                refseq,
                assembly_unit,
                length,
                ucsc_name,
            ) = fields[:10]

            if assigned_molecule != "15":
                continue

            for value in (
                sequence_name,
                genbank,
                refseq,
                ucsc_name,
            ):
                if value and value.lower() != "na":
                    ids.add(value)

    if not ids:
        raise RuntimeError(
            f"Could not identify chromosome 15 in assembly report: {report}"
        )

    return ids


def find_cortex(gff: Path, chr15_ids: set[str]) -> pd.DataFrame:
    opener = gzip.open if gff.suffix == ".gz" else open
    records = []

    with opener(gff, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attrs = fields

            if feature != "gene":
                continue

            attr = parse_attributes(attrs)
            searchable = " ".join(
                [
                    attrs,
                    attr.get("Name", ""),
                    attr.get("gene", ""),
                    attr.get("gene_synonym", ""),
                    attr.get("description", ""),
                ]
            ).lower()

            if "cortex" not in searchable:
                continue

            records.append(
                {
                    "seqid": seqid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "source": source,
                    "attributes": attrs,
                    "on_chr15": seqid in chr15_ids,
                }
            )

    data = pd.DataFrame(records)

    if data.empty:
        raise RuntimeError(
            "No GFF gene record containing 'cortex' was found."
        )

    on15 = data.loc[data["on_chr15"]].copy()

    if len(on15) != 1:
        raise RuntimeError(
            "Expected exactly one cortex gene on chromosome 15; found "
            f"{len(on15)}. Candidates:\n"
            + data.to_string(index=False)
        )

    return on15


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly-report", required=True, type=Path)
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--reference-fai", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--flank", type=int, default=125_000)
    args = parser.parse_args()

    chr15 = chromosome15_ids(args.assembly_report)
    cortex = find_cortex(args.gff, chr15).iloc[0]

    fai = pd.read_csv(
        args.reference_fai,
        sep="\t",
        header=None,
        names=["seqid", "length", "offset", "line_bases", "line_width"],
    )

    match = fai.loc[fai["seqid"].astype(str).eq(str(cortex["seqid"]))]

    if len(match) != 1:
        raise RuntimeError(
            f"Cortex seqid {cortex['seqid']} not uniquely present in FASTA index."
        )

    chrom_length = int(match.iloc[0]["length"])
    region_start = max(1, int(cortex["start"]) - args.flank)
    region_end = min(chrom_length, int(cortex["end"]) + args.flank)

    out = pd.DataFrame(
        [
            {
                "assembly": "GCF_056956725.1",
                "assembly_name": "ASM5695672v1",
                "chromosome": "15",
                "seqid": cortex["seqid"],
                "chromosome_length": chrom_length,
                "cortex_start": int(cortex["start"]),
                "cortex_end": int(cortex["end"]),
                "cortex_strand": cortex["strand"],
                "plot_start": region_start,
                "plot_end": region_end,
                "plot_flank": args.flank,
                "historical_Hmel2_5_Krista_region":
                    "Hmel215003o:1408879-1429719",
                "historical_Yb_d":
                    "Hmel215003o:1403500-1405500",
                "historical_Yb_v":
                    "Hmel215003o:1429500-1430500",
                "note":
                    "Historical intervals are recorded for provenance only; "
                    "they are not treated as coordinates on ASM5695672v1.",
            }
        ]
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)

    print("HMROS_YB_COORDINATE_AUDIT")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
