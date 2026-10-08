#!/usr/bin/env python3
"""Extract one Hmel2.5 scaffold from a bigWig and compute its true mean depth.

Uncovered bases are counted as zero when calculating the full-scaffold mean,
matching a whole-scaffold depth denominator rather than the mean of only
bigWig-covered bases.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pyBigWig


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bigwig", required=True, type=Path)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--bedgraph", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()

    if not args.bigwig.is_file() or args.bigwig.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty bigWig: {args.bigwig}")

    bw = pyBigWig.open(str(args.bigwig))
    try:
        chroms = bw.chroms()

        if args.chrom not in chroms:
            raise RuntimeError(
                f"{args.chrom} is absent from {args.bigwig}; "
                f"available chromosomes={len(chroms)}"
            )

        chrom_length = int(chroms[args.chrom])
        intervals = bw.intervals(args.chrom, 0, chrom_length)

        if intervals is None:
            intervals = []

        numerator = 0.0

        args.bedgraph.parent.mkdir(parents=True, exist_ok=True)

        with args.bedgraph.open("w") as out:
            for start, end, value in intervals:
                start = int(start)
                end = int(end)
                value = float(value)

                if end <= start:
                    continue

                numerator += (end - start) * value
                out.write(
                    f"{args.chrom}\t{start}\t{end}\t{value:.10g}\n"
                )

        scaffold_mean = numerator / chrom_length

        if scaffold_mean <= 0:
            raise RuntimeError(
                f"Non-positive full-scaffold mean depth: {scaffold_mean}"
            )

        args.stats.parent.mkdir(parents=True, exist_ok=True)
        args.stats.write_text(
            "sample\tscaffold_length\tscaffold_mean_depth\tbedgraph\n"
            f"{args.sample}\t{chrom_length}\t{scaffold_mean:.10f}\t"
            f"{args.bedgraph.resolve()}\n"
        )

        print(
            "BIGWIG_SCAFFOLD_EXTRACTED",
            args.sample,
            args.chrom,
            chrom_length,
            f"{scaffold_mean:.10f}",
            sep="\t",
        )

    finally:
        bw.close()


if __name__ == "__main__":
    main()
