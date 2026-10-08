#!/usr/bin/env python3
"""Map exact Hmel2.5 Yb feature coordinates onto HmRos from minimap2 PAF."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


FEATURES = [
    ("Yb_d", 1_403_500, 1_405_500),
    ("Krista_2024_color_region", 1_408_879, 1_429_719),
    ("cortex_proximal_TSS", 1_424_680, 1_424_680),
    ("Yb_v", 1_429_500, 1_430_500),
    ("ivory_delta78k", 1_494_902, 1_573_177),
]

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_paf(path: Path):
    rows = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            tags = {}
            for tag in fields[12:]:
                parts = tag.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]

            if "cg" not in tags:
                continue

            rows.append(
                {
                    "qname": fields[0],
                    "qlen": int(fields[1]),
                    "qstart": int(fields[2]),
                    "qend": int(fields[3]),
                    "strand": fields[4],
                    "tname": fields[5],
                    "tlen": int(fields[6]),
                    "tstart": int(fields[7]),
                    "tend": int(fields[8]),
                    "matches": int(fields[9]),
                    "alnlen": int(fields[10]),
                    "mapq": int(fields[11]),
                    "cigar": tags["cg"],
                }
            )

    if not rows:
        raise RuntimeError(f"No PAF alignments with cg tag in {path}")

    rows.sort(
        key=lambda x: (x["mapq"], x["matches"], x["alnlen"]),
        reverse=True,
    )
    return rows


def build_map(
    aln,
    old_region_start: int,
    new_region_start: int,
):
    q_to_t = {}
    target_insertions = []
    query_insertions = []

    t = aln["tstart"]

    if aln["strand"] == "+":
        q = aln["qstart"]
        q_step = 1
    else:
        q = aln["qend"] - 1
        q_step = -1

    for length_text, op in CIGAR_RE.findall(aln["cigar"]):
        n = int(length_text)

        if op in ("M", "=", "X"):
            for _ in range(n):
                old_pos = old_region_start + q
                new_pos = new_region_start + t
                q_to_t[old_pos] = new_pos
                q += q_step
                t += 1

        elif op == "I":
            q_first = q
            q_last = q + q_step * (n - 1)
            query_insertions.append(
                {
                    "length": n,
                    "old_start": old_region_start + min(q_first, q_last),
                    "old_end": old_region_start + max(q_first, q_last),
                    "new_anchor": new_region_start + t,
                }
            )
            q += q_step * n

        elif op in ("D", "N"):
            target_insertions.append(
                {
                    "length": n,
                    "old_anchor": old_region_start + q,
                    "new_start": new_region_start + t,
                    "new_end": new_region_start + t + n - 1,
                }
            )
            t += n

        elif op == "S":
            q += q_step * n

        elif op in ("H", "P"):
            continue

        else:
            raise RuntimeError(f"Unhandled CIGAR op {op}")

    return q_to_t, target_insertions, query_insertions


def nearest_map(mapping, position, max_distance=1000):
    if position in mapping:
        return mapping[position], 0

    for distance in range(1, max_distance + 1):
        for candidate in (position - distance, position + distance):
            if candidate in mapping:
                return mapping[candidate], candidate - position

    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True, type=Path)
    ap.add_argument("--old-region-start", required=True, type=int)
    ap.add_argument("--new-region-start", required=True, type=int)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    alignments = parse_paf(args.paf)
    best = alignments[0]

    q_to_t, target_insertions, query_insertions = build_map(
        best,
        args.old_region_start,
        args.new_region_start,
    )

    feature_rows = []
    for name, start, end in FEATURES:
        mapped_start, delta_start = nearest_map(q_to_t, start)
        mapped_end, delta_end = nearest_map(q_to_t, end)

        if mapped_start is not None and mapped_end is not None:
            hmros_start = min(mapped_start, mapped_end)
            hmros_end = max(mapped_start, mapped_end)
            orientation = "+" if mapped_end >= mapped_start else "-"
        else:
            hmros_start = None
            hmros_end = None
            orientation = best["strand"]

        feature_rows.append(
            {
                "feature": name,
                "Hmel2_5_start": start,
                "Hmel2_5_end": end,
                "HmRos_seqid": "NC_140403.1",
                "HmRos_start": hmros_start,
                "HmRos_end": hmros_end,
                "orientation": orientation,
                "start_shift_bp": (
                    mapped_start - start
                    if mapped_start is not None
                    else None
                ),
                "end_shift_bp": (
                    mapped_end - end
                    if mapped_end is not None
                    else None
                ),
                "start_nearest_offset_bp": delta_start,
                "end_nearest_offset_bp": delta_end,
                "alignment_mapq": best["mapq"],
            }
        )

    feature_table = pd.DataFrame(feature_rows)
    feature_path = args.output_dir / "Hmel2.5_to_HmRos_Yb_features.tsv"
    feature_table.to_csv(feature_path, sep="\t", index=False)

    gaps = []
    for row in target_insertions:
        if row["length"] >= 20:
            gaps.append(
                {
                    "type": "sequence_present_in_HmRos_not_Hmel2.5",
                    **row,
                }
            )
    for row in query_insertions:
        if row["length"] >= 20:
            gaps.append(
                {
                    "type": "sequence_present_in_Hmel2.5_not_HmRos",
                    **row,
                }
            )

    gap_table = pd.DataFrame(gaps)
    gap_path = args.output_dir / "Hmel2.5_to_HmRos_structural_gaps_ge20bp.tsv"
    gap_table.to_csv(gap_path, sep="\t", index=False)

    sample_rows = []
    old_min = min(q_to_t)
    old_max = max(q_to_t)
    for old_pos in range(
        ((old_min + 999) // 1000) * 1000,
        old_max + 1,
        1000,
    ):
        if old_pos not in q_to_t:
            continue
        new_pos = q_to_t[old_pos]
        sample_rows.append(
            {
                "Hmel2_5_position": old_pos,
                "HmRos_position": new_pos,
                "shift_bp": new_pos - old_pos,
            }
        )

    shift_table = pd.DataFrame(sample_rows)
    shift_path = args.output_dir / "Hmel2.5_to_HmRos_shift_by_1kb.tsv"
    shift_table.to_csv(shift_path, sep="\t", index=False)

    aln_table = pd.DataFrame(alignments)
    aln_path = args.output_dir / "Hmel2.5_to_HmRos_all_alignments.tsv"
    aln_table.to_csv(aln_path, sep="\t", index=False)

    print("YB_EXACT_COORDINATE_TRANSFER_COMPLETE")
    print("BEST_ALIGNMENT")
    print(pd.DataFrame([best]).to_string(index=False))
    print()
    print("FEATURE_MAPPING")
    print(feature_table.to_string(index=False))
    print()
    print("STRUCTURAL_GAPS_GE20BP")
    if gap_table.empty:
        print("none")
    else:
        print(gap_table.to_string(index=False))
    print()
    print(f"FEATURES\t{feature_path.resolve()}")
    print(f"GAPS\t{gap_path.resolve()}")
    print(f"SHIFTS\t{shift_path.resolve()}")
    print(f"ALIGNMENTS\t{aln_path.resolve()}")


if __name__ == "__main__":
    main()
