#!/usr/bin/env python3
"""Map Hmel2.5 Yb coordinates onto HmRos from a fragmented minimap2 PAF.

Important: this locus is repeat-rich and minimap2 returns many collinear
alignment blocks.  Coordinate transfer must therefore combine all compatible
blocks rather than use only the single longest alignment.
"""

from __future__ import annotations

import argparse
import bisect
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

    return rows


def alignment_map_and_indels(
    aln,
    old_region_start: int,
    new_region_start: int,
):
    mapping = {}
    events = []

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
                mapping[old_region_start + q] = new_region_start + t
                q += q_step
                t += 1

        elif op == "I":
            q_first = q
            q_last = q + q_step * (n - 1)
            if n >= 20:
                events.append(
                    {
                        "event_basis": "within_alignment",
                        "type": "sequence_present_in_Hmel2.5_not_HmRos",
                        "length_bp": n,
                        "Hmel2_5_start": old_region_start
                        + min(q_first, q_last),
                        "Hmel2_5_end": old_region_start
                        + max(q_first, q_last),
                        "HmRos_start": None,
                        "HmRos_end": None,
                        "HmRos_anchor": new_region_start + t,
                        "alignment_qstart": aln["qstart"],
                        "alignment_qend": aln["qend"],
                        "alignment_mapq": aln["mapq"],
                    }
                )
            q += q_step * n

        elif op in ("D", "N"):
            if n >= 20:
                events.append(
                    {
                        "event_basis": "within_alignment",
                        "type": "sequence_present_in_HmRos_not_Hmel2.5",
                        "length_bp": n,
                        "Hmel2_5_start": None,
                        "Hmel2_5_end": None,
                        "Hmel2_5_anchor": old_region_start + q,
                        "HmRos_start": new_region_start + t,
                        "HmRos_end": new_region_start + t + n - 1,
                        "alignment_qstart": aln["qstart"],
                        "alignment_qend": aln["qend"],
                        "alignment_mapq": aln["mapq"],
                    }
                )
            t += n

        elif op == "S":
            q += q_step * n

        elif op in ("H", "P"):
            continue

        else:
            raise RuntimeError(f"Unhandled CIGAR op {op}")

    return mapping, events


def merge_alignment_maps(
    alignments,
    old_region_start: int,
    new_region_start: int,
):
    # Higher MAPQ / more matching sequence gets first claim on any query base
    # represented by overlapping alignment blocks.
    ranked = sorted(
        alignments,
        key=lambda x: (x["mapq"], x["matches"], x["alnlen"]),
        reverse=True,
    )

    merged = {}
    events = []

    for aln in ranked:
        local_map, local_events = alignment_map_and_indels(
            aln,
            old_region_start,
            new_region_start,
        )

        for old_pos, new_pos in local_map.items():
            merged.setdefault(old_pos, new_pos)

        events.extend(local_events)

    return merged, events


def endpoint_mapping(mapping, sorted_old, position):
    if position in mapping:
        new_pos = mapping[position]
        return {
            "status": "exact",
            "mapped": new_pos,
            "shift": new_pos - position,
            "lower_old": position,
            "lower_new": new_pos,
            "upper_old": position,
            "upper_new": new_pos,
        }

    i = bisect.bisect_left(sorted_old, position)
    lower_old = sorted_old[i - 1] if i > 0 else None
    upper_old = sorted_old[i] if i < len(sorted_old) else None

    return {
        "status": "bracketed",
        "mapped": None,
        "shift": None,
        "lower_old": lower_old,
        "lower_new": mapping[lower_old] if lower_old is not None else None,
        "upper_old": upper_old,
        "upper_new": mapping[upper_old] if upper_old is not None else None,
    }


def collinear_block_gaps(
    alignments,
    old_region_start: int,
    new_region_start: int,
):
    # The current Yb comparison is collinear and plus-strand.  Keep the
    # dominant strand and report sequence between adjacent non-overlapping
    # blocks as unresolved inter-block gaps rather than pretending they are
    # exact insertion/deletion calls.
    strand_score = {}
    for aln in alignments:
        strand_score[aln["strand"]] = (
            strand_score.get(aln["strand"], 0) + aln["matches"]
        )

    dominant = max(strand_score, key=strand_score.get)
    blocks = sorted(
        [x for x in alignments if x["strand"] == dominant],
        key=lambda x: (x["qstart"], x["qend"]),
    )

    # Resolve any query-overlapping blocks by retaining the one with the
    # strongest evidence before constructing the collinear chain.
    chain = []
    for aln in blocks:
        if not chain or aln["qstart"] >= chain[-1]["qend"]:
            chain.append(aln)
            continue

        prev = chain[-1]
        prev_score = (prev["mapq"], prev["matches"], prev["alnlen"])
        new_score = (aln["mapq"], aln["matches"], aln["alnlen"])
        if new_score > prev_score:
            chain[-1] = aln

    rows = []
    for left, right in zip(chain, chain[1:]):
        q_gap = right["qstart"] - left["qend"]

        if dominant == "+":
            t_gap = right["tstart"] - left["tend"]
        else:
            # Not expected for this locus, but preserve a usable diagnostic.
            t_gap = left["tstart"] - right["tend"]

        if q_gap < 0 or t_gap < 0:
            continue

        rows.append(
            {
                "strand": dominant,
                "left_Hmel2_5_last_aligned": (
                    old_region_start + left["qend"] - 1
                ),
                "right_Hmel2_5_first_aligned": (
                    old_region_start + right["qstart"]
                ),
                "Hmel2_5_unaligned_gap_bp": q_gap,
                "left_HmRos_last_aligned": (
                    new_region_start + left["tend"] - 1
                ),
                "right_HmRos_first_aligned": (
                    new_region_start + right["tstart"]
                ),
                "HmRos_unaligned_gap_bp": t_gap,
                "HmRos_minus_Hmel2_5_gap_bp": t_gap - q_gap,
                "left_mapq": left["mapq"],
                "right_mapq": right["mapq"],
            }
        )

    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paf", required=True, type=Path)
    ap.add_argument("--old-region-start", required=True, type=int)
    ap.add_argument("--new-region-start", required=True, type=int)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    alignments = parse_paf(args.paf)
    q_to_t, within_events = merge_alignment_maps(
        alignments,
        args.old_region_start,
        args.new_region_start,
    )

    if not q_to_t:
        raise RuntimeError("PAF contained no aligned query bases")

    sorted_old = sorted(q_to_t)

    feature_rows = []
    for name, start, end in FEATURES:
        s = endpoint_mapping(q_to_t, sorted_old, start)
        e = endpoint_mapping(q_to_t, sorted_old, end)

        exact_start = s["mapped"]
        exact_end = e["mapped"]

        if exact_start is not None and exact_end is not None:
            orientation = "+" if exact_end >= exact_start else "-"
        else:
            orientation = None

        feature_rows.append(
            {
                "feature": name,
                "Hmel2_5_start": start,
                "Hmel2_5_end": end,
                "HmRos_seqid": "NC_140403.1",
                "HmRos_start_exact": exact_start,
                "HmRos_end_exact": exact_end,
                "orientation_if_both_exact": orientation,
                "start_status": s["status"],
                "end_status": e["status"],
                "start_shift_bp": s["shift"],
                "end_shift_bp": e["shift"],
                "start_lower_Hmel2_5": s["lower_old"],
                "start_lower_HmRos": s["lower_new"],
                "start_upper_Hmel2_5": s["upper_old"],
                "start_upper_HmRos": s["upper_new"],
                "end_lower_Hmel2_5": e["lower_old"],
                "end_lower_HmRos": e["lower_new"],
                "end_upper_Hmel2_5": e["upper_old"],
                "end_upper_HmRos": e["upper_new"],
            }
        )

    feature_table = pd.DataFrame(feature_rows)
    feature_path = args.output_dir / "Hmel2.5_to_HmRos_Yb_features.tsv"
    feature_table.to_csv(feature_path, sep="\t", index=False)

    event_table = pd.DataFrame(within_events)
    event_path = (
        args.output_dir
        / "Hmel2.5_to_HmRos_structural_gaps_ge20bp.tsv"
    )
    event_table.to_csv(event_path, sep="\t", index=False)

    block_gap_table = pd.DataFrame(
        collinear_block_gaps(
            alignments,
            args.old_region_start,
            args.new_region_start,
        )
    )
    block_gap_path = (
        args.output_dir
        / "Hmel2.5_to_HmRos_interblock_gaps.tsv"
    )
    block_gap_table.to_csv(block_gap_path, sep="\t", index=False)

    sample_rows = []
    old_min = min(sorted_old)
    old_max = max(sorted_old)

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

    aln_table = pd.DataFrame(
        sorted(
            alignments,
            key=lambda x: (x["qstart"], x["qend"]),
        )
    )
    aln_path = args.output_dir / "Hmel2.5_to_HmRos_all_alignments.tsv"
    aln_table.to_csv(aln_path, sep="\t", index=False)

    print("YB_EXACT_COORDINATE_TRANSFER_COMPLETE")
    print(
        "alignment_blocks="
        f"{len(alignments)} merged_aligned_bases={len(q_to_t)}"
    )
    print()

    print("FEATURE_MAPPING")
    print(feature_table.to_string(index=False))
    print()

    print("WITHIN_ALIGNMENT_STRUCTURAL_GAPS_GE20BP")
    if event_table.empty:
        print("none")
    else:
        print(event_table.to_string(index=False))
    print()

    print("INTERBLOCK_GAPS")
    if block_gap_table.empty:
        print("none")
    else:
        print(block_gap_table.to_string(index=False))
    print()

    print(f"FEATURES\t{feature_path.resolve()}")
    print(f"GAPS\t{event_path.resolve()}")
    print(f"INTERBLOCK_GAPS\t{block_gap_path.resolve()}")
    print(f"SHIFTS\t{shift_path.resolve()}")
    print(f"ALIGNMENTS\t{aln_path.resolve()}")


if __name__ == "__main__":
    main()
