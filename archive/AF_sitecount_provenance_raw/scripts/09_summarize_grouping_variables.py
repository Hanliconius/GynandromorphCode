#!/usr/bin/env python3

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
METADATA = PROJECT / "data/derived/region_panel_metadata.tsv"
OUTPUT = PROJECT / "analysis/checkpoints/grouping_variable_examples.tsv"

FOCAL = {
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
}

GROUP_COLUMNS = (
    "Grouping-witihin-species",
    "Species_Group",
)

TRAITS = (
    "B",
    "R",
    "D",
    "Yb-dorsal",
    "Yb-ventral",
)


def clean(value):
    if value is None:
        return ""

    value = re.sub(r"\s+", " ", str(value)).strip()
    return value


def state(value):
    try:
        number = float(clean(value))
    except ValueError:
        return None

    return int(number) if number in (0.0, 1.0, 2.0) else None


def figure(value):
    return state(value)


def trait_counts(rows, trait):
    counts = Counter(state(row.get(trait, "")) for row in rows)

    return (
        counts.get(0, 0),
        counts.get(1, 0),
        counts.get(None, 0),
    )


def top_taxa(rows, limit=5):
    counts = Counter(
        clean(row.get("Taxon Name")) or "<missing>"
        for row in rows
    )

    return "; ".join(
        f"{taxon} ({count})"
        for taxon, count in counts.most_common(limit)
    )


def main():
    if not METADATA.is_file():
        print(f"ERROR: missing metadata: {METADATA}")
        return 1

    with METADATA.open(newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if clean(row.get("Accession")) not in FOCAL
        ]

    output_rows = []

    for grouping in GROUP_COLUMNS:
        grouped = defaultdict(list)

        for row in rows:
            group = clean(row.get(grouping)) or "<missing>"
            grouped[group].append(row)

        for group, group_rows in sorted(grouped.items()):
            data_types = Counter(
                clean(row.get("Data Type")) or "<missing>"
                for row in group_rows
            )

            optix_figures = Counter(
                figure(row.get("Optix Figure"))
                for row in group_rows
            )

            ivory_figures = Counter(
                figure(row.get("Ivory Figure"))
                for row in group_rows
            )

            trait_summary = {}

            for trait in TRAITS:
                zeros, ones, missing = trait_counts(group_rows, trait)
                trait_summary[trait] = (zeros, ones, missing)

            output_rows.append(
                [
                    grouping,
                    group,
                    len(group_rows),
                    top_taxa(group_rows),
                    data_types.get("CS", 0),
                    data_types.get("WGS", 0),
                    optix_figures.get(0, 0),
                    optix_figures.get(1, 0),
                    optix_figures.get(2, 0),
                    ivory_figures.get(0, 0),
                    ivory_figures.get(1, 0),
                    ivory_figures.get(2, 0),
                    *trait_summary["B"],
                    *trait_summary["R"],
                    *trait_summary["D"],
                    *trait_summary["Yb-dorsal"],
                    *trait_summary["Yb-ventral"],
                ]
            )

    header = [
        "grouping_variable",
        "group",
        "n",
        "largest_taxa",
        "CS_n",
        "WGS_n",
        "optix_fig0_n",
        "optix_fig1_n",
        "optix_fig2_n",
        "ivory_fig0_n",
        "ivory_fig1_n",
        "ivory_fig2_n",
        "B0",
        "B1",
        "B_missing",
        "R0",
        "R1",
        "R_missing",
        "D0",
        "D1",
        "D_missing",
        "Yb_dorsal0",
        "Yb_dorsal1",
        "Yb_dorsal_missing",
        "Yb_ventral0",
        "Yb_ventral1",
        "Yb_ventral_missing",
    ]

    with OUTPUT.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(output_rows)

    print("GROUPING VARIABLE EXAMPLES")
    print(
        "grouping_variable\tgroup\tn\tlargest_taxa\t"
        "CS\tWGS\tB0\tB1\tR0\tR1\tD0\tD1\t"
        "Ybd0\tYbd1\tYbv0\tYbv1"
    )

    for row in output_rows:
        print(
            f"{row[0]}\t{row[1]}\t{row[2]}\t{row[3]}\t"
            f"{row[4]}\t{row[5]}\t"
            f"{row[12]}\t{row[13]}\t"
            f"{row[15]}\t{row[16]}\t"
            f"{row[18]}\t{row[19]}\t"
            f"{row[21]}\t{row[22]}\t"
            f"{row[24]}\t{row[25]}"
        )

    print(f"\nFull grouping table: {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
