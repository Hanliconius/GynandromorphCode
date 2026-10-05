#!/usr/bin/env python3

import csv
from collections import Counter, defaultdict
from pathlib import Path


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
METADATA = PROJECT / "data/derived/region_panel_metadata.tsv"
CHECKPOINTS = PROJECT / "analysis/checkpoints"

GROUP_DETAILS = CHECKPOINTS / "phenotype_group_structure.tsv"
FIGURE2_DETAILS = CHECKPOINTS / "nonfocal_figure2_samples.tsv"

FOCAL = {
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
}

OPTIX_TRAITS = ("B", "R", "D")
IVORY_TRAITS = ("Yb-dorsal", "Yb-ventral")

GROUP_COLUMNS = (
    "Taxon Name",
    "Grouping-witihin-species",
    "Species_Group",
)


def clean(value):
    if value is None:
        return ""
    return str(value).strip()


def numeric_state(value):
    value = clean(value)

    if not value:
        return None

    try:
        number = float(value)
    except ValueError:
        return None

    if number in (0.0, 1.0, 2.0):
        return int(number)

    return None


def binary_state(row, trait):
    state = numeric_state(row.get(trait, ""))
    return state if state in (0, 1) else None


def figure_state(row, column):
    return numeric_state(row.get(column, ""))


def complete(row, traits):
    return all(binary_state(row, trait) is not None for trait in traits)


def label(value):
    return clean(value) if clean(value) else "<missing>"


def print_counter(counter, header):
    print(header)
    for key, count in sorted(
        counter.items(),
        key=lambda item: (str(item[0]), str(item[1])),
    ):
        if not isinstance(key, tuple):
            key = (key,)
        print("\t".join(str(value) for value in (*key, count)))


def main():
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)

    if not METADATA.is_file():
        print(f"ERROR: missing metadata table: {METADATA}")
        return 1

    with METADATA.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    nonfocal = [
        row
        for row in rows
        if clean(row.get("Accession")) not in FOCAL
    ]

    focal = [
        row
        for row in rows
        if clean(row.get("Accession")) in FOCAL
    ]

    print("ANALYSIS DESIGN AUDIT")
    print(f"All metadata rows: {len(rows)}")
    print(f"Non-focal rows: {len(nonfocal)}")
    print(f"Focal rows withheld from scoring: {len(focal)}")

    optix_strata = Counter()
    ivory_strata = Counter()

    for row in nonfocal:
        data_type = label(row.get("Data Type"))
        optix_figure = figure_state(row, "Optix Figure")
        ivory_figure = figure_state(row, "Ivory Figure")

        optix_status = (
            "complete_BRD"
            if complete(row, OPTIX_TRAITS)
            else "incomplete_BRD"
        )

        ivory_status = (
            "complete_Yb_pair"
            if complete(row, IVORY_TRAITS)
            else "incomplete_Yb_pair"
        )

        optix_strata[
            (
                optix_figure if optix_figure is not None else "missing",
                data_type,
                optix_status,
            )
        ] += 1

        ivory_strata[
            (
                ivory_figure if ivory_figure is not None else "missing",
                data_type,
                ivory_status,
            )
        ] += 1

    print("\nOPTIX COHORT STRATA")
    print("figure\tdata_type\tphenotype_status\tn")
    for key, count in sorted(
        optix_strata.items(),
        key=lambda item: tuple(str(value) for value in item[0]),
    ):
        print("\t".join(str(value) for value in (*key, count)))

    print("\nIVORY COHORT STRATA")
    print("figure\tdata_type\tphenotype_status\tn")
    for key, count in sorted(
        ivory_strata.items(),
        key=lambda item: tuple(str(value) for value in item[0]),
    ):
        print("\t".join(str(value) for value in (*key, count)))

    optix_scopes = {
        "all_nonfocal_complete": [
            row
            for row in nonfocal
            if complete(row, OPTIX_TRAITS)
        ],
        "optix_figure1_complete": [
            row
            for row in nonfocal
            if figure_state(row, "Optix Figure") == 1
            and complete(row, OPTIX_TRAITS)
        ],
    }

    print("\nJOINT B/R/D PHENOTYPE COMBINATIONS")
    print("scope\tB\tR\tD\tn")

    for scope, scope_rows in optix_scopes.items():
        combinations = Counter(
            (
                binary_state(row, "B"),
                binary_state(row, "R"),
                binary_state(row, "D"),
            )
            for row in scope_rows
        )

        for combination, count in sorted(combinations.items()):
            print(
                f"{scope}\t{combination[0]}\t"
                f"{combination[1]}\t{combination[2]}\t{count}"
            )

    print("\nOPTIX TRAIT COMPLEMENTARITY")
    print(
        "scope\tn\tB_not_equal_1_minus_R\t"
        "B_not_equal_1_minus_D\tR_not_equal_D"
    )

    for scope, scope_rows in optix_scopes.items():
        b_not_inverse_r = sum(
            binary_state(row, "B")
            != 1 - binary_state(row, "R")
            for row in scope_rows
        )

        b_not_inverse_d = sum(
            binary_state(row, "B")
            != 1 - binary_state(row, "D")
            for row in scope_rows
        )

        r_not_d = sum(
            binary_state(row, "R")
            != binary_state(row, "D")
            for row in scope_rows
        )

        print(
            f"{scope}\t{len(scope_rows)}\t"
            f"{b_not_inverse_r}\t"
            f"{b_not_inverse_d}\t"
            f"{r_not_d}"
        )

    ivory_scopes = {
        "all_nonfocal_complete": [
            row
            for row in nonfocal
            if complete(row, IVORY_TRAITS)
        ],
        "ivory_figure1_complete": [
            row
            for row in nonfocal
            if figure_state(row, "Ivory Figure") == 1
            and complete(row, IVORY_TRAITS)
        ],
    }

    print("\nJOINT Yb PHENOTYPE COMBINATIONS")
    print("scope\tYb_dorsal\tYb_ventral\tn")

    for scope, scope_rows in ivory_scopes.items():
        combinations = Counter(
            (
                binary_state(row, "Yb-dorsal"),
                binary_state(row, "Yb-ventral"),
            )
            for row in scope_rows
        )

        for combination, count in sorted(combinations.items()):
            print(
                f"{scope}\t{combination[0]}\t"
                f"{combination[1]}\t{count}"
            )

    print("\nYb DORSAL/VENTRAL DISCORDANCE")
    print("scope\tn\tdorsal_not_equal_ventral")

    for scope, scope_rows in ivory_scopes.items():
        discordant = sum(
            binary_state(row, "Yb-dorsal")
            != binary_state(row, "Yb-ventral")
            for row in scope_rows
        )

        print(f"{scope}\t{len(scope_rows)}\t{discordant}")

    figure2_rows = [
        row
        for row in nonfocal
        if (
            figure_state(row, "Optix Figure") == 2
            or figure_state(row, "Ivory Figure") == 2
        )
    ]

    figure2_fields = [
        "Accession",
        "Taxon Name",
        "Country",
        "Data Type",
        "B",
        "R",
        "D",
        "Yb-dorsal",
        "Yb-ventral",
        "Optix Figure",
        "Ivory Figure",
        "note",
    ]

    with FIGURE2_DETAILS.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=figure2_fields,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(figure2_rows)

    print("\nNON-FOCAL FIGURE-2 SAMPLES")
    print(
        "accession\ttaxon\tdata_type\tB\tR\tD\t"
        "Yb_dorsal\tYb_ventral\toptix_figure\tivory_figure"
    )

    for row in figure2_rows:
        print(
            f"{clean(row.get('Accession'))}\t"
            f"{label(row.get('Taxon Name'))}\t"
            f"{label(row.get('Data Type'))}\t"
            f"{label(row.get('B'))}\t"
            f"{label(row.get('R'))}\t"
            f"{label(row.get('D'))}\t"
            f"{label(row.get('Yb-dorsal'))}\t"
            f"{label(row.get('Yb-ventral'))}\t"
            f"{label(row.get('Optix Figure'))}\t"
            f"{label(row.get('Ivory Figure'))}"
        )

    group_detail_rows = []

    print("\nGROUP-STRUCTURE SUMMARY")
    print(
        "cohort\tgrouping\tgroups\tsingleton_groups\t"
        "trait\tgroups_with_0\tgroups_with_1\t"
        "groups_with_both_states"
    )

    cohorts = {
        "optix_figure1": optix_scopes["optix_figure1_complete"],
        "ivory_figure1": ivory_scopes["ivory_figure1_complete"],
    }

    cohort_traits = {
        "optix_figure1": OPTIX_TRAITS,
        "ivory_figure1": IVORY_TRAITS,
    }

    for cohort_name, cohort_rows in cohorts.items():
        for group_column in GROUP_COLUMNS:
            grouped = defaultdict(list)

            for row in cohort_rows:
                group = clean(row.get(group_column))

                if group:
                    grouped[group].append(row)

            singleton_groups = sum(
                len(group_rows) == 1
                for group_rows in grouped.values()
            )

            for trait in cohort_traits[cohort_name]:
                groups_with_0 = 0
                groups_with_1 = 0
                groups_with_both = 0

                for group, group_rows in grouped.items():
                    states = {
                        binary_state(row, trait)
                        for row in group_rows
                        if binary_state(row, trait) is not None
                    }

                    groups_with_0 += int(0 in states)
                    groups_with_1 += int(1 in states)
                    groups_with_both += int(states == {0, 1})

                    group_detail_rows.append(
                        [
                            cohort_name,
                            group_column,
                            group,
                            len(group_rows),
                            trait,
                            sum(
                                binary_state(row, trait) == 0
                                for row in group_rows
                            ),
                            sum(
                                binary_state(row, trait) == 1
                                for row in group_rows
                            ),
                        ]
                    )

                print(
                    f"{cohort_name}\t{group_column}\t"
                    f"{len(grouped)}\t{singleton_groups}\t"
                    f"{trait}\t{groups_with_0}\t"
                    f"{groups_with_1}\t{groups_with_both}"
                )

    with GROUP_DETAILS.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "cohort",
                "grouping_variable",
                "group",
                "group_n",
                "trait",
                "trait_0_n",
                "trait_1_n",
            ]
        )
        writer.writerows(group_detail_rows)

    print(f"\nDetailed group table: {GROUP_DETAILS}")
    print(f"Non-focal Figure-2 table: {FIGURE2_DETAILS}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
