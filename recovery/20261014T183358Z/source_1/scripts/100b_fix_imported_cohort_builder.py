#!/usr/bin/env python3

from pathlib import Path


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
SCRIPTS = PROJECT / "scripts"

COHORT_SOURCE = (
    SCRIPTS
    / "94b_compare_frozen_snp_scores_WGS_vs_all.py"
)

COHORT_DESTINATION = (
    SCRIPTS
    / "100b_cohort_builder_drop_problematic.py"
)

BASE_SOURCE = (
    SCRIPTS
    / "94d_rescore_all_original502_snps_WGS_vs_WGSplusCS.py"
)

BASE_DESTINATION = (
    SCRIPTS
    / "100b_base_scorer_drop_problematic.py"
)

OLD_COHORT_SCRIPT = (
    "94b_compare_frozen_snp_scores_WGS_vs_all.py"
)

NEW_COHORT_SCRIPT = (
    "100b_cohort_builder_drop_problematic.py"
)


def main():
    required = [
        COHORT_SOURCE,
        BASE_SOURCE,
    ]

    missing = [
        str(path)
        for path in required
        if not path.is_file()
    ]

    if missing:
        print("ERROR: required source scripts are missing:")
        for path in missing:
            print(path)
        return

    cohort_text = COHORT_SOURCE.read_text()

    metadata_count = cohort_text.count(
        "region_panel_metadata.tsv"
    )

    if metadata_count < 1:
        print(
            "ERROR: metadata filename was not found "
            "in the 94b cohort builder."
        )
        return

    cohort_text = cohort_text.replace(
        "region_panel_metadata.tsv",
        "region_panel_metadata.drop_problematic.tsv",
    )

    compile(
        cohort_text,
        str(COHORT_DESTINATION),
        "exec",
    )

    COHORT_DESTINATION.write_text(cohort_text)
    COHORT_DESTINATION.chmod(0o755)

    base_text = BASE_SOURCE.read_text()

    cohort_import_count = base_text.count(
        OLD_COHORT_SCRIPT
    )

    if cohort_import_count < 1:
        print(
            "ERROR: the 94b cohort-builder import was not "
            "found in the base scorer."
        )
        return

    base_text = base_text.replace(
        OLD_COHORT_SCRIPT,
        NEW_COHORT_SCRIPT,
    )

    base_text = base_text.replace(
        "region_panel_metadata.tsv",
        "region_panel_metadata.drop_problematic.tsv",
    )

    compile(
        base_text,
        str(BASE_DESTINATION),
        "exec",
    )

    BASE_DESTINATION.write_text(base_text)
    BASE_DESTINATION.chmod(0o755)

    print("CORRECTED IMPORT CHAIN")
    print()
    print(f"Cohort builder: {COHORT_DESTINATION}")
    print(f"Base scorer:    {BASE_DESTINATION}")
    print()
    print("Cohort builder metadata references:")
    for line_number, line in enumerate(
        cohort_text.splitlines(),
        start=1,
    ):
        if "region_panel_metadata" in line:
            print(f"{line_number}\t{line.strip()}")

    print()
    print("Base scorer cohort-builder references:")
    for line_number, line in enumerate(
        base_text.splitlines(),
        start=1,
    ):
        if "cohort_builder" in line:
            print(f"{line_number}\t{line.strip()}")


if __name__ == "__main__":
    main()
