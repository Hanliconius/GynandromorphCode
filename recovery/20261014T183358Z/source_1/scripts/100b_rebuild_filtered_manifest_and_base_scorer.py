#!/usr/bin/env python3

from pathlib import Path
import py_compile
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

OLD_MANIFEST = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)

EXCLUSION_MANIFEST = (
    PROJECT
    / "analysis/drop_problematic_reanalysis/checkpoints"
    / "excluded_comparative_samples.tsv"
)

NEW_MANIFEST = (
    PROJECT
    / "analysis/drop_problematic_reanalysis/checkpoints"
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)

BASE_SOURCE = (
    PROJECT
    / "scripts"
    / "94d_rescore_all_original502_snps_WGS_vs_WGSplusCS.py"
)

BASE_DESTINATION = (
    PROJECT
    / "scripts"
    / "100b_base_scorer_drop_problematic.py"
)


def truthy(value):
    return str(value).strip().lower() in {
        "true", "1", "yes", "t"
    }


def retained_counts(table):
    retained = table[table["keep"].map(truthy)]

    return (
        retained.groupby(
            ["cohort", "module"],
            dropna=False,
        )["sample"]
        .nunique()
        .reset_index(name="retained_samples")
        .sort_values(["cohort", "module"])
    )


def main():
    required = [
        OLD_MANIFEST,
        EXCLUSION_MANIFEST,
        BASE_SOURCE,
    ]

    missing = [
        str(path)
        for path in required
        if not path.is_file()
    ]

    if missing:
        print("ERROR: required inputs are missing:")
        for path in missing:
            print(path)
        return

    NEW_MANIFEST.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    old = pd.read_csv(
        OLD_MANIFEST,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    exclusions = pd.read_csv(
        EXCLUSION_MANIFEST,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    excluded_samples = sorted(
        exclusions["sample"].unique()
    )

    if len(excluded_samples) != 6:
        print(
            "ERROR: expected six excluded samples; found "
            f"{len(excluded_samples)}"
        )
        print(",".join(excluded_samples))
        return

    before = retained_counts(old)

    updated = old.copy()

    target_mask = updated["sample"].isin(
        excluded_samples
    )

    previously_kept = (
        target_mask
        & updated["keep"].map(truthy)
    )

    changed_samples = sorted(
        updated.loc[
            previously_kept,
            "sample",
        ].unique()
    )

    # Preserve every manifest row and all previous filters.
    # Only force the six requested samples to keep=False.
    updated.loc[target_mask, "keep"] = "False"

    updated.to_csv(
        NEW_MANIFEST,
        sep="\t",
        index=False,
    )

    after = retained_counts(updated)

    base_text = BASE_SOURCE.read_text()

    if "region_panel_metadata.tsv" not in base_text:
        print(
            "ERROR: canonical metadata reference was not found "
            "in the base scorer."
        )
        return

    base_text = base_text.replace(
        "region_panel_metadata.tsv",
        "region_panel_metadata.drop_problematic.tsv",
    )

    base_text = base_text.replace(
        '"original_502"',
        '"drop_problematic_reanalysis"',
    )

    base_text = base_text.replace(
        "'original_502'",
        "'drop_problematic_reanalysis'",
    )

    compile(
        base_text,
        str(BASE_DESTINATION),
        "exec",
    )

    BASE_DESTINATION.write_text(base_text)
    BASE_DESTINATION.chmod(0o755)

    print("CORRECTED FILTER-MANIFEST REANALYSIS")
    print()
    print("Requested exclusions:")
    print(",".join(excluded_samples))

    print()
    print("Samples newly changing keep=True to keep=False:")
    if changed_samples:
        print(",".join(changed_samples))
    else:
        print("None")

    print()
    print("BEFORE")
    print(
        before.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print()
    print("AFTER")
    print(
        after.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print()
    print(f"New filtered manifest: {NEW_MANIFEST}")
    print(f"Corrected base scorer: {BASE_DESTINATION}")

    print()
    print("INHERITED TAXON EXCLUSIONS")
    for term in (
        "cydno",
        "numata",
        "pachinus",
        "hecale",
        "ismenius",
        "contigua",
        "linaresi",
    ):
        present = f'"{term}"' in base_text
        print(f"{term}\t{present}")


if __name__ == "__main__":
    main()
