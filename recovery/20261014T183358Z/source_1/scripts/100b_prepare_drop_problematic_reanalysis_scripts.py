#!/usr/bin/env python3

from pathlib import Path
import py_compile


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
SCRIPTS = PROJECT / "scripts"

BASE_SOURCE = (
    SCRIPTS
    / "94d_rescore_all_original502_snps_WGS_vs_WGSplusCS.py"
)

BASE_DESTINATION = (
    SCRIPTS
    / "100b_base_scorer_drop_problematic.py"
)

SCANNER_MAP = {
    "95d_scan_long_band_original502_WGS_vs_WGSplusCS.py":
        "100c_scan_band_drop_problematic.py",

    "96a_scan_broad_optix_Dennis.py":
        "100d_scan_Dennis_drop_problematic.py",

    "96a_scan_broad_optix_Ray.py":
        "100d_scan_Ray_drop_problematic.py",
}

OLD_BASE_NAME = (
    "94d_rescore_all_original502_snps_WGS_vs_WGSplusCS.py"
)

NEW_BASE_NAME = (
    "100b_base_scorer_drop_problematic.py"
)


def main():
    required = [BASE_SOURCE] + [
        SCRIPTS / source
        for source in SCANNER_MAP
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

    # Clone the established base scorer. Only its metadata input changes.
    # Its taxon, callability, focal, piano-key and phenotype filters remain.
    base_text = BASE_SOURCE.read_text()

    metadata_occurrences = base_text.count(
        "region_panel_metadata.tsv"
    )

    if metadata_occurrences < 1:
        print(
            "ERROR: metadata filename was not found "
            "in the base scorer."
        )
        return

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

    created = [BASE_DESTINATION]

    for source_name, destination_name in SCANNER_MAP.items():
        source = SCRIPTS / source_name
        destination = SCRIPTS / destination_name

        text = source.read_text()

        base_occurrences = text.count(OLD_BASE_NAME)

        if base_occurrences < 1:
            print(
                f"ERROR: {source_name} does not reference "
                f"{OLD_BASE_NAME}"
            )
            return

        text = text.replace(
            OLD_BASE_NAME,
            NEW_BASE_NAME,
        )

        text = text.replace(
            '"original_502"',
            '"drop_problematic_reanalysis"',
        )

        text = text.replace(
            "'original_502'",
            "'drop_problematic_reanalysis'",
        )

        if NEW_BASE_NAME not in text:
            print(
                f"ERROR: cloned base scorer was not installed "
                f"in {destination_name}"
            )
            return

        if "drop_problematic_reanalysis" not in text:
            print(
                f"ERROR: isolated output directory was not installed "
                f"in {destination_name}"
            )
            return

        compile(
            text,
            str(destination),
            "exec",
        )

        destination.write_text(text)
        destination.chmod(0o755)
        created.append(destination)

    print("CORRECTED DROP-PROBLEMATIC SCANNERS")
    print()
    print("Inherited taxon exclusions:")
    print("cydno")
    print("numata")
    print("pachinus")
    print()
    print("Metadata:")
    print(
        PROJECT
        / "data/derived"
        / "region_panel_metadata.drop_problematic.tsv"
    )
    print()
    print("Created scripts:")

    for path in created:
        print(path)


if __name__ == "__main__":
    main()
