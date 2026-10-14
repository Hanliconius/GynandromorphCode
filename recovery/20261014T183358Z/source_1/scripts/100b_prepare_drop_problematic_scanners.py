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

    # Clone the established scoring implementation. Only replace the
    # canonical metadata input; retain its taxon, callability, focal,
    # piano-key and phenotype filters.
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

    required_taxon_terms = [
        '"cydno"',
        '"numata"',
        '"pachinus"',
    ]

    absent_terms = [
        term for term in required_taxon_terms
        if term not in base_text
    ]

    if absent_terms:
        print(
            "ERROR: inherited taxon filters were not found: "
            + ",".join(absent_terms)
        )
        return

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

        if OLD_BASE_NAME not in text:
            print(
                "ERROR: imported base-scorer filename was not found "
                f"in {source}"
            )
            return

        text = text.replace(
            OLD_BASE_NAME,
            NEW_BASE_NAME,
        )

        # Isolate new outputs from the original analysis.
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
                f"ERROR: new base scorer is absent from {destination}"
            )
            return

        if "drop_problematic_reanalysis" not in text:
            print(
                f"ERROR: isolated output path is absent from {destination}"
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

    print("DROP-PROBLEMATIC FULL-VCF SCANNERS")
    print()

    for path in created:
        print(path)

    print()
    print("INHERITED TAXON EXCLUSIONS")
    print("cydno")
    print("numata")
    print("pachinus")

    print()
    print("NEW METADATA INPUT")
    print(
        PROJECT
        / "data/derived"
        / "region_panel_metadata.drop_problematic.tsv"
    )

    print()
    print("NEW OUTPUT DIRECTORY")
    print(
        PROJECT
        / "analysis/drop_problematic_reanalysis/checkpoints"
    )


if __name__ == "__main__":
    main()
