#!/usr/bin/env python3

from pathlib import Path
import py_compile


SCRIPT = Path(
    "/scratch/martinlab/joe/GynandroCalling/"
    "scripts/99h_audit_comparative_B_D_R_breakers.py"
)

BACKUP = SCRIPT.with_suffix(
    ".py.before_diagnostic_group_fix"
)

replacements = [
    (
'''    genotypes = pd.read_csv(
        GENOTYPES,
        sep="\\t",
        low_memory=False,
    )

    # Select SNPs entirely from the comparative classification.
''',
'''    genotypes = pd.read_csv(
        GENOTYPES,
        sep="\\t",
        low_memory=False,
    )

    # Convert the display labels used in the genotype table to the
    # internal diagnostic-group labels used by this audit.
    diagnostic_group_map = {
        "H. besckei": "besckei_B",
        "H. m. meriana": "meriana_D",
        "H. t. contigua": "contigua_R",
    }

    genotypes["diagnostic_group"] = (
        genotypes["group"].map(diagnostic_group_map)
    )

    # Select SNPs entirely from the comparative classification.
'''
    ),
    (
'''            & genotypes["group"].isin(EXPECTED_GROUP_SIZES)
''',
'''            & genotypes["diagnostic_group"].isin(
                EXPECTED_GROUP_SIZES
            )
'''
    ),
    (
'''            genotypes[
                genotypes["group"].isin(EXPECTED_GROUP_SIZES)
            ][["sample", "group", "sample_label"]]
            .drop_duplicates("sample")
''',
'''            genotypes[
                genotypes["diagnostic_group"].isin(
                    EXPECTED_GROUP_SIZES
                )
            ][["sample", "diagnostic_group", "sample_label"]]
            .drop_duplicates("sample")
'''
    ),
    (
'''            on=["sample", "group", "sample_label"],
''',
'''            on=["sample", "diagnostic_group", "sample_label"],
'''
    ),
    (
'''            group = str(call["group"])
''',
'''            group = str(call["diagnostic_group"])
'''
    ),
]

if not SCRIPT.is_file():
    print(f"ERROR: script is missing: {SCRIPT}")
else:
    text = SCRIPT.read_text()
    original = text
    failures = []

    for old, new in replacements:
        count = text.count(old)

        if count != 1:
            failures.append(
                f"Expected one occurrence but found {count}: "
                f"{old.splitlines()[0]}"
            )
        else:
            text = text.replace(old, new, 1)

    if failures:
        print("ERROR: patch was not applied.")
        for failure in failures:
            print(failure)
    else:
        if not BACKUP.exists():
            BACKUP.write_text(original)
            print(f"Backup created: {BACKUP}")

        SCRIPT.write_text(text)
        py_compile.compile(str(SCRIPT), doraise=True)

        print(f"Patched and syntax-checked: {SCRIPT}")
