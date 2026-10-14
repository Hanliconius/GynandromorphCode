#!/usr/bin/env python3

from pathlib import Path
import py_compile


SCRIPT = Path(
    "/scratch/martinlab/joe/GynandroCalling/"
    "scripts/101a_identify_best_B_D_R_site_breakers.py"
)

BACKUP = SCRIPT.with_suffix(
    ".py.before_dynamic_cohort_size_fix"
)

old = '''        expected_n = EXPECTED_COHORT_SIZES[module]

        if len(samples) != expected_n:
'''

new = '''        expected_n = int(
            number(top["comparative_cohort_n"])
        )

        if len(samples) != expected_n:
'''

if not SCRIPT.is_file():
    print(f"ERROR: script is missing: {SCRIPT}")
else:
    text = SCRIPT.read_text()

    if old not in text:
        print("ERROR: expected cohort-size block was not found.")
        print("No file was changed.")
    else:
        if not BACKUP.exists():
            BACKUP.write_text(text)
            print(f"Backup created: {BACKUP}")

        SCRIPT.write_text(text.replace(old, new, 1))
        py_compile.compile(str(SCRIPT), doraise=True)

        print(f"Patched and syntax-checked: {SCRIPT}")
