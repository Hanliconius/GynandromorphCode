#!/usr/bin/env python3

from pathlib import Path


project = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

scripts = [
    project
    / "scripts/43_make_aligned_optix_genoplot_score_panels.py",
    project
    / "scripts/45_make_aligned_yb_genoplot_score_panels.py",
]

palette_replacements = [
    (
        '"#F2F2F2"',
        '"#B12A90"',
    ),
    (
        '"#56B4E9"',
        '"#FCA636"',
    ),
    (
        '"#005A8D"',
        '"#2C9CCB"',
    ),
    (
        '"#BDBDBD"',
        '"#E6E6E6"',
    ),
    (
        '"#0072B2"',
        '"#7A0177"',
    ),
    (
        '"Inclusive focal-compatible site"',
        '"Comparative-supported + focal-matching"',
    ),
]

optix_absent_old = '''        (
            absent_n
            / 2
            / n_samples
        ),
'''

optix_absent_new = '''        (
            1
            - absent_n
            / 2
            / n_samples
        ),
'''

optix_present_old = '''        (
            absent_n
            + present_n / 2
        )
        / n_samples,
'''

optix_present_new = '''        (
            1
            - (
                absent_n
                + present_n / 2
            )
            / n_samples
        ),
'''

yb_absent_old = '''        (
            absent_n
            / 2
            / number_of_samples
        ),
'''

yb_absent_new = '''        (
            1
            - absent_n
            / 2
            / number_of_samples
        ),
'''

yb_present_old = '''        (
            absent_n
            + present_n / 2
        )
        / number_of_samples,
'''

yb_present_new = '''        (
            1
            - (
                absent_n
                + present_n / 2
            )
            / number_of_samples
        ),
'''

for script in scripts:
    if not script.is_file():
        print(
            f"ERROR: script is missing: {script}"
        )

        continue

    original = script.read_text()
    updated = original

    backup = script.with_suffix(
        script.suffix + ".before_mago_restyle"
    )

    if not backup.exists():
        backup.write_text(original)

        print(
            f"Backup created: {backup}"
        )

    for old, new in palette_replacements:
        updated = updated.replace(
            old,
            new,
        )

    if script.name.startswith("43_"):
        updated = updated.replace(
            optix_absent_old,
            optix_absent_new,
            1,
        )

        updated = updated.replace(
            optix_present_old,
            optix_present_new,
            1,
        )

    if script.name.startswith("45_"):
        updated = updated.replace(
            yb_absent_old,
            yb_absent_new,
            1,
        )

        updated = updated.replace(
            yb_present_old,
            yb_present_new,
            1,
        )

    if updated == original:
        print(
            f"No changes required: {script}"
        )

    else:
        script.write_text(updated)

        print(
            f"Restyled: {script}"
        )

print()
print("GENOTYPE COLOURS")
print("0 associated copies: #B12A90  magenta")
print("1 associated copy:   #FCA636  orange")
print("2 associated copies: #2C9CCB  blue")
print("Missing/failed QC:   #E6E6E6  pale gray")
print()
print("Eligible QC score points: #7A0177")
print("Current triangle selection was not changed.")
