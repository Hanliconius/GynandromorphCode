#!/usr/bin/env python3
"""Create analysis 141 from the validated analysis-140 implementation."""
import hashlib
import re
from pathlib import Path

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
SOURCE = PROJECT / "scripts/140_B_western_and_adjusted_LR.py"
TARGET = PROJECT / "scripts/141_B_western_plus_heurippa_adjusted_LR.py"

if not SOURCE.is_file():
    raise RuntimeError(f"Missing validated source: {SOURCE}")

text = SOURCE.read_text()
source_sha = hashlib.sha256(text.encode()).hexdigest()

literal_changes = {
    "B_western_geographic_reanalysis":
        "B_western_plus_heurippa_reanalysis",

    "140_plot_B_western_adjusted_LR.R":
        "141_plot_B_western_plus_heurippa_adjusted_LR.R",

    "Band association in Pacific-west / Central American samples":
        "Band association in western samples plus heurippa",

    "Pacific-west plus Central America":
        "Pacific-west plus Central America plus explicit heurippa sensitivity group",

    'version="140.1"':
        'version="141.1"',

    '"version": "140.1"':
        '"version": "141.1"',
}

for old, new in literal_changes.items():
    if old not in text:
        if old in (
            'version="140.1"',
            '"version": "140.1"',
        ):
            continue

        raise RuntimeError(
            f"Expected source text was not found: {old}"
        )

    text = text.replace(old, new)

pattern = re.compile(
    r'(?ms)^(?P<indent>\s*)'
    r'if key\.startswith\("heurippa"\):\n'
    r'.*?'
    r'(?=^(?P=indent)'
    r'if key\.startswith\("timareta"\):)'
)

matches = list(pattern.finditer(text))

if len(matches) != 1:
    raise RuntimeError(
        "Expected one heurippa geography block; "
        f"found {len(matches)}"
    )

indent = matches[0].group("indent")

replacement = (
    f'{indent}if key.startswith("heurippa"):\n'
    f'{indent}    return (\n'
    f'{indent}        True,\n'
    f'{indent}        "explicit_heurippa_sensitivity_group",\n'
    f'{indent}        "eastern_Andes_sensitivity_group",\n'
    f'{indent}        "Salazar2010",\n'
    f'{indent}    )\n\n'
)

text = pattern.sub(
    replacement,
    text,
    count=1,
)

message = (
    "Pachinus included geographically; "
    "Cauca-valley taxa excluded from this strict definition."
)

if message in text:
    text = text.replace(
        message,
        message
        + "\\nHeurippa explicitly retained "
        + "as a B=1 sensitivity group.",
    )

compile(
    text,
    str(TARGET),
    "exec",
)

TARGET.write_text(text)
target_sha = hashlib.sha256(text.encode()).hexdigest()

print("WESTERN PLUS HEURIPPA SCRIPT PREPARATION")
print(f"Source: {SOURCE}")
print(f"Source SHA256: {source_sha}")
print(f"Created: {TARGET}")
print(f"Created SHA256: {target_sha}")
print(
    "Only cohort/output labels and the explicit "
    "heurippa geography decision changed."
)
