#!/usr/bin/env python3

from pathlib import Path


source_python = Path(
    "scripts/128_extract_first5000bp_broad_optix_controls.py"
)

target_python = Path(
    "scripts/130_extract_first20000bp_broad_optix_controls.py"
)

source_r = Path(
    "scripts/129_plot_first5000bp_focal_vs_WGS_allele_balance.R"
)

target_r = Path(
    "scripts/131_plot_first20000bp_focal_vs_WGS_allele_balance.R"
)

for required in [source_python, source_r]:
    if not required.is_file():
        raise RuntimeError(f"Missing source script: {required}")


python_text = source_python.read_text()

python_replacements = [
    ("first5000bp", "first20000bp"),
    ("interval_start + 4999", "interval_start + 19999"),
    ('"span_bp": 5000,', '"span_bp": 20000,'),
    (
        "first_5000bp_from_first_broad_VCF_record",
        "first_20000bp_from_first_broad_VCF_record",
    ),
    (
        "FIRST 5000 BP OF BROAD OPTIX VCF",
        "FIRST 20000 BP OF BROAD OPTIX VCF",
    ),
]

for old, new in python_replacements:
    if old not in python_text:
        raise RuntimeError(
            f"Expected Python text was not found: {old}"
        )

    python_text = python_text.replace(old, new)

compile(
    python_text,
    str(target_python),
    "exec",
)

target_python.write_text(python_text)


r_text = source_r.read_text()

r_replacements = [
    ("first5000bp", "first20000bp"),
    ("first 5000-bp", "first 20000-bp"),
]

for old, new in r_replacements:
    if old not in r_text:
        raise RuntimeError(
            f"Expected R text was not found: {old}"
        )

    r_text = r_text.replace(old, new)

target_r.write_text(r_text)

print(f"Created and Python-syntax-checked: {target_python}")
print(f"Created: {target_r}")
