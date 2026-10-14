#!/usr/bin/env python3

from pathlib import Path


source_python = Path(
    "scripts/126_extract_first500bp_broad_optix_controls.py"
)

target_python = Path(
    "scripts/128_extract_first5000bp_broad_optix_controls.py"
)

source_r = Path(
    "scripts/127_plot_first500bp_focal_vs_WGS_allele_balance.R"
)

target_r = Path(
    "scripts/129_plot_first5000bp_focal_vs_WGS_allele_balance.R"
)

for required in [source_python, source_r]:
    if not required.is_file():
        raise RuntimeError(f"Missing source script: {required}")


# Python extraction script
python_text = source_python.read_text()

python_replacements = [
    ("first500bp", "first5000bp"),
    ("interval_start + 499", "interval_start + 4999"),
    ('"span_bp": 500,', '"span_bp": 5000,'),
    (
        "first_500bp_from_first_broad_VCF_record",
        "first_5000bp_from_first_broad_VCF_record",
    ),
    (
        "FIRST 500 BP OF BROAD OPTIX VCF",
        "FIRST 5000 BP OF BROAD OPTIX VCF",
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


# R plotting script
r_text = source_r.read_text()

r_replacements = [
    ("first500bp", "first5000bp"),
    ("first 500-bp", "first 5000-bp"),
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
