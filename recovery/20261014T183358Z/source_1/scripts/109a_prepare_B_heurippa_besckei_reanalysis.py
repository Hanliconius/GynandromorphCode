#!/usr/bin/env python3

from pathlib import Path
import re
import py_compile

import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

SOURCE_MANIFEST = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis"
    / "B_heurippa_besckei_reanalysis"
    / "checkpoints"
)

OUTPUT_MANIFEST = (
    OUTDIR
    / "original502_filtered_cohort_manifest."
    "B_heurippa_besckei.WGS_vs_WGSplusCS.tsv"
)

BASE_SOURCE = (
    PROJECT
    / "scripts"
    / "100b_base_scorer_drop_problematic.py"
)

SCAN_SOURCE = (
    PROJECT
    / "scripts"
    / "100c_scan_band_drop_problematic.py"
)

PLOT_SOURCE = (
    PROJECT
    / "scripts"
    / "108_plot_B_BA_with_separate_focal_genotype_tracks.R"
)

BASE_OUTPUT = (
    PROJECT
    / "scripts"
    / "109b_base_scorer_B_heurippa_besckei.py"
)

SCAN_OUTPUT = (
    PROJECT
    / "scripts"
    / "109c_scan_B_heurippa_besckei.py"
)

PLOT_OUTPUT = (
    PROJECT
    / "scripts"
    / "109d_plot_B_heurippa_besckei_BA_and_focal.R"
)


# Cydno is deliberately absent.
EXCLUDED_TAXON_TERMS = (
    "numata",
    "pachinus",
    "hecale",
    "ismenius",
    "timareta",
)

MANUAL_EXCLUSIONS = {
    "ERR3653118",
    "ERR3653120",
    "ERR3653298",
    "ERR3653306",
    "ERR3653308",
    "ERR3653309",
    "ERR3650481",
    "ERR3652402",
}


def clean_series(series):
    return (
        series
        .fillna("")
        .astype(str)
        .str.strip()
    )


def as_bool(series):
    return (
        clean_series(series)
        .str.lower()
        .isin(
            {
                "true",
                "t",
                "1",
                "yes",
                "y",
            }
        )
    )


required_files = [
    SOURCE_MANIFEST,
    BASE_SOURCE,
    SCAN_SOURCE,
    PLOT_SOURCE,
]

missing_files = [
    str(path)
    for path in required_files
    if not path.exists()
]

if missing_files:
    raise RuntimeError(
        "Missing required source files:\n"
        + "\n".join(missing_files)
    )


OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


manifest = pd.read_csv(
    SOURCE_MANIFEST,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

required_columns = {
    "cohort",
    "module",
    "sample",
    "taxon_name",
    "data_type",
    "passes_callability",
    "excluded_taxon",
    "keep",
}

missing_columns = (
    required_columns
    - set(manifest.columns)
)

if missing_columns:
    raise RuntimeError(
        "Source manifest is missing columns: "
        + ", ".join(
            sorted(missing_columns)
        )
    )


taxon_lower = (
    clean_series(
        manifest["taxon_name"]
    )
    .str.lower()
)

sample = clean_series(
    manifest["sample"]
)

data_type = clean_series(
    manifest["data_type"]
)

cohort = clean_series(
    manifest["cohort"]
)

passes_callability = as_bool(
    manifest["passes_callability"]
)


excluded_taxon = pd.Series(
    False,
    index=manifest.index,
)

for term in EXCLUDED_TAXON_TERMS:
    excluded_taxon = (
        excluded_taxon
        | taxon_lower.str.contains(
            term,
            regex=False,
            na=False,
        )
    )


manual_excluded = sample.isin(
    MANUAL_EXCLUSIONS
)

cohort_data_type_ok = (
    (
        (cohort == "WGS_only")
        & (data_type == "WGS")
    )
    |
    (
        (cohort == "WGS_plus_CS")
        & data_type.isin(
            {
                "WGS",
                "CS",
            }
        )
    )
)

new_keep = (
    cohort_data_type_ok
    & passes_callability
    & ~excluded_taxon
    & ~manual_excluded
)


manifest["excluded_taxon"] = (
    excluded_taxon
)

manifest["keep"] = new_keep

manifest.to_csv(
    OUTPUT_MANIFEST,
    sep="\t",
    index=False,
)


# ------------------------------------------------------------
# Clone and redirect the existing corrected base scorer.
# ------------------------------------------------------------

base_text = BASE_SOURCE.read_text()

base_text = base_text.replace(
    '"drop_problematic_reanalysis"',
    '"B_heurippa_besckei_reanalysis"',
)

base_text = base_text.replace(
    "original502_filtered_cohort_manifest."
    "WGS_vs_WGSplusCS.tsv",
    "original502_filtered_cohort_manifest."
    "B_heurippa_besckei.WGS_vs_WGSplusCS.tsv",
)

new_exclusion_block = '''EXCLUDED_TAXON_TERMS = (
    "numata",
    "pachinus",
    "hecale",
    "ismenius",
    "timareta",
)'''

base_text, exclusion_replacements = re.subn(
    r"EXCLUDED_TAXON_TERMS\s*=\s*\(.*?\n\)",
    new_exclusion_block,
    base_text,
    count=1,
    flags=re.DOTALL,
)

if exclusion_replacements != 1:
    raise RuntimeError(
        "Could not replace the base scorer's "
        "taxon-exclusion tuple."
    )

BASE_OUTPUT.write_text(
    base_text
)


# ------------------------------------------------------------
# Clone the broad Band scanner and redirect its base scorer
# and output directory.
# ------------------------------------------------------------

scan_text = SCAN_SOURCE.read_text()

scan_text = scan_text.replace(
    "100b_base_scorer_drop_problematic.py",
    "109b_base_scorer_B_heurippa_besckei.py",
)

scan_text = scan_text.replace(
    '"drop_problematic_reanalysis"',
    '"B_heurippa_besckei_reanalysis"',
)

SCAN_OUTPUT.write_text(
    scan_text
)


# ------------------------------------------------------------
# Clone the exact-only Band plot into the new analysis branch.
# ------------------------------------------------------------

plot_text = PLOT_SOURCE.read_text()

plot_text = plot_text.replace(
    '"drop_problematic_reanalysis"',
    '"B_heurippa_besckei_reanalysis"',
)

plot_text = plot_text.replace(
    '"broad_B_BA_manhattan"',
    '"broad_B_heurippa_besckei_BA_manhattan"',
)

plot_text = plot_text.replace(
    '"Band association and pooled focal genotypes"',
    '"Band association with cydno restored and timareta excluded"',
)

PLOT_OUTPUT.write_text(
    plot_text
)


for script in [
    BASE_OUTPUT,
    SCAN_OUTPUT,
]:
    py_compile.compile(
        str(script),
        doraise=True,
    )


# ------------------------------------------------------------
# Cohort audit
# ------------------------------------------------------------

band = manifest[
    (
        manifest["module"] == "Band"
    )
].copy()

retained = band[
    as_bool(
        band["keep"]
    )
].copy()

print(
    "BAND HEURIPPA/BESCKEI COHORT PREPARATION"
)
print()

print("FILTERING RULES")
print("Cydno restored: True")
print("All timareta excluded: True")
print(
    "Other taxon exclusions: "
    "numata,pachinus,hecale,ismenius"
)
print(
    "Manual exclusions: "
    + ",".join(
        sorted(
            MANUAL_EXCLUSIONS
        )
    )
)

print()
print("RETAINED BAND COHORTS")
print(
    retained.groupby(
        [
            "cohort",
            "data_type",
        ]
    )["sample"]
    .nunique()
    .reset_index(
        name="samples"
    )
    .to_string(
        index=False
    )
)

audit_groups = {
    "heurippa": taxon_lower.str.contains(
        "heurippa",
        regex=False,
        na=False,
    ),
    "besckei": taxon_lower.str.contains(
        "besckei",
        regex=False,
        na=False,
    ),
    "cydno": taxon_lower.str.contains(
        "cydno",
        regex=False,
        na=False,
    ),
    "timareta": taxon_lower.str.contains(
        "timareta",
        regex=False,
        na=False,
    ),
    "numata": taxon_lower.str.contains(
        "numata",
        regex=False,
        na=False,
    ),
    "pachinus": taxon_lower.str.contains(
        "pachinus",
        regex=False,
        na=False,
    ),
}

print()
print(
    "TARGET-TAXON RETENTION IN "
    "WGS_PLUS_CS BAND COHORT"
)
print(
    "group\ttotal_samples\tretained_samples"
)

band_wgscs = (
    (manifest["module"] == "Band")
    & (manifest["cohort"] == "WGS_plus_CS")
)

for group, group_mask in audit_groups.items():
    subset_mask = (
        band_wgscs
        & group_mask
    )

    total_n = manifest.loc[
        subset_mask,
        "sample",
    ].nunique()

    retained_n = manifest.loc[
        subset_mask
        & as_bool(
            manifest["keep"]
        ),
        "sample",
    ].nunique()

    print(
        f"{group}\t"
        f"{total_n}\t"
        f"{retained_n}"
    )


manual_survivors = sorted(
    set(
        manifest.loc[
            band_wgscs
            & as_bool(
                manifest["keep"]
            ),
            "sample",
        ]
    )
    & MANUAL_EXCLUSIONS
)

cydno_retained = manifest.loc[
    band_wgscs
    & as_bool(
        manifest["keep"]
    )
    & audit_groups["cydno"],
    "sample",
].nunique()

timareta_retained = manifest.loc[
    band_wgscs
    & as_bool(
        manifest["keep"]
    )
    & audit_groups["timareta"],
    "sample",
].nunique()

if cydno_retained == 0:
    raise RuntimeError(
        "No cydno samples survived the new Band cohort."
    )

if timareta_retained != 0:
    raise RuntimeError(
        f"{timareta_retained} timareta samples "
        "incorrectly survived."
    )

if manual_survivors:
    raise RuntimeError(
        "Manual exclusions survived: "
        + ",".join(
            manual_survivors
        )
    )


print()
print(
    f"Manifest: {OUTPUT_MANIFEST}"
)
print(
    f"Base scorer: {BASE_OUTPUT}"
)
print(
    f"Band scanner: {SCAN_OUTPUT}"
)
print(
    f"Band plotter: {PLOT_OUTPUT}"
)
print()
print(
    "The existing D/R analysis and outputs "
    "were not modified."
)
