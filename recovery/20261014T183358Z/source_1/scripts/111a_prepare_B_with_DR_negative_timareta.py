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

METADATA = (
    PROJECT
    / "data"
    / "derived"
    / "region_panel_metadata.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis"
    / "B_cydno_DRnegative_timareta_reanalysis"
    / "checkpoints"
)

OUTPUT_MANIFEST = (
    OUTDIR
    / "original502_filtered_cohort_manifest."
    "B_cydno_DRnegative_timareta."
    "WGS_vs_WGSplusCS.tsv"
)

DROPPED_DR_TABLE = (
    OUTDIR
    / "B_samples_dropped_for_D_or_R_positive.tsv"
)

TIMARETA_AUDIT_TABLE = (
    OUTDIR
    / "B_timareta_retention_audit.tsv"
)

BASE_SOURCE = (
    PROJECT
    / "scripts"
    / "109b_base_scorer_B_heurippa_besckei.py"
)

SCAN_SOURCE = (
    PROJECT
    / "scripts"
    / "109c_scan_B_heurippa_besckei.py"
)

PLOT_SOURCE = (
    PROJECT
    / "scripts"
    / "109d_plot_B_heurippa_besckei_BA_and_focal.R"
)

BASE_OUTPUT = (
    PROJECT
    / "scripts"
    / "111b_base_scorer_B_DRnegative_timareta.py"
)

SCAN_OUTPUT = (
    PROJECT
    / "scripts"
    / "111c_scan_B_DRnegative_timareta.py"
)

PLOT_OUTPUT = (
    PROJECT
    / "scripts"
    / "111d_plot_B_DRnegative_timareta_BA_and_focal.R"
)


# Cydno and timareta are deliberately absent from this list.
BASE_EXCLUDED_TAXON_TERMS = (
    "numata",
    "pachinus",
    "hecale",
    "ismenius",
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


def numeric_phenotype(series):
    return pd.to_numeric(
        clean_series(series),
        errors="coerce",
    )


required_files = [
    SOURCE_MANIFEST,
    METADATA,
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
        "Missing required files:\n"
        + "\n".join(
            missing_files
        )
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

metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


required_manifest_columns = {
    "cohort",
    "module",
    "sample",
    "taxon_name",
    "data_type",
    "passes_callability",
    "excluded_taxon",
    "keep",
}

missing_manifest_columns = (
    required_manifest_columns
    - set(manifest.columns)
)

if missing_manifest_columns:
    raise RuntimeError(
        "Manifest is missing columns: "
        + ", ".join(
            sorted(
                missing_manifest_columns
            )
        )
    )


required_metadata_columns = {
    "Accession",
    "Taxon Name",
    "B",
    "D",
    "R",
}

missing_metadata_columns = (
    required_metadata_columns
    - set(metadata.columns)
)

if missing_metadata_columns:
    raise RuntimeError(
        "Metadata is missing columns: "
        + ", ".join(
            sorted(
                missing_metadata_columns
            )
        )
    )


metadata = metadata[
    [
        "Accession",
        "Taxon Name",
        "B",
        "D",
        "R",
    ]
].copy()

metadata["Accession"] = clean_series(
    metadata["Accession"]
)

metadata = metadata[
    metadata["Accession"] != ""
].copy()

metadata = metadata.drop_duplicates(
    subset=[
        "Accession",
    ],
    keep="first",
)


manifest = manifest.merge(
    metadata,
    how="left",
    left_on="sample",
    right_on="Accession",
    validate="many_to_one",
)


sample = clean_series(
    manifest["sample"]
)

taxon_lower = (
    clean_series(
        manifest["taxon_name"]
    )
    .str.lower()
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

D_numeric = numeric_phenotype(
    manifest["D"]
)

R_numeric = numeric_phenotype(
    manifest["R"]
)

B_numeric = numeric_phenotype(
    manifest["B"]
)

D_positive = (
    D_numeric == 1
)

R_positive = (
    R_numeric == 1
)

DR_positive = (
    D_positive
    | R_positive
)

explicit_DR_negative = (
    (D_numeric == 0)
    & (R_numeric == 0)
)

is_timareta = taxon_lower.str.contains(
    "timareta",
    regex=False,
    na=False,
)


base_excluded_taxon = pd.Series(
    False,
    index=manifest.index,
)

for term in BASE_EXCLUDED_TAXON_TERMS:
    base_excluded_taxon = (
        base_excluded_taxon
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


# Timareta are eligible only when explicitly D=0 and R=0.
timareta_state_ok = (
    ~is_timareta
    | explicit_DR_negative
)


new_keep = (
    cohort_data_type_ok
    & passes_callability
    & ~base_excluded_taxon
    & ~manual_excluded
    & ~DR_positive
    & timareta_state_ok
)


manifest["excluded_taxon"] = (
    base_excluded_taxon
)

manifest["keep"] = new_keep


original_columns = list(
    pd.read_csv(
        SOURCE_MANIFEST,
        sep="\t",
        nrows=0,
    ).columns
)

manifest[
    original_columns
].to_csv(
    OUTPUT_MANIFEST,
    sep="\t",
    index=False,
)


band_wgscs = (
    (manifest["module"] == "Band")
    & (
        manifest["cohort"]
        == "WGS_plus_CS"
    )
)

new_retained = manifest.loc[
    band_wgscs
    & as_bool(
        manifest["keep"]
    )
].drop_duplicates(
    subset=[
        "sample",
    ]
)


dropped_DR = manifest.loc[
    band_wgscs
    & cohort_data_type_ok
    & passes_callability
    & ~base_excluded_taxon
    & ~manual_excluded
    & DR_positive
].drop_duplicates(
    subset=[
        "sample",
    ]
).copy()

dropped_DR[
    [
        "sample",
        "taxon_name",
        "data_type",
        "B",
        "D",
        "R",
    ]
].sort_values(
    [
        "taxon_name",
        "sample",
    ]
).to_csv(
    DROPPED_DR_TABLE,
    sep="\t",
    index=False,
)


timareta_audit = manifest.loc[
    band_wgscs
    & is_timareta
].drop_duplicates(
    subset=[
        "sample",
    ]
).copy()

timareta_audit["D_numeric"] = D_numeric.loc[
    timareta_audit.index
]

timareta_audit["R_numeric"] = R_numeric.loc[
    timareta_audit.index
]

timareta_audit["explicit_DR_negative"] = (
    explicit_DR_negative.loc[
        timareta_audit.index
    ]
)

timareta_audit["manual_excluded"] = (
    manual_excluded.loc[
        timareta_audit.index
    ]
)

timareta_audit["retained"] = as_bool(
    timareta_audit["keep"]
)

timareta_audit[
    [
        "sample",
        "taxon_name",
        "data_type",
        "B",
        "D",
        "R",
        "explicit_DR_negative",
        "manual_excluded",
        "passes_callability",
        "retained",
    ]
].sort_values(
    [
        "taxon_name",
        "sample",
    ]
).to_csv(
    TIMARETA_AUDIT_TABLE,
    sep="\t",
    index=False,
)


surviving_DR_positive = new_retained.loc[
    (
        numeric_phenotype(
            new_retained["D"]
        ) == 1
    )
    |
    (
        numeric_phenotype(
            new_retained["R"]
        ) == 1
    ),
    "sample",
].tolist()

if surviving_DR_positive:
    raise RuntimeError(
        "D-positive or R-positive samples survived: "
        + ",".join(
            sorted(
                surviving_DR_positive
            )
        )
    )


retained_timareta = new_retained.loc[
    clean_series(
        new_retained["taxon_name"]
    )
    .str.lower()
    .str.contains(
        "timareta",
        regex=False,
        na=False,
    )
].copy()

invalid_retained_timareta = retained_timareta.loc[
    ~(
        (
            numeric_phenotype(
                retained_timareta["D"]
            ) == 0
        )
        & (
            numeric_phenotype(
                retained_timareta["R"]
            ) == 0
        )
    ),
    "sample",
].tolist()

if invalid_retained_timareta:
    raise RuntimeError(
        "Timareta without explicit D=0,R=0 survived: "
        + ",".join(
            sorted(
                invalid_retained_timareta
            )
        )
    )


retained_taxon_lower = (
    clean_series(
        new_retained["taxon_name"]
    )
    .str.lower()
)

target_counts = {}

for term in [
    "heurippa",
    "besckei",
    "cydno",
    "timareta",
]:
    target_counts[term] = new_retained.loc[
        retained_taxon_lower.str.contains(
            term,
            regex=False,
            na=False,
        ),
        "sample",
    ].nunique()


for required_group in [
    "heurippa",
    "besckei",
    "cydno",
]:
    if target_counts[required_group] == 0:
        raise RuntimeError(
            f"No {required_group} samples survived."
        )


B_counts = (
    new_retained.assign(
        B_state=numeric_phenotype(
            new_retained["B"]
        )
    )
    .groupby(
        "B_state",
        dropna=False,
    )["sample"]
    .nunique()
    .reset_index(
        name="samples"
    )
)

explicit_B_states = set(
    B_counts["B_state"]
    .dropna()
    .tolist()
)

if not {
    0.0,
    1.0,
}.issubset(
    explicit_B_states
):
    raise RuntimeError(
        "Both B=0 and B=1 classes were not retained."
    )


# ------------------------------------------------------------
# Clone scorer, scanner and exact-only Band plotter.
# ------------------------------------------------------------

base_text = BASE_SOURCE.read_text()

base_text = base_text.replace(
    '"B_heurippa_besckei_reanalysis"',
    '"B_cydno_DRnegative_timareta_reanalysis"',
)

base_text = base_text.replace(
    "original502_filtered_cohort_manifest."
    "B_heurippa_besckei.WGS_vs_WGSplusCS.tsv",
    "original502_filtered_cohort_manifest."
    "B_cydno_DRnegative_timareta."
    "WGS_vs_WGSplusCS.tsv",
)

new_exclusion_block = '''EXCLUDED_TAXON_TERMS = (
    "numata",
    "pachinus",
    "hecale",
    "ismenius",
)'''

base_text, replacement_n = re.subn(
    r"EXCLUDED_TAXON_TERMS\s*=\s*\(.*?\n\)",
    new_exclusion_block,
    base_text,
    count=1,
    flags=re.DOTALL,
)

if replacement_n != 1:
    raise RuntimeError(
        "Could not update scorer taxon exclusions."
    )

BASE_OUTPUT.write_text(
    base_text
)


scan_text = SCAN_SOURCE.read_text()

scan_text = scan_text.replace(
    "109b_base_scorer_B_heurippa_besckei.py",
    "111b_base_scorer_B_DRnegative_timareta.py",
)

scan_text = scan_text.replace(
    '"B_heurippa_besckei_reanalysis"',
    '"B_cydno_DRnegative_timareta_reanalysis"',
)

SCAN_OUTPUT.write_text(
    scan_text
)


plot_text = PLOT_SOURCE.read_text()

plot_text = plot_text.replace(
    '"B_heurippa_besckei_reanalysis"',
    '"B_cydno_DRnegative_timareta_reanalysis"',
)

plot_text = plot_text.replace(
    '"broad_B_heurippa_besckei_BA_manhattan"',
    '"broad_B_cydno_DRnegative_timareta_BA_manhattan"',
)

plot_text = plot_text.replace(
    '"Band association with cydno restored and timareta excluded"',
    '"Band association with D/R-negative timareta restored"',
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


print(
    "BAND COHORT WITH D/R-NEGATIVE TIMARETA RESTORED"
)
print()

print(
    "Rule: all explicit D=1 or R=1 samples are excluded."
)
print(
    "Timareta are retained only when D=0 and R=0."
)
print(
    "Timareta with missing D or R are excluded."
)

print()
print(
    "Retained WGS+CS Band samples: "
    f"{new_retained['sample'].nunique()}"
)

print(
    "D-positive or R-positive samples removed: "
    f"{dropped_DR['sample'].nunique()}"
)

print()
print("RETAINED TARGET GROUPS")
for term in [
    "heurippa",
    "besckei",
    "cydno",
    "timareta",
]:
    print(
        f"{term}\t{target_counts[term]}"
    )

print()
print("RETAINED TIMARETA BY TAXON")

if retained_timareta.empty:
    print("None")
else:
    print(
        retained_timareta.groupby(
            "taxon_name"
        )["sample"]
        .nunique()
        .reset_index(
            name="samples"
        )
        .to_string(
            index=False
        )
    )

print()
print("RETAINED BAND PHENOTYPE COUNTS")
print(
    B_counts.to_string(
        index=False
    )
)

print()
print(f"Manifest: {OUTPUT_MANIFEST}")
print(f"Dropped D/R-positive samples: {DROPPED_DR_TABLE}")
print(f"Timareta audit: {TIMARETA_AUDIT_TABLE}")
print(f"Base scorer: {BASE_OUTPUT}")
print(f"Band scanner: {SCAN_OUTPUT}")
print(f"Plotter: {PLOT_OUTPUT}")
print()
print(
    "All earlier D/R and Band branches remain unchanged."
)
