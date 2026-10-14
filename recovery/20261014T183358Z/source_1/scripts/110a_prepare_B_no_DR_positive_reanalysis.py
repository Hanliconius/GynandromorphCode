#!/usr/bin/env python3

from pathlib import Path
import py_compile

import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

SOURCE_DIR = (
    PROJECT
    / "analysis"
    / "B_heurippa_besckei_reanalysis"
    / "checkpoints"
)

SOURCE_MANIFEST = (
    SOURCE_DIR
    / "original502_filtered_cohort_manifest."
    "B_heurippa_besckei.WGS_vs_WGSplusCS.tsv"
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
    / "B_heurippa_besckei_no_DR_positive_reanalysis"
    / "checkpoints"
)

OUTPUT_MANIFEST = (
    OUTDIR
    / "original502_filtered_cohort_manifest."
    "B_heurippa_besckei_no_DR_positive."
    "WGS_vs_WGSplusCS.tsv"
)

DROPPED_TABLE = (
    OUTDIR
    / "B_samples_dropped_for_D_or_R_positive.tsv"
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
    / "110b_base_scorer_B_no_DR_positive.py"
)

SCAN_OUTPUT = (
    PROJECT
    / "scripts"
    / "110c_scan_B_no_DR_positive.py"
)

PLOT_OUTPUT = (
    PROJECT
    / "scripts"
    / "110d_plot_B_no_DR_positive_BA_and_focal.R"
)


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


def phenotype_positive(series):
    numeric = pd.to_numeric(
        clean_series(series),
        errors="coerce",
    )

    return numeric == 1


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


duplicate_accessions = metadata.loc[
    metadata["Accession"].duplicated(
        keep=False
    ),
    "Accession",
].unique()

if len(duplicate_accessions) > 0:
    conflicting = []

    for accession in duplicate_accessions:
        subset = metadata.loc[
            metadata["Accession"] == accession,
            [
                "B",
                "D",
                "R",
            ],
        ].drop_duplicates()

        if len(subset) > 1:
            conflicting.append(
                accession
            )

    if conflicting:
        raise RuntimeError(
            "Conflicting duplicate metadata accessions: "
            + ",".join(
                sorted(
                    conflicting
                )
            )
        )

    metadata = metadata.drop_duplicates(
        subset=[
            "Accession",
        ],
        keep="first",
    )


metadata["D_positive"] = phenotype_positive(
    metadata["D"]
)

metadata["R_positive"] = phenotype_positive(
    metadata["R"]
)

metadata["DR_positive"] = (
    metadata["D_positive"]
    | metadata["R_positive"]
)

metadata["B_numeric"] = pd.to_numeric(
    clean_series(
        metadata["B"]
    ),
    errors="coerce",
)


manifest = manifest.merge(
    metadata,
    how="left",
    left_on="sample",
    right_on="Accession",
    validate="many_to_one",
)

manifest["previous_keep"] = as_bool(
    manifest["keep"]
)

manifest["DR_positive"] = (
    manifest["DR_positive"]
    .fillna(False)
    .astype(bool)
)

manifest["keep"] = (
    manifest["previous_keep"]
    & ~manifest["DR_positive"]
)


output_columns = [
    column
    for column in pd.read_csv(
        SOURCE_MANIFEST,
        sep="\t",
        nrows=0,
    ).columns
]

manifest[
    output_columns
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

previous_retained = manifest.loc[
    band_wgscs
    & manifest["previous_keep"]
].drop_duplicates(
    subset=[
        "sample",
    ]
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

dropped = manifest.loc[
    band_wgscs
    & manifest["previous_keep"]
    & manifest["DR_positive"]
].drop_duplicates(
    subset=[
        "sample",
    ]
).copy()

dropped[
    [
        "sample",
        "taxon_name",
        "data_type",
        "B",
        "D",
        "R",
        "D_positive",
        "R_positive",
    ]
].sort_values(
    [
        "taxon_name",
        "sample",
    ]
).to_csv(
    DROPPED_TABLE,
    sep="\t",
    index=False,
)


# ------------------------------------------------------------
# Verify that no explicit D-positive or R-positive sample
# survives in the new Band cohort.
# ------------------------------------------------------------

surviving_DR_positive = new_retained.loc[
    new_retained["DR_positive"],
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


taxon_lower = (
    clean_series(
        new_retained["taxon_name"]
    )
    .str.lower()
)

heurippa_n = new_retained.loc[
    taxon_lower.str.contains(
        "heurippa",
        regex=False,
        na=False,
    ),
    "sample",
].nunique()

besckei_n = new_retained.loc[
    taxon_lower.str.contains(
        "besckei",
        regex=False,
        na=False,
    ),
    "sample",
].nunique()

cydno_n = new_retained.loc[
    taxon_lower.str.contains(
        "cydno",
        regex=False,
        na=False,
    ),
    "sample",
].nunique()

timareta_n = new_retained.loc[
    taxon_lower.str.contains(
        "timareta",
        regex=False,
        na=False,
    ),
    "sample",
].nunique()

if heurippa_n == 0:
    raise RuntimeError(
        "No heurippa samples survived."
    )

if besckei_n == 0:
    raise RuntimeError(
        "No besckei samples survived."
    )

if cydno_n == 0:
    raise RuntimeError(
        "No cydno samples survived."
    )

if timareta_n != 0:
    raise RuntimeError(
        f"{timareta_n} timareta samples survived."
    )


B_counts = (
    new_retained.assign(
        B_state=pd.to_numeric(
            clean_series(
                new_retained["B"]
            ),
            errors="coerce",
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
        "The no-D/R-positive cohort does not retain "
        "both B=0 and B=1 phenotype classes."
    )


# ------------------------------------------------------------
# Clone the working scorer into another independent branch.
# ------------------------------------------------------------

base_text = BASE_SOURCE.read_text()

base_text = base_text.replace(
    '"B_heurippa_besckei_reanalysis"',
    '"B_heurippa_besckei_no_DR_positive_reanalysis"',
)

base_text = base_text.replace(
    "original502_filtered_cohort_manifest."
    "B_heurippa_besckei.WGS_vs_WGSplusCS.tsv",
    "original502_filtered_cohort_manifest."
    "B_heurippa_besckei_no_DR_positive."
    "WGS_vs_WGSplusCS.tsv",
)

BASE_OUTPUT.write_text(
    base_text
)


scan_text = SCAN_SOURCE.read_text()

scan_text = scan_text.replace(
    "109b_base_scorer_B_heurippa_besckei.py",
    "110b_base_scorer_B_no_DR_positive.py",
)

scan_text = scan_text.replace(
    '"B_heurippa_besckei_reanalysis"',
    '"B_heurippa_besckei_no_DR_positive_reanalysis"',
)

SCAN_OUTPUT.write_text(
    scan_text
)


plot_text = PLOT_SOURCE.read_text()

plot_text = plot_text.replace(
    '"B_heurippa_besckei_reanalysis"',
    '"B_heurippa_besckei_no_DR_positive_reanalysis"',
)

plot_text = plot_text.replace(
    '"broad_B_heurippa_besckei_BA_manhattan"',
    '"broad_B_heurippa_besckei_no_DR_positive_BA_manhattan"',
)

plot_text = plot_text.replace(
    '"Band association with cydno restored and timareta excluded"',
    '"Band association after excluding D-positive and R-positive samples"',
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
    "BAND NO-D/R-POSITIVE COHORT PREPARATION"
)
print()

print(
    "Interpretation: DR-positive means D=1 OR R=1."
)
print()

print(
    "Previous WGS+CS retained samples: "
    f"{previous_retained['sample'].nunique()}"
)

print(
    "New WGS+CS retained samples: "
    f"{new_retained['sample'].nunique()}"
)

print(
    "D-positive or R-positive samples removed: "
    f"{dropped['sample'].nunique()}"
)

print()
print("REMOVED PHENOTYPE COMBINATIONS")

if dropped.empty:
    print("None")
else:
    print(
        dropped.groupby(
            [
                "D",
                "R",
            ],
            dropna=False,
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
print("RETAINED TARGET GROUPS")
print(
    f"heurippa\t{heurippa_n}"
)
print(
    f"besckei\t{besckei_n}"
)
print(
    f"cydno\t{cydno_n}"
)
print(
    f"timareta\t{timareta_n}"
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
print(f"Dropped samples: {DROPPED_TABLE}")
print(f"Base scorer: {BASE_OUTPUT}")
print(f"Band scanner: {SCAN_OUTPUT}")
print(f"Plotter: {PLOT_OUTPUT}")
print()
print(
    "Existing D/R and previous Band branches "
    "were not modified."
)
