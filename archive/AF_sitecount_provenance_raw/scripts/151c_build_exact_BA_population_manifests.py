#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

OUTDIR = (
    PROJECT
    / "analysis/figures"
    / "key_regions_population_context_genoplots"
)

DR_LONG = (
    PROJECT
    / "analysis/drop_problematic_reanalysis/checkpoints"
    / "perfect_D_R_all_comparative_genotypes.long.tsv"
)

BAND_MANIFEST = (
    PROJECT
    / "analysis/B_western_postman_vs_western_Amazon_DR_reanalysis"
    / "checkpoints"
    / "B_western_retained_sample_manifest.tsv"
)

OUT_SAMPLES = (
    OUTDIR
    / "exact_BA_comparative_sample_manifest.tsv"
)

OUT_POPS = (
    OUTDIR
    / "exact_BA_population_summary.tsv"
)


def clean_text(x):
    x = x.fillna("").astype(str)
    return (
        x.str.replace("\n", " ", regex=False)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )


def normalise_country(x):
    x = clean_text(x)

    # Display-only repair of the UTF-8 mojibake seen in the
    # terminal output. This does not alter sample selection.
    x = x.replace(
        {
            "PanamÃ¡": "Panamá",
            "Panama": "Panamá",
        }
    )

    return x


def population_label(taxon, country):
    if country:
        return f"{taxon} | {country}"
    return f"{taxon} | country unavailable"


OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# DENNIS / RAY
# ============================================================

if not DR_LONG.is_file():
    raise RuntimeError(
        f"Missing final D/R comparative genotype table: {DR_LONG}"
    )

dr = pd.read_csv(
    DR_LONG,
    sep="\t",
    dtype=str,
    keep_default_na=False,
    low_memory=False,
)

required_dr = {
    "module",
    "sample",
    "phenotype",
    "taxon",
    "country",
    "data_type",
}

missing_dr = required_dr - set(dr.columns)

if missing_dr:
    raise RuntimeError(
        "D/R table is missing columns: "
        + ", ".join(sorted(missing_dr))
    )

dr["module"] = clean_text(dr["module"])
dr["sample"] = clean_text(dr["sample"])
dr["taxon"] = clean_text(dr["taxon"])
dr["country"] = normalise_country(dr["country"])
dr["data_type"] = clean_text(dr["data_type"])

dr_samples = (
    dr[
        dr["module"].isin(["Dennis", "Ray"])
    ][
        [
            "module",
            "sample",
            "phenotype",
            "taxon",
            "country",
            "data_type",
        ]
    ]
    .drop_duplicates()
    .copy()
)

# A sample must have one unambiguous metadata row within module.
dup = (
    dr_samples
    .groupby(
        ["module", "sample"],
        dropna=False,
    )
    .size()
)

bad_dup = dup[dup != 1]

if not bad_dup.empty:
    raise RuntimeError(
        "D/R samples have conflicting metadata rows:\n"
        + bad_dup.to_string()
    )

dr_samples["BA_pattern"] = dr_samples["module"]

dr_samples = dr_samples[
    [
        "BA_pattern",
        "sample",
        "phenotype",
        "taxon",
        "country",
        "data_type",
    ]
]


# ============================================================
# BAND
# ============================================================

if not BAND_MANIFEST.is_file():
    raise RuntimeError(
        f"Missing final Band retained manifest: {BAND_MANIFEST}"
    )

band = pd.read_csv(
    BAND_MANIFEST,
    sep="\t",
    dtype=str,
    keep_default_na=False,
    low_memory=False,
)

required_band = {
    "cohort",
    "sample",
    "taxon_name",
    "country",
    "data_type",
}

missing_band = required_band - set(band.columns)

if missing_band:
    raise RuntimeError(
        "Band manifest is missing columns: "
        + ", ".join(sorted(missing_band))
    )

band["cohort"] = clean_text(band["cohort"])
band["sample"] = clean_text(band["sample"])
band["taxon_name"] = clean_text(band["taxon_name"])
band["country"] = normalise_country(band["country"])
band["data_type"] = clean_text(band["data_type"])

band_samples = (
    band[
        band["cohort"].eq("WGS_plus_CS")
    ][
        [
            "sample",
            "taxon_name",
            "country",
            "data_type",
        ]
    ]
    .drop_duplicates()
    .copy()
)

band_samples["BA_pattern"] = "Band"
band_samples["phenotype"] = ""

band_samples = band_samples.rename(
    columns={
        "taxon_name": "taxon",
    }
)

band_samples = band_samples[
    [
        "BA_pattern",
        "sample",
        "phenotype",
        "taxon",
        "country",
        "data_type",
    ]
]


# ============================================================
# COMBINE + HARD VALIDATION
# ============================================================

samples = pd.concat(
    [
        band_samples,
        dr_samples,
    ],
    ignore_index=True,
)

counts = (
    samples
    .groupby("BA_pattern")["sample"]
    .nunique()
    .reindex(
        ["Band", "Ray", "Dennis"]
    )
)

expected = {
    "Band": 163,
    "Ray": 368,
    "Dennis": 369,
}

problems = []

for pattern, expected_n in expected.items():
    observed_n = int(counts.loc[pattern])

    if observed_n != expected_n:
        problems.append(
            f"{pattern}: expected {expected_n}, "
            f"recovered {observed_n}"
        )

if problems:
    raise RuntimeError(
        "Exact BA cohort validation failed:\n"
        + "\n".join(problems)
    )

focal_names = {
    "TN_Left",
    "TN_Right",
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
}

bad_focal = samples[
    samples["sample"].isin(focal_names)
]

if not bad_focal.empty:
    raise RuntimeError(
        "Focal samples leaked into a comparative BA cohort:\n"
        + bad_focal.to_string(index=False)
    )


# Population identity is assigned only AFTER exact BA
# sample selection.
samples["population"] = [
    population_label(taxon, country)
    for taxon, country in zip(
        samples["taxon"],
        samples["country"],
    )
]

samples = samples.sort_values(
    [
        "BA_pattern",
        "taxon",
        "country",
        "sample",
    ]
).reset_index(drop=True)

samples.to_csv(
    OUT_SAMPLES,
    sep="\t",
    index=False,
    na_rep="NA",
)


# ============================================================
# POPULATION SUMMARY
# ============================================================

pops = (
    samples
    .groupby(
        [
            "BA_pattern",
            "population",
            "taxon",
            "country",
        ],
        dropna=False,
    )
    .agg(
        samples=("sample", "nunique"),
        WGS=("data_type", lambda x: (x == "WGS").sum()),
        CS=("data_type", lambda x: (x == "CS").sum()),
    )
    .reset_index()
    .sort_values(
        [
            "BA_pattern",
            "taxon",
            "country",
        ]
    )
)

pops.to_csv(
    OUT_POPS,
    sep="\t",
    index=False,
    na_rep="NA",
)


print("EXACT BA COMPARATIVE SAMPLE COUNTS")
print(
    counts.to_string()
)

print(
    "\nPOPULATIONS REPRESENTED IN EACH BA COHORT"
)

for pattern in ["Band", "Ray", "Dennis"]:
    print(
        f"\n{pattern}"
    )

    print(
        pops[
            pops["BA_pattern"].eq(pattern)
        ][
            [
                "population",
                "samples",
                "WGS",
                "CS",
            ]
        ].to_string(
            index=False
        )
    )

print(
    "\nSample manifest:",
    OUT_SAMPLES,
)

print(
    "Population summary:",
    OUT_POPS,
)
