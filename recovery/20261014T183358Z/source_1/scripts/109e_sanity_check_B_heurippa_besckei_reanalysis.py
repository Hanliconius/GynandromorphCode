#!/usr/bin/env python3

from pathlib import Path
import hashlib

import numpy as np
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

OLD_DIR = (
    PROJECT
    / "analysis"
    / "drop_problematic_reanalysis"
    / "checkpoints"
)

NEW_DIR = (
    PROJECT
    / "analysis"
    / "B_heurippa_besckei_reanalysis"
    / "checkpoints"
)

OLD_MANIFEST = (
    OLD_DIR
    / "original502_filtered_cohort_manifest."
    "WGS_vs_WGSplusCS.tsv"
)

NEW_MANIFEST = (
    NEW_DIR
    / "original502_filtered_cohort_manifest."
    "B_heurippa_besckei.WGS_vs_WGSplusCS.tsv"
)

SCORE_NAME = (
    "longBand_original502_allSNP_scores."
    "filtered.WGS_vs_WGSplusCS.tsv"
)

OLD_SCORES = (
    OLD_DIR
    / SCORE_NAME
)

NEW_SCORES = (
    NEW_DIR
    / SCORE_NAME
)

OLD_PDF = (
    PROJECT
    / "analysis"
    / "figures"
    / "drop_problematic_reanalysis"
    / "broad_B_BA_manhattan"
    / "broad_B_DP2_GQ20_BA_and_focal_genotypes.pdf"
)

NEW_PDF = (
    PROJECT
    / "analysis"
    / "figures"
    / "B_heurippa_besckei_reanalysis"
    / "broad_B_heurippa_besckei_BA_manhattan"
    / "broad_B_DP2_GQ20_BA_and_focal_genotypes.pdf"
)

PLOT_SCRIPT = (
    PROJECT
    / "scripts"
    / "109d_plot_B_heurippa_besckei_BA_and_focal.R"
)

OUTPUT_DIR = NEW_DIR

ADDED_FILE = (
    OUTPUT_DIR
    / "B_heurippa_besckei_samples_added_vs_previous.tsv"
)

REMOVED_FILE = (
    OUTPUT_DIR
    / "B_heurippa_besckei_samples_removed_vs_previous.tsv"
)

SCORE_DIFFERENCE_FILE = (
    OUTPUT_DIR
    / "B_heurippa_besckei_BA_differences_vs_previous.tsv"
)


def as_bool(series):
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
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


def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(
                1024 * 1024
            )

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


required_files = [
    OLD_MANIFEST,
    NEW_MANIFEST,
    OLD_SCORES,
    NEW_SCORES,
    PLOT_SCRIPT,
]

missing = [
    str(path)
    for path in required_files
    if not path.exists()
]

if missing:
    raise RuntimeError(
        "Missing required files:\n"
        + "\n".join(missing)
    )


old_manifest = pd.read_csv(
    OLD_MANIFEST,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

new_manifest = pd.read_csv(
    NEW_MANIFEST,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)


def retained_band_samples(manifest):
    subset = manifest.loc[
        (
            manifest["module"] == "Band"
        )
        & (
            manifest["cohort"]
            == "WGS_plus_CS"
        )
        & as_bool(
            manifest["keep"]
        )
    ].copy()

    return subset.drop_duplicates(
        subset=[
            "sample",
        ]
    )


old_retained = retained_band_samples(
    old_manifest
)

new_retained = retained_band_samples(
    new_manifest
)

old_samples = set(
    old_retained["sample"]
)

new_samples = set(
    new_retained["sample"]
)

added_samples = sorted(
    new_samples - old_samples
)

removed_samples = sorted(
    old_samples - new_samples
)


added_table = new_retained[
    new_retained["sample"].isin(
        added_samples
    )
][
    [
        "sample",
        "taxon_name",
        "data_type",
        "callable_fraction",
    ]
].sort_values(
    [
        "taxon_name",
        "sample",
    ]
)

removed_table = old_retained[
    old_retained["sample"].isin(
        removed_samples
    )
][
    [
        "sample",
        "taxon_name",
        "data_type",
        "callable_fraction",
    ]
].sort_values(
    [
        "taxon_name",
        "sample",
    ]
)

added_table.to_csv(
    ADDED_FILE,
    sep="\t",
    index=False,
)

removed_table.to_csv(
    REMOVED_FILE,
    sep="\t",
    index=False,
)


def count_taxon_term(table, term):
    taxon = (
        table["taxon_name"]
        .fillna("")
        .astype(str)
        .str.lower()
    )

    return table.loc[
        taxon.str.contains(
            term,
            regex=False,
            na=False,
        ),
        "sample",
    ].nunique()


print(
    "BAND HEURIPPA/BESCKEI REANALYSIS SANITY CHECK"
)
print()

print("COHORT COMPARISON")
print(
    f"Previous retained WGS+CS Band samples: "
    f"{len(old_samples)}"
)
print(
    f"New retained WGS+CS Band samples: "
    f"{len(new_samples)}"
)
print(
    f"Samples added: {len(added_samples)}"
)
print(
    f"Samples removed: {len(removed_samples)}"
)

print()
print(
    "group\tprevious_retained\tnew_retained"
)

for term in [
    "heurippa",
    "besckei",
    "cydno",
    "timareta",
    "numata",
    "pachinus",
    "hecale",
    "ismenius",
]:
    print(
        f"{term}\t"
        f"{count_taxon_term(old_retained, term)}\t"
        f"{count_taxon_term(new_retained, term)}"
    )


print()
print("ADDED SAMPLES BY TAXON")

if added_table.empty:
    print("None")
else:
    print(
        added_table.groupby(
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
print("REMOVED SAMPLES BY TAXON")

if removed_table.empty:
    print("None")
else:
    print(
        removed_table.groupby(
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


def load_primary_scores(path):
    d = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "analysis_cohort",
        "scheme",
        "chrom",
        "pos",
        "ref",
        "alt",
        "balanced_accuracy",
        "score_eligible",
    }

    missing_columns = (
        required
        - set(d.columns)
    )

    if missing_columns:
        raise RuntimeError(
            f"{path.name} is missing columns: "
            + ",".join(
                sorted(
                    missing_columns
                )
            )
        )

    d = d.loc[
        (
            d["analysis_cohort"]
            == "WGS_plus_CS"
        )
        & (
            d["scheme"]
            == "DP2_GQ20"
        )
    ].copy()

    d["balanced_accuracy"] = pd.to_numeric(
        d["balanced_accuracy"],
        errors="coerce",
    )

    d["score_eligible_bool"] = as_bool(
        d["score_eligible"]
    )

    d["site_key"] = (
        d["chrom"]
        + ":"
        + d["pos"]
        + ":"
        + d["ref"]
        + ":"
        + d["alt"]
    )

    return d.drop_duplicates(
        subset=[
            "site_key",
        ]
    )


old_scores = load_primary_scores(
    OLD_SCORES
)

new_scores = load_primary_scores(
    NEW_SCORES
)

comparison = old_scores[
    [
        "site_key",
        "chrom",
        "pos",
        "ref",
        "alt",
        "balanced_accuracy",
        "score_eligible_bool",
    ]
].merge(
    new_scores[
        [
            "site_key",
            "balanced_accuracy",
            "score_eligible_bool",
        ]
    ],
    on="site_key",
    how="outer",
    suffixes=(
        "_previous",
        "_new",
    ),
    indicator=True,
)

comparison["BA_change"] = (
    comparison[
        "balanced_accuracy_new"
    ]
    - comparison[
        "balanced_accuracy_previous"
    ]
)

comparison["absolute_BA_change"] = (
    comparison["BA_change"].abs()
)

comparison["eligibility_changed"] = (
    comparison[
        "score_eligible_bool_previous"
    ]
    != comparison[
        "score_eligible_bool_new"
    ]
)

comparison = comparison.sort_values(
    [
        "absolute_BA_change",
        "site_key",
    ],
    ascending=[
        False,
        True,
    ],
    na_position="last",
)

comparison.to_csv(
    SCORE_DIFFERENCE_FILE,
    sep="\t",
    index=False,
)


shared = comparison.loc[
    comparison["_merge"] == "both"
].copy()

changed_BA = shared.loc[
    shared["absolute_BA_change"]
    > 1e-12
]

print()
print("SCORE-TABLE COMPARISON")
print(
    f"Previous score rows: {len(old_scores)}"
)
print(
    f"New score rows: {len(new_scores)}"
)
print(
    f"Shared physical SNPs: {len(shared)}"
)
print(
    f"SNPs with BA change > 1e-12: "
    f"{len(changed_BA)}"
)
print(
    f"SNPs changing score eligibility: "
    f"{int(shared['eligibility_changed'].sum())}"
)

if len(shared):
    print(
        "Median absolute BA change: "
        f"{shared['absolute_BA_change'].median():.8f}"
    )
    print(
        "Maximum absolute BA change: "
        f"{shared['absolute_BA_change'].max():.8f}"
    )

print(
    "Previous maximum eligible BA: "
    f"{old_scores.loc[old_scores['score_eligible_bool'], 'balanced_accuracy'].max():.8f}"
)

print(
    "New maximum eligible BA: "
    f"{new_scores.loc[new_scores['score_eligible_bool'], 'balanced_accuracy'].max():.8f}"
)

print(
    "Previous eligible BA >= 0.80: "
    f"{int((old_scores['score_eligible_bool'] & (old_scores['balanced_accuracy'] >= 0.80)).sum())}"
)

print(
    "New eligible BA >= 0.80: "
    f"{int((new_scores['score_eligible_bool'] & (new_scores['balanced_accuracy'] >= 0.80)).sum())}"
)


print()
print("TEN LARGEST BA CHANGES")

display_columns = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "balanced_accuracy_previous",
    "balanced_accuracy_new",
    "BA_change",
    "score_eligible_bool_previous",
    "score_eligible_bool_new",
]

print(
    comparison[
        display_columns
    ]
    .head(10)
    .to_string(
        index=False
    )
)


print()
print("FILE IDENTITY")

old_score_hash = sha256(
    OLD_SCORES
)

new_score_hash = sha256(
    NEW_SCORES
)

print(
    f"Previous score SHA256: {old_score_hash}"
)
print(
    f"New score SHA256:      {new_score_hash}"
)
print(
    "Score tables byte-identical: "
    f"{old_score_hash == new_score_hash}"
)


if OLD_PDF.exists() and NEW_PDF.exists():
    old_pdf_hash = sha256(
        OLD_PDF
    )

    new_pdf_hash = sha256(
        NEW_PDF
    )

    print()
    print(
        f"Previous PDF SHA256: {old_pdf_hash}"
    )
    print(
        f"New PDF SHA256:      {new_pdf_hash}"
    )
    print(
        "PDFs byte-identical: "
        f"{old_pdf_hash == new_pdf_hash}"
    )

    print(
        f"Previous PDF bytes: {OLD_PDF.stat().st_size}"
    )
    print(
        f"New PDF bytes:      {NEW_PDF.stat().st_size}"
    )
else:
    print()
    print("PDF comparison unavailable.")

    if not OLD_PDF.exists():
        print(
            f"Previous PDF missing: {OLD_PDF}"
        )

    if not NEW_PDF.exists():
        print(
            f"New PDF missing: {NEW_PDF}"
        )


plot_text = PLOT_SCRIPT.read_text()

print()
print("PLOTTER ROUTING")
print(
    "Uses new analysis directory: "
    f"{'B_heurippa_besckei_reanalysis' in plot_text}"
)
print(
    "Still contains old analysis directory: "
    f"{'drop_problematic_reanalysis' in plot_text}"
)

print()
print(f"Added-sample table: {ADDED_FILE}")
print(f"Removed-sample table: {REMOVED_FILE}")
print(
    f"Per-SNP BA differences: {SCORE_DIFFERENCE_FILE}"
)
