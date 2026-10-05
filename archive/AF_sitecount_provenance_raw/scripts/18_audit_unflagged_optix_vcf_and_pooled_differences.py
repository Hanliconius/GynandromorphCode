#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
GENOTYPE_FILE = (
    PROJECT / "data/derived/region_genotypes/optix.genotypes.tsv.gz"
)
METADATA_FILE = (
    PROJECT / "data/derived/region_panel_metadata.tsv"
)
CHECKPOINT_DIR = PROJECT / "analysis/checkpoints"

FOCAL = {
    "TN2L_S1", "TN2R_S2", "TN3L_S3",
    "TN3R_S4", "TN_Left", "TN_Right",
}

BAND_START = 1_112_000
BAND_END = 1_116_000


def normalize_name(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def resolve_column(dataframe, requested):
    target = normalize_name(requested)
    matches = [
        column for column in dataframe.columns
        if normalize_name(column) == target
    ]

    if len(matches) != 1:
        raise ValueError(
            f"Could not resolve {requested!r}; matches={matches}"
        )

    return matches[0]


def numeric(series):
    return pd.to_numeric(
        series.replace({"": np.nan, ".": np.nan}),
        errors="coerce",
    )


def canonical_gt(value):
    value = str(value)

    if value in {"", ".", "./.", ".|.", "nan"}:
        return None

    fields = re.split(r"[/|]", value)

    if any(field == "." for field in fields):
        return None

    try:
        alleles = sorted(int(field) for field in fields)
    except ValueError:
        return None

    return "/".join(str(allele) for allele in alleles)


def variant_type(ref, alt):
    alternate_alleles = str(alt).split(",")

    if len(alternate_alleles) > 1:
        return "multiallelic"

    if len(str(ref)) == 1 and len(str(alt)) == 1:
        return "biallelic_SNP"

    return "biallelic_indel_or_complex"


def main():
    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    accession = resolve_column(metadata, "Accession")
    taxon = resolve_column(metadata, "Taxon Name")

    metadata["sample"] = metadata[accession].astype(str)
    metadata["taxon_qc"] = metadata[taxon].astype(str)
    metadata["role_qc"] = "comparative"

    metadata.loc[
        metadata["taxon_qc"].str.contains(
            "piano keys",
            case=False,
            na=False,
        ),
        "role_qc",
    ] = "piano_keys"

    metadata.loc[
        metadata["sample"].isin(FOCAL),
        "role_qc",
    ] = "focal"

    metadata_by_sample = metadata.set_index("sample")

    calls = pd.read_csv(
        GENOTYPE_FILE,
        sep="\t",
        dtype=str,
        compression="gzip",
        keep_default_na=False,
    )

    calls["role"] = calls["sample"].map(
        metadata_by_sample["role_qc"]
    )
    calls["pos_numeric"] = numeric(calls["pos"])
    calls["DP_numeric"] = numeric(calls["DP"])
    calls["GQ_numeric"] = numeric(calls["GQ"])
    calls["canonical_GT"] = calls["GT"].map(canonical_gt)

    calls["variant_type"] = [
        variant_type(ref, alt)
        for ref, alt in zip(calls["ref"], calls["alt"])
    ]

    site_columns = [
        "chrom", "pos", "ref", "alt", "variant_type"
    ]

    site_table = calls[site_columns].drop_duplicates()

    variant_summary = (
        site_table["variant_type"]
        .value_counts()
        .rename_axis("variant_type")
        .reset_index(name="records")
    )

    band_records = site_table.loc[
        numeric(site_table["pos"]).between(
            BAND_START,
            BAND_END,
        )
    ].copy()

    pooled = calls.loc[
        calls["sample"].isin(["TN_Left", "TN_Right"])
    ].copy()

    pooled_wide = pooled.pivot_table(
        index=["chrom", "pos", "ref", "alt", "variant_type"],
        columns="sample",
        values=["canonical_GT", "DP_numeric", "GQ_numeric", "AD"],
        aggfunc="first",
    )

    pooled_wide.columns = [
        f"{field}__{sample}"
        for field, sample in pooled_wide.columns
    ]
    pooled_wide = pooled_wide.reset_index()

    for required in [
        "canonical_GT__TN_Left",
        "canonical_GT__TN_Right",
    ]:
        if required not in pooled_wide.columns:
            pooled_wide[required] = None

    pooled_wide["both_GT_called"] = (
        pooled_wide["canonical_GT__TN_Left"].notna()
        & pooled_wide["canonical_GT__TN_Right"].notna()
    )

    pooled_wide["GT_different"] = (
        pooled_wide["both_GT_called"]
        & (
            pooled_wide["canonical_GT__TN_Left"]
            != pooled_wide["canonical_GT__TN_Right"]
        )
    )

    pooled_wide["either_GT_missing"] = (
        pooled_wide["canonical_GT__TN_Left"].isna()
        | pooled_wide["canonical_GT__TN_Right"].isna()
    )

    differing = pooled_wide.loc[
        pooled_wide["GT_different"]
    ].sort_values(["chrom", "pos"])

    one_sided = pooled_wide.loc[
        pooled_wide["either_GT_missing"]
        & (
            pooled_wide["canonical_GT__TN_Left"].notna()
            | pooled_wide["canonical_GT__TN_Right"].notna()
        )
    ].sort_values(["chrom", "pos"])

    band_pooled = pooled_wide.loc[
        numeric(pooled_wide["pos"]).between(
            BAND_START,
            BAND_END,
        )
    ].sort_values(["chrom", "pos"])

    differing_file = (
        CHECKPOINT_DIR
        / "optix_pooled_left_right_GT_differences_all_records.tsv"
    )
    one_sided_file = (
        CHECKPOINT_DIR
        / "optix_pooled_one_sided_GT_calls_all_records.tsv"
    )
    band_file = (
        CHECKPOINT_DIR
        / "optix_band_1.112_1.116Mb_all_variant_records.tsv"
    )

    differing.to_csv(
        differing_file,
        sep="\t",
        index=False,
    )
    one_sided.to_csv(
        one_sided_file,
        sep="\t",
        index=False,
    )
    band_pooled.to_csv(
        band_file,
        sep="\t",
        index=False,
    )

    print("OPTIX VCF RECORD TYPES")
    print(
        variant_summary.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print("\nPOOLED LEFT/RIGHT HARD-GT AUDIT")
    print(f"Total records: {len(pooled_wide)}")
    print(
        "Both pooled GTs called: "
        f"{int(pooled_wide['both_GT_called'].sum())}"
    )
    print(
        "Different pooled GTs: "
        f"{int(pooled_wide['GT_different'].sum())}"
    )
    print(
        "Only one pooled GT called: "
        f"{len(one_sided)}"
    )

    print(
        f"\nALL RECORDS IN BAND INTERVAL "
        f"{BAND_START}-{BAND_END}: {len(band_records)}"
    )

    if not band_pooled.empty:
        print(
            band_pooled.to_csv(
                sep="\t",
                index=False,
            ).rstrip()
        )
    else:
        print("No VCF records of any type occur in this interval.")

    print("\nPOOLED SITES WITH DIFFERENT HARD GENOTYPES")

    if differing.empty:
        print("None")
    else:
        display_columns = [
            "chrom", "pos", "ref", "alt", "variant_type",
            "canonical_GT__TN_Left",
            "canonical_GT__TN_Right",
            "DP_numeric__TN_Left",
            "DP_numeric__TN_Right",
            "GQ_numeric__TN_Left",
            "GQ_numeric__TN_Right",
        ]

        print(
            differing[display_columns].to_csv(
                sep="\t",
                index=False,
            ).rstrip()
        )

    print("\nDetailed outputs:")
    print(differing_file)
    print(one_sided_file)
    print(band_file)


if __name__ == "__main__":
    main()
