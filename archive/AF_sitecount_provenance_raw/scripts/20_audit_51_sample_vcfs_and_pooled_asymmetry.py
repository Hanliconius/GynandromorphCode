#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
GENOTYPE_DIR = PROJECT / "data/derived/region_genotypes"
METADATA_FILE = PROJECT / "data/derived/region_panel_metadata.tsv"
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
    alternates = str(alt).split(",")

    if len(alternates) > 1:
        return "multiallelic"

    if len(str(ref)) == 1 and len(str(alt)) == 1:
        return "biallelic_SNP"

    return "biallelic_indel_or_complex"


def site_set(table):
    return set(
        zip(
            table["chrom"].astype(str),
            table["pos"].astype(str),
            table["ref"].astype(str),
            table["alt"].astype(str),
        )
    )


def main():
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    accession = resolve_column(metadata, "Accession")
    taxon = resolve_column(metadata, "Taxon Name")
    data_type = resolve_column(metadata, "Data Type")

    metadata["sample"] = metadata[accession].astype(str)
    metadata["taxon_qc"] = metadata[taxon].astype(str)
    metadata["data_type_qc"] = metadata[data_type].astype(str)
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

    phenotype_columns = {
        "B": resolve_column(metadata, "B"),
        "R": resolve_column(metadata, "R"),
        "D": resolve_column(metadata, "D"),
        "Yb_dorsal": resolve_column(metadata, "Yb-dorsal"),
        "Yb_ventral": resolve_column(metadata, "Yb-ventral"),
    }

    metadata_by_sample = metadata.set_index("sample")
    metadata_samples = set(metadata["sample"])

    audit_rows = []
    overlap_rows = []
    panel_tables = {}

    for panel in ["optix", "ivory"]:
        detailed_name = f"{panel}_51sample"
        detailed = pd.read_csv(
            GENOTYPE_DIR / f"{detailed_name}.genotypes.tsv.gz",
            sep="\t",
            dtype=str,
            compression="gzip",
            keep_default_na=False,
        )

        comparison = pd.read_csv(
            GENOTYPE_DIR / f"{panel}.genotypes.tsv.gz",
            sep="\t",
            dtype=str,
            compression="gzip",
            keep_default_na=False,
        )

        panel_tables[panel] = detailed

        detailed_samples = set(detailed["sample"])
        matched_samples = detailed_samples & metadata_samples
        unmatched_samples = detailed_samples - metadata_samples

        detailed_sites = detailed[
            ["chrom", "pos", "ref", "alt"]
        ].drop_duplicates()

        comparison_sites = comparison[
            ["chrom", "pos", "ref", "alt"]
        ].drop_duplicates()

        detailed_site_set = site_set(detailed_sites)
        comparison_site_set = site_set(comparison_sites)

        exact_overlap = detailed_site_set & comparison_site_set
        position_overlap = set(
            zip(
                detailed_sites["chrom"],
                detailed_sites["pos"],
            )
        ) & set(
            zip(
                comparison_sites["chrom"],
                comparison_sites["pos"],
            )
        )

        focal_present = sorted(detailed_samples & FOCAL)

        matched_metadata = metadata.loc[
            metadata["sample"].isin(matched_samples)
        ]

        role_counts = (
            matched_metadata["role_qc"]
            .value_counts()
            .to_dict()
        )

        data_type_counts = (
            matched_metadata["data_type_qc"]
            .value_counts()
            .to_dict()
        )

        audit_rows.append({
            "panel": panel,
            "samples": len(detailed_samples),
            "metadata_matches": len(matched_samples),
            "metadata_unmatched": len(unmatched_samples),
            "comparative_samples": int(
                role_counts.get("comparative", 0)
            ),
            "piano_key_samples": int(
                role_counts.get("piano_keys", 0)
            ),
            "focal_samples": int(
                role_counts.get("focal", 0)
            ),
            "CS_samples": int(data_type_counts.get("CS", 0)),
            "WGS_samples": int(data_type_counts.get("WGS", 0)),
            "variants": len(detailed_sites),
            "exact_sites_shared_with_502_VCF":
                len(exact_overlap),
            "positions_shared_with_502_VCF":
                len(position_overlap),
            "focal_names_present": ",".join(focal_present),
            "unmatched_names": ",".join(sorted(unmatched_samples)),
        })

        overlap_rows.append({
            "panel": panel,
            "detailed_51sample_sites": len(detailed_sites),
            "comparative_502sample_sites": len(comparison_sites),
            "exact_site_overlap": len(exact_overlap),
            "position_overlap": len(position_overlap),
            "detailed_only_sites": (
                len(detailed_site_set - comparison_site_set)
            ),
            "502sample_only_sites": (
                len(comparison_site_set - detailed_site_set)
            ),
        })

    audit = pd.DataFrame(audit_rows)
    overlap = pd.DataFrame(overlap_rows)

    optix = panel_tables["optix"].copy()
    optix["pos_numeric"] = numeric(optix["pos"])
    optix["DP_numeric"] = numeric(optix["DP"])
    optix["GQ_numeric"] = numeric(optix["GQ"])
    optix["canonical_GT"] = optix["GT"].map(canonical_gt)
    optix["variant_type"] = [
        variant_type(ref, alt)
        for ref, alt in zip(optix["ref"], optix["alt"])
    ]

    optix["passes_primary_qc"] = (
        optix["canonical_GT"].notna()
        & optix["DP_numeric"].ge(2)
        & optix["GQ_numeric"].ge(20)
    )

    site_columns = [
        "chrom", "pos", "ref", "alt", "variant_type"
    ]

    optix_sites = optix[site_columns].drop_duplicates()

    band_sites = optix_sites.loc[
        numeric(optix_sites["pos"]).between(
            BAND_START,
            BAND_END,
        )
    ]

    pooled = optix.loc[
        optix["sample"].isin(["TN_Left", "TN_Right"])
    ].copy()

    pooled_wide = pooled.pivot_table(
        index=site_columns,
        columns="sample",
        values=[
            "canonical_GT",
            "DP_numeric",
            "GQ_numeric",
            "AD",
            "passes_primary_qc",
        ],
        aggfunc="first",
    )

    pooled_wide.columns = [
        f"{field}__{sample}"
        for field, sample in pooled_wide.columns
    ]
    pooled_wide = pooled_wide.reset_index()
    pooled_wide["pos_numeric"] = numeric(pooled_wide["pos"])

    for column in [
        "canonical_GT__TN_Left",
        "canonical_GT__TN_Right",
        "passes_primary_qc__TN_Left",
        "passes_primary_qc__TN_Right",
    ]:
        if column not in pooled_wide:
            pooled_wide[column] = np.nan

    pooled_wide["both_called"] = (
        pooled_wide["canonical_GT__TN_Left"].notna()
        & pooled_wide["canonical_GT__TN_Right"].notna()
    )

    pooled_wide["GT_different"] = (
        pooled_wide["both_called"]
        & (
            pooled_wide["canonical_GT__TN_Left"]
            != pooled_wide["canonical_GT__TN_Right"]
        )
    )

    pooled_wide["both_pass_primary_qc"] = (
        pooled_wide["passes_primary_qc__TN_Left"].eq(True)
        & pooled_wide["passes_primary_qc__TN_Right"].eq(True)
    )

    pooled_wide["GT_different_primary_qc"] = (
        pooled_wide["GT_different"]
        & pooled_wide["both_pass_primary_qc"]
    )

    pooled_wide["left_homo_right_heterozygous"] = (
        pooled_wide["canonical_GT__TN_Left"].isin(
            ["0/0", "1/1"]
        )
        & pooled_wide["canonical_GT__TN_Right"].eq("0/1")
    )

    pooled_wide["left_alt_homo_right_heterozygous"] = (
        pooled_wide["canonical_GT__TN_Left"].eq("1/1")
        & pooled_wide["canonical_GT__TN_Right"].eq("0/1")
    )

    pooled_differences = pooled_wide.loc[
        pooled_wide["GT_different"]
    ].copy()

    band_pooled = pooled_wide.loc[
        pooled_wide["pos_numeric"].between(
            BAND_START,
            BAND_END,
        )
    ].copy()

    band_differences = band_pooled.loc[
        band_pooled["GT_different"]
    ].copy()

    audit_file = (
        CHECKPOINT_DIR / "region_51sample_panel_audit.tsv"
    )
    overlap_file = (
        CHECKPOINT_DIR / "region_51sample_vs_502_site_overlap.tsv"
    )
    pooled_file = (
        CHECKPOINT_DIR
        / "optix_51sample_pooled_left_right_differences.tsv"
    )
    band_file = (
        CHECKPOINT_DIR
        / "optix_51sample_band_interval_pooled_genotypes.tsv"
    )

    audit.to_csv(audit_file, sep="\t", index=False)
    overlap.to_csv(overlap_file, sep="\t", index=False)
    pooled_differences.to_csv(
        pooled_file,
        sep="\t",
        index=False,
    )
    band_pooled.to_csv(
        band_file,
        sep="\t",
        index=False,
    )

    print("51-SAMPLE PANEL COMPOSITION")
    print(audit.to_csv(sep="\t", index=False).rstrip())

    print("\nSITE OVERLAP WITH 502-SAMPLE FILTERED VCF")
    print(overlap.to_csv(sep="\t", index=False).rstrip())

    print("\nOPTIX 51-SAMPLE POOLED ASYMMETRY")
    print(f"All optix records: {len(pooled_wide)}")
    print(f"Band-interval records: {len(band_sites)}")
    print(
        "Different pooled hard GTs, raw: "
        f"{int(pooled_wide['GT_different'].sum())}"
    )
    print(
        "Different pooled hard GTs, both DP>=2/GQ>=20: "
        f"{int(pooled_wide['GT_different_primary_qc'].sum())}"
    )
    print(
        "Left homozygous/right heterozygous, raw: "
        f"{int(pooled_wide['left_homo_right_heterozygous'].sum())}"
    )
    print(
        "Left ALT-homozygous/right heterozygous, raw: "
        f"{int(pooled_wide['left_alt_homo_right_heterozygous'].sum())}"
    )
    print(
        "Band-interval pooled differences, raw: "
        f"{len(band_differences)}"
    )

    print("\nBAND-INTERVAL POOLED DIFFERENCES")

    if band_differences.empty:
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
            "both_pass_primary_qc",
        ]

        print(
            band_differences[display_columns]
            .head(40)
            .to_csv(
                sep="\t",
                index=False,
            )
            .rstrip()
        )

        if len(band_differences) > 40:
            print(
                f"... {len(band_differences) - 40} additional "
                "Band-interval differences saved to file."
            )

    print("\nDetailed outputs:")
    print(audit_file)
    print(overlap_file)
    print(pooled_file)
    print(band_file)


if __name__ == "__main__":
    main()
