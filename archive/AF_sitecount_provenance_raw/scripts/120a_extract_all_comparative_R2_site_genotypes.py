#!/usr/bin/env python3

from pathlib import Path
import re

import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")


STRICT_SITES = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "broad_optix_strict_comparative_SNPs.tsv"
)

COHORT_MANIFEST = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)

METADATA = (
    PROJECT
    / "data/derived/region_panel_metadata.tsv"
)

VCF = (
    PROJECT
    / "external/krista_genotype_plots_original/Genotype_Plots"
    / "ray2-1174547-1187833.vcf"
)

OUTDIR = (
    PROJECT
    / "analysis/original_502/checkpoints"
)

OUT_LONG = (
    OUTDIR
    / "optix_R2_all_comparative_selected_genotypes.long.tsv"
)

OUT_SITES = (
    OUTDIR
    / "optix_R2_all_comparative_selected_sites.tsv"
)

OUT_SUMMARY = (
    OUTDIR
    / "optix_R2_all_comparative_selected_genotypes.summary.tsv"
)

R2_START = 1174547
R2_END = 1187833


def clean(value):
    return re.sub(r"\s+", " ", str(value)).strip()


def as_bool(value):
    return clean(value).lower() in {
        "true", "t", "1", "yes", "y"
    }


def normalize_module(value):
    value = clean(value).lower()

    if value.startswith("ray"):
        return "Ray"

    if value.startswith("dennis"):
        return "Dennis"

    if value.startswith("band"):
        return "Band"

    return clean(value)


def numeric(value):
    try:
        if value in {"", ".", None}:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_gt(gt):
    if gt in {"", ".", "./.", ".|."}:
        return None

    alleles = re.split(r"[/|]", gt)

    if len(alleles) != 2:
        return None

    if any(allele not in {"0", "1"} for allele in alleles):
        return None

    return tuple(int(allele) for allele in alleles)


for required in [
    STRICT_SITES,
    COHORT_MANIFEST,
    METADATA,
    VCF,
]:
    if not required.is_file():
        raise RuntimeError(f"Missing required input: {required}")

OUTDIR.mkdir(parents=True, exist_ok=True)

strict = pd.read_csv(
    STRICT_SITES,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

strict["module"] = strict["module"].map(normalize_module)
strict["pos"] = pd.to_numeric(strict["pos"], errors="raise").astype(int)

sites = strict.loc[
    (strict["module"] == "Ray")
    & (strict["chrom"] == "Hmr18")
    & strict["pos"].between(R2_START, R2_END)
].copy()

sites = sites.drop_duplicates(
    subset=["chrom", "pos", "ref", "alt"]
)

sites = sites.sort_values(
    ["pos", "ref", "alt"]
).reset_index(drop=True)

if sites.empty:
    raise RuntimeError(
        "No comparative-selected Ray sites were found in R2."
    )

# For Ray, the reference genome carries the non-pattern state.
# Therefore ALT is the Ray-pattern-associated allele.
sites["pattern_allele"] = "ALT"
sites["pattern_allele_index"] = 1

sites.to_csv(
    OUT_SITES,
    sep="\t",
    index=False,
)

manifest = pd.read_csv(
    COHORT_MANIFEST,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

manifest["normalized_module"] = manifest["module"].map(
    normalize_module
)

manifest["keep_bool"] = manifest["keep"].map(as_bool)

cohort = manifest.loc[
    (manifest["cohort"] == "WGS_plus_CS")
    & (manifest["normalized_module"] == "Ray")
    & manifest["keep_bool"]
].copy()

cohort = cohort.drop_duplicates(subset=["sample"])

metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

metadata.columns = [clean(column) for column in metadata.columns]

accession_column = next(
    column
    for column in metadata.columns
    if column.lower() == "accession"
)

metadata[accession_column] = metadata[accession_column].map(clean)

metadata = metadata.drop_duplicates(
    subset=[accession_column]
).set_index(accession_column)

sample_information = {}

for row in cohort.itertuples(index=False):
    phenotype = "NA"

    if row.sample in metadata.index and "R" in metadata.columns:
        phenotype = clean(metadata.loc[row.sample, "R"])

    if phenotype in {"0", "0.0"}:
        phenotype = "0"
    elif phenotype in {"1", "1.0"}:
        phenotype = "1"
    else:
        phenotype = "NA"

    sample_information[row.sample] = {
        "sample_type": "comparative",
        "phenotype": phenotype,
        "taxon_name": clean(row.taxon_name),
        "data_type": clean(row.data_type),
    }

sample_information["TN_Left"] = {
    "sample_type": "focal",
    "phenotype": "Left",
    "taxon_name": "gynandromorph pooled left",
    "data_type": "WGS",
}

sample_information["TN_Right"] = {
    "sample_type": "focal",
    "phenotype": "Right",
    "taxon_name": "gynandromorph pooled right",
    "data_type": "WGS",
}

target_sites = {
    (
        row.chrom,
        int(row.pos),
        row.ref,
        row.alt,
    ): row
    for row in sites.itertuples(index=False)
}

output_rows = []
found_sites = set()

with VCF.open("r", encoding="utf-8", errors="replace") as handle:
    sample_index = None

    for line in handle:
        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            fields = line.rstrip("\n").split("\t")
            vcf_samples = fields[9:]

            sample_index = {
                sample: index
                for index, sample in enumerate(vcf_samples)
            }

            missing_samples = set(sample_information).difference(
                sample_index
            )

            if missing_samples:
                raise RuntimeError(
                    "Samples missing from R2 VCF: "
                    + ",".join(sorted(missing_samples))
                )

            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 10:
            continue

        key = (
            fields[0],
            int(fields[1]),
            fields[3],
            fields[4],
        )

        if key not in target_sites:
            continue

        found_sites.add(key)
        site = target_sites[key]

        format_names = fields[8].split(":")
        format_index = {
            name: index
            for index, name in enumerate(format_names)
        }

        for sample, sample_info in sample_information.items():
            sample_values = fields[
                9 + sample_index[sample]
            ].split(":")

            def get_field(name):
                index = format_index.get(name)

                if index is None or index >= len(sample_values):
                    return ""

                return sample_values[index]

            raw_gt = get_field("GT")
            ad = get_field("AD")
            dp = numeric(get_field("DP"))
            gq = numeric(get_field("GQ"))
            gt = parse_gt(raw_gt)

            qc_pass = (
                gt is not None
                and dp is not None
                and gq is not None
                and dp >= 2
                and gq >= 20
            )

            if gt is None:
                pattern_dosage = None
                biological_gt = "N/QC"
            else:
                pattern_count = sum(allele == 1 for allele in gt)
                pattern_dosage = pattern_count / 2

                biological_gt = {
                    0: "A/A",
                    1: "P/A",
                    2: "P/P",
                }[pattern_count]

            display_gt = biological_gt if qc_pass else "N/QC"

            if sample_info["sample_type"] == "focal":
                row_group = "Focal gynandromorph"
                sample_label = (
                    f'{sample_info["phenotype"]} | {sample}'
                )
            else:
                if sample_info["phenotype"] == "1":
                    state = "Ray present"
                elif sample_info["phenotype"] == "0":
                    state = "Ray absent"
                else:
                    state = "Phenotype unavailable"

                row_group = state

                sample_label = (
                    f'{sample} | '
                    f'{sample_info["taxon_name"]} | '
                    f'R={sample_info["phenotype"]} | '
                    f'{sample_info["data_type"]}'
                )

            output_rows.append(
                {
                    "module": "Ray",
                    "region": "R2",
                    "chrom": key[0],
                    "pos": key[1],
                    "ref": key[2],
                    "alt": key[3],
                    "site_label": f"{key[1]:,}",
                    "minimum_BA_all_analyses": getattr(
                        site,
                        "minimum_BA_all_analyses",
                        ""
                    ),
                    "sample": sample,
                    "sample_type": sample_info["sample_type"],
                    "phenotype": sample_info["phenotype"],
                    "taxon_name": sample_info["taxon_name"],
                    "data_type": sample_info["data_type"],
                    "row_group": row_group,
                    "sample_label": sample_label,
                    "raw_GT": raw_gt,
                    "AD": ad,
                    "DP": dp,
                    "GQ": gq,
                    "QC_pass": qc_pass,
                    "pattern_dosage": pattern_dosage,
                    "biological_GT": biological_gt,
                    "display_GT": display_gt,
                }
            )

missing_sites = set(target_sites).difference(found_sites)

if missing_sites:
    preview = sorted(missing_sites)[:20]

    raise RuntimeError(
        "Comparative-selected sites missing from regional VCF: "
        + ";".join(
            f"{chrom}:{pos}:{ref}:{alt}"
            for chrom, pos, ref, alt in preview
        )
    )

long_table = pd.DataFrame(output_rows)

# Classify focal L/R outcome without filtering any site.
focal = long_table.loc[
    long_table["sample_type"] == "focal"
].pivot_table(
    index=["chrom", "pos", "ref", "alt"],
    columns="phenotype",
    values=["display_GT", "QC_pass"],
    aggfunc="first",
)

focal.columns = [
    f"{field}_{side}"
    for field, side in focal.columns
]

focal = focal.reset_index()

focal["focal_outcome"] = "one_or_both_fail_QC"

both_qc = (
    focal["QC_pass_Left"].map(as_bool)
    & focal["QC_pass_Right"].map(as_bool)
)

identical = (
    both_qc
    & (
        focal["display_GT_Left"]
        == focal["display_GT_Right"]
    )
)

expected = (
    both_qc
    & (focal["display_GT_Left"] == "P/P")
    & (focal["display_GT_Right"] == "P/A")
)

other_difference = (
    both_qc
    & ~identical
    & ~expected
)

focal.loc[identical, "focal_outcome"] = "identical_left_right"
focal.loc[expected, "focal_outcome"] = "expected_P_P_to_P_A"
focal.loc[other_difference, "focal_outcome"] = "other_L_R_difference"

long_table = long_table.merge(
    focal[
        [
            "chrom",
            "pos",
            "ref",
            "alt",
            "focal_outcome",
        ]
    ],
    on=["chrom", "pos", "ref", "alt"],
    how="left",
    validate="many_to_one",
)

long_table.to_csv(
    OUT_LONG,
    sep="\t",
    index=False,
)

outcome_counts = (
    focal["focal_outcome"]
    .value_counts()
    .to_dict()
)

summary = pd.DataFrame(
    [
        {
            "region": "R2",
            "comparative_selected_SNPs": len(sites),
            "comparative_samples": len(cohort),
            "expected_P_P_to_P_A": outcome_counts.get(
                "expected_P_P_to_P_A", 0
            ),
            "identical_left_right": outcome_counts.get(
                "identical_left_right", 0
            ),
            "other_L_R_difference": outcome_counts.get(
                "other_L_R_difference", 0
            ),
            "one_or_both_fail_QC": outcome_counts.get(
                "one_or_both_fail_QC", 0
            ),
        }
    ]
)

summary.to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False,
)

print("ALL COMPARATIVE-SELECTED R2 GENOTYPES")
print()
print(summary.to_string(index=False))
print()
print(f"Long table: {OUT_LONG}")
print(f"Sites: {OUT_SITES}")
print(f"Summary: {OUT_SUMMARY}")
