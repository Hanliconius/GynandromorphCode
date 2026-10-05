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

FOCAL_TABLE = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "core_explanatory_optix_genotypes.long.tsv"
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
    / "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

OUTDIR = (
    PROJECT
    / "analysis/original_502/checkpoints"
)

OUT_LONG = (
    OUTDIR
    / "optix_D_R_expected_focal_site_genotypes.long.tsv"
)

OUT_SITES = (
    OUTDIR
    / "optix_D_R_expected_focal_sites.tsv"
)

OUT_SAMPLES = (
    OUTDIR
    / "optix_D_R_expected_site_genotype_samples.tsv"
)

OUT_SUMMARY = (
    OUTDIR
    / "optix_D_R_expected_site_genotypes.summary.tsv"
)


def as_bool(value):
    return str(value).strip().lower() in {
        "true", "t", "1", "yes", "y"
    }


def clean_text(value):
    return re.sub(r"\s+", " ", str(value)).strip()


def normalize_module(value):
    value = clean_text(value).lower()

    if value.startswith("dennis"):
        return "Dennis"

    if value.startswith("ray"):
        return "Ray"

    return clean_text(value)


def number(value):
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
    FOCAL_TABLE,
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

focal = pd.read_csv(
    FOCAL_TABLE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

manifest = pd.read_csv(
    COHORT_MANIFEST,
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

metadata.columns = [
    clean_text(column)
    for column in metadata.columns
]

required_strict = {
    "module", "chrom", "pos", "ref", "alt",
    "minimum_BA_all_analyses",
}

missing = required_strict.difference(strict.columns)

if missing:
    raise RuntimeError(
        "Strict-site table lacks columns: "
        + ", ".join(sorted(missing))
    )

required_focal = {
    "module", "chrom", "pos", "ref", "alt",
    "sample", "canonical_GT", "biological_GT",
    "DP2_GQ20_pass",
}

missing = required_focal.difference(focal.columns)

if missing:
    raise RuntimeError(
        "Focal table lacks columns: "
        + ", ".join(sorted(missing))
    )

strict["module"] = strict["module"].map(normalize_module)
focal["module"] = focal["module"].map(normalize_module)

strict = strict.loc[
    strict["module"].isin(["Dennis", "Ray"])
].copy()

key_columns = [
    "module", "chrom", "pos", "ref", "alt"
]

strict = strict.drop_duplicates(
    subset=key_columns
)

focal_lr = focal.loc[
    focal["sample"].isin(["TN_Left", "TN_Right"])
].copy()

focal_lr = focal_lr.loc[
    focal_lr["module"].isin(["Dennis", "Ray"])
].copy()

left = focal_lr.loc[
    focal_lr["sample"] == "TN_Left"
].copy()

right = focal_lr.loc[
    focal_lr["sample"] == "TN_Right"
].copy()

left = left.rename(
    columns={
        "canonical_GT": "left_canonical_GT",
        "biological_GT": "left_biological_GT",
        "DP2_GQ20_pass": "left_QC",
    }
)

right = right.rename(
    columns={
        "canonical_GT": "right_canonical_GT",
        "biological_GT": "right_biological_GT",
        "DP2_GQ20_pass": "right_QC",
    }
)

site_table = strict.merge(
    left[
        key_columns
        + [
            "left_canonical_GT",
            "left_biological_GT",
            "left_QC",
        ]
    ],
    on=key_columns,
    how="left",
    validate="one_to_one",
)

site_table = site_table.merge(
    right[
        key_columns
        + [
            "right_canonical_GT",
            "right_biological_GT",
            "right_QC",
        ]
    ],
    on=key_columns,
    how="left",
    validate="one_to_one",
)

site_table["left_QC_bool"] = site_table["left_QC"].map(as_bool)
site_table["right_QC_bool"] = site_table["right_QC"].map(as_bool)

site_table["exact_expected_focal_pattern"] = (
    site_table["left_QC_bool"]
    & site_table["right_QC_bool"]
    & (site_table["left_biological_GT"] == "P/P")
    & (site_table["right_biological_GT"] == "P/A")
)

sites = site_table.loc[
    site_table["exact_expected_focal_pattern"]
].copy()

if sites.empty:
    raise RuntimeError(
        "No exact expected focal sites were recovered."
    )

pattern_indices = []

for row in sites.itertuples(index=False):
    gt = parse_gt(row.left_canonical_GT)

    if gt == (0, 0):
        pattern_indices.append(0)
    elif gt == (1, 1):
        pattern_indices.append(1)
    else:
        raise RuntimeError(
            "Expected a homozygous left genotype at "
            f"{row.module} {row.chrom}:{row.pos}; "
            f"found {row.left_canonical_GT}"
        )

sites["pattern_allele_index"] = pattern_indices
sites["pos"] = sites["pos"].astype(int)
sites["minimum_BA_all_analyses"] = pd.to_numeric(
    sites["minimum_BA_all_analyses"],
    errors="coerce",
)

sites = sites.sort_values(
    ["module", "pos", "ref", "alt"]
).reset_index(drop=True)

sites.to_csv(
    OUT_SITES,
    sep="\t",
    index=False,
)

required_manifest = {
    "cohort", "module", "sample",
    "taxon_name", "data_type", "keep",
}

missing = required_manifest.difference(manifest.columns)

if missing:
    raise RuntimeError(
        "Cohort manifest lacks columns: "
        + ", ".join(sorted(missing))
    )

manifest["normalized_module"] = manifest["module"].map(
    normalize_module
)

manifest["keep_bool"] = manifest["keep"].map(as_bool)

cohort = manifest.loc[
    (manifest["cohort"] == "WGS_plus_CS")
    & manifest["keep_bool"]
    & manifest["normalized_module"].isin(["Dennis", "Ray"])
].copy()

cohort["module"] = cohort["normalized_module"]

cohort = cohort.drop_duplicates(
    subset=["module", "sample"]
)

accession_column = next(
    (
        column for column in metadata.columns
        if column.lower() == "accession"
    ),
    None,
)

if accession_column is None:
    raise RuntimeError(
        "Could not identify the metadata accession column."
    )

metadata[accession_column] = metadata[
    accession_column
].map(clean_text)

metadata_lookup = (
    metadata.drop_duplicates(subset=[accession_column])
    .set_index(accession_column)
)

sample_rows = []

for row in cohort.itertuples(index=False):
    trait = "D" if row.module == "Dennis" else "R"

    phenotype = ""

    if (
        row.sample in metadata_lookup.index
        and trait in metadata_lookup.columns
    ):
        phenotype = clean_text(
            metadata_lookup.loc[row.sample, trait]
        )

    if phenotype in {"0.0", "0"}:
        phenotype = "0"
    elif phenotype in {"1.0", "1"}:
        phenotype = "1"
    else:
        phenotype = "NA"

    sample_rows.append(
        {
            "module": row.module,
            "sample": row.sample,
            "sample_type": "comparative",
            "trait": trait,
            "phenotype": phenotype,
            "taxon_name": clean_text(row.taxon_name),
            "data_type": clean_text(row.data_type),
        }
    )

for module in ["Dennis", "Ray"]:
    trait = "D" if module == "Dennis" else "R"

    for sample, side in [
        ("TN_Left", "Left"),
        ("TN_Right", "Right"),
    ]:
        sample_rows.append(
            {
                "module": module,
                "sample": sample,
                "sample_type": "focal",
                "trait": trait,
                "phenotype": side,
                "taxon_name": "gynandromorph pooled sample",
                "data_type": "WGS",
            }
        )

samples = pd.DataFrame(sample_rows)

samples = samples.drop_duplicates(
    subset=["module", "sample"]
)

samples.to_csv(
    OUT_SAMPLES,
    sep="\t",
    index=False,
)

target_modules = {}

for row in sites.itertuples(index=False):
    key = (
        row.chrom,
        int(row.pos),
        row.ref,
        row.alt,
    )

    target_modules.setdefault(key, []).append(
        {
            "module": row.module,
            "pattern_allele_index": int(
                row.pattern_allele_index
            ),
            "comparative_BA": row.minimum_BA_all_analyses,
        }
    )

wanted_samples = set(samples["sample"])

records_found = set()
output_rows = []

with VCF.open("r", encoding="utf-8", errors="replace") as handle:
    vcf_samples = None

    for line_number, line in enumerate(handle, start=1):
        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            fields = line.rstrip("\n").split("\t")
            vcf_samples = fields[9:]
            sample_index = {
                sample: index
                for index, sample in enumerate(vcf_samples)
            }

            missing_samples = wanted_samples.difference(sample_index)

            if missing_samples:
                raise RuntimeError(
                    "Requested samples absent from VCF: "
                    + ",".join(sorted(missing_samples))
                )

            continue

        if vcf_samples is None:
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 10:
            continue

        chrom = fields[0]

        try:
            pos = int(fields[1])
        except ValueError:
            continue

        ref = fields[3]
        alt = fields[4]

        key = (chrom, pos, ref, alt)

        if key not in target_modules:
            continue

        records_found.add(key)

        format_fields = fields[8].split(":")
        format_index = {
            name: index
            for index, name in enumerate(format_fields)
        }

        for site_info in target_modules[key]:
            module = site_info["module"]
            pattern_index = site_info["pattern_allele_index"]

            module_samples = samples.loc[
                samples["module"] == module
            ]

            for sample_row in module_samples.itertuples(index=False):
                sample_field = fields[
                    9 + sample_index[sample_row.sample]
                ]

                values = sample_field.split(":")

                def field_value(name):
                    index = format_index.get(name)

                    if index is None or index >= len(values):
                        return ""

                    return values[index]

                raw_gt = field_value("GT")
                ad = field_value("AD")
                dp = number(field_value("DP"))
                gq = number(field_value("GQ"))
                parsed_gt = parse_gt(raw_gt)

                qc_pass = (
                    parsed_gt is not None
                    and dp is not None
                    and gq is not None
                    and dp >= 2
                    and gq >= 20
                )

                biological_gt = "N/QC"
                pattern_dosage = None

                if parsed_gt is not None:
                    pattern_count = sum(
                        allele == pattern_index
                        for allele in parsed_gt
                    )

                    pattern_dosage = pattern_count / 2

                    biological_gt = {
                        0: "A/A",
                        1: "P/A",
                        2: "P/P",
                    }[pattern_count]

                if not qc_pass:
                    display_gt = "N/QC"
                else:
                    display_gt = biological_gt

                if sample_row.sample_type == "focal":
                    row_group = "Focal gynandromorph"
                    sample_label = (
                        f"{sample_row.phenotype} | "
                        f"{sample_row.sample}"
                    )
                else:
                    state = (
                        "present"
                        if sample_row.phenotype == "1"
                        else "absent"
                        if sample_row.phenotype == "0"
                        else "phenotype NA"
                    )

                    row_group = (
                        f"{module} {state}"
                    )

                    sample_label = (
                        f"{sample_row.sample} | "
                        f"{sample_row.taxon_name} | "
                        f"{sample_row.trait}="
                        f"{sample_row.phenotype} | "
                        f"{sample_row.data_type}"
                    )

                output_rows.append(
                    {
                        "module": module,
                        "trait": sample_row.trait,
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "site_label": f"{pos:,}",
                        "comparative_BA": (
                            site_info["comparative_BA"]
                        ),
                        "pattern_allele_index": pattern_index,
                        "sample": sample_row.sample,
                        "sample_type": sample_row.sample_type,
                        "phenotype": sample_row.phenotype,
                        "taxon_name": sample_row.taxon_name,
                        "data_type": sample_row.data_type,
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

missing_records = set(target_modules).difference(records_found)

if missing_records:
    preview = sorted(missing_records)[:20]

    raise RuntimeError(
        "Required SNP records absent from VCF. First entries: "
        + "; ".join(
            f"{chrom}:{pos}:{ref}:{alt}"
            for chrom, pos, ref, alt in preview
        )
    )

long_table = pd.DataFrame(output_rows)

long_table.to_csv(
    OUT_LONG,
    sep="\t",
    index=False,
)

summary_rows = []

for module in ["Dennis", "Ray"]:
    d = long_table.loc[
        long_table["module"] == module
    ].copy()

    comparative = d.loc[
        d["sample_type"] == "comparative"
    ]

    focal = d.loc[
        d["sample_type"] == "focal"
    ]

    summary_rows.append(
        {
            "module": module,
            "sites": d[
                ["chrom", "pos", "ref", "alt"]
            ].drop_duplicates().shape[0],
            "comparative_samples": comparative[
                "sample"
            ].nunique(),
            "focal_samples": focal["sample"].nunique(),
            "comparative_cells": len(comparative),
            "comparative_QC_cells": int(
                comparative["QC_pass"].sum()
            ),
            "comparative_call_rate": (
                comparative["QC_pass"].mean()
            ),
            "focal_cells": len(focal),
            "focal_QC_cells": int(
                focal["QC_pass"].sum()
            ),
        }
    )

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False,
)

print("D/R EXPECTED-SITE GENOTYPE EXTRACTION")
print()
print(summary.to_string(index=False))
print()
print(f"Long table: {OUT_LONG}")
print(f"Site table: {OUT_SITES}")
print(f"Sample table: {OUT_SAMPLES}")
print(f"Summary: {OUT_SUMMARY}")
