#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

VCF = (
    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
    / "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

SELECTED_SITES = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
    / "broad_optix_strict_comparative_SNPs.tsv"
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
    / "original_502"
    / "checkpoints"
)

OUT_SAMPLES = (
    OUTDIR
    / "core_explanatory_optix_samples.tsv"
)

OUT_SITES = (
    OUTDIR
    / "core_explanatory_optix_sites.tsv"
)

OUT_LONG = (
    OUTDIR
    / "core_explanatory_optix_genotypes.long.tsv"
)

OUT_MISSING = (
    OUTDIR
    / "core_explanatory_optix_missing_sites.tsv"
)


def clean_text(value):

    if pd.isna(value):
        return ""

    return re.sub(
        r"\s+",
        " ",
        str(value),
    ).strip()


def phenotype_text(value):

    value = pd.to_numeric(
        pd.Series([value]),
        errors="coerce",
    ).iloc[0]

    if pd.isna(value):
        return "NA"

    return str(
        int(value)
    )


def numeric_or_nan(value):

    try:
        if value in {
            "",
            ".",
            None,
        }:
            return np.nan

        return float(value)

    except (TypeError, ValueError):
        return np.nan


def canonical_gt(value):

    if value in {
        "",
        ".",
        "./.",
        ".|.",
        None,
    }:
        return None

    alleles = re.split(
        r"[/|]",
        value,
    )

    if (
        len(alleles) != 2
        or
        any(
            allele == "."
            for allele in alleles
        )
    ):
        return None

    try:
        numeric = sorted(
            int(allele)
            for allele in alleles
        )

    except ValueError:
        return None

    if any(
        allele not in {0, 1}
        for allele in numeric
    ):
        return None

    return "/".join(
        str(allele)
        for allele in numeric
    )


def parse_ad(value):

    if value in {
        "",
        ".",
        None,
    }:
        return (
            np.nan,
            np.nan,
        )

    fields = value.split(",")

    if len(fields) < 2:
        return (
            np.nan,
            np.nan,
        )

    return (
        numeric_or_nan(fields[0]),
        numeric_or_nan(fields[1]),
    )


def pattern_dosage(gt, module):

    if gt is None:
        return np.nan

    alleles = [
        int(value)
        for value in gt.split("/")
    ]

    # Reference genome carries Band but lacks Ray and Dennis.
    pattern_allele = (
        0
        if module == "Band"
        else 1
    )

    return sum(
        allele == pattern_allele
        for allele in alleles
    ) / 2


def biological_gt_label(dosage):

    if not np.isfinite(dosage):
        return "missing"

    if dosage == 0:
        return "A/A"

    if dosage == 0.5:
        return "P/A"

    if dosage == 1:
        return "P/P"

    return "other"


OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)

sites = pd.read_csv(
    SELECTED_SITES,
    sep="\t",
)

sites["pos"] = pd.to_numeric(
    sites["pos"],
    errors="raise",
).astype(int)

sites = sites.sort_values(
    [
        "module",
        "pos",
        "ref",
        "alt",
    ]
).reset_index(
    drop=True
)

metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

required_metadata = {
    "Accession",
    "Taxon Name",
    "Data Type",
    "B",
    "R",
    "D",
}

missing_metadata = (
    required_metadata
    -
    set(metadata.columns)
)

if missing_metadata:
    raise RuntimeError(
        "Metadata is missing columns: "
        + ", ".join(
            sorted(
                missing_metadata
            )
        )
    )

metadata["taxon_clean"] = (
    metadata[
        "Taxon Name"
    ]
    .map(
        clean_text
    )
)

metadata["taxon_lower"] = (
    metadata[
        "taxon_clean"
    ]
    .str.lower()
)

focal_samples = {
    "TN_Left":
        "Gynandromorph left — expected B0 R1 D1",

    "TN_Right":
        "Gynandromorph right — expected B1 R1 D1",
}

sample_rows = []

for _, row in metadata.iterrows():

    sample = clean_text(
        row[
            "Accession"
        ]
    )

    taxon = row[
        "taxon_clean"
    ]

    taxon_lower = row[
        "taxon_lower"
    ]

    group = None
    biological_note = ""

    if sample in focal_samples:

        group = "Focal gynandromorph"
        biological_note = focal_samples[
            sample
        ]

    elif "besckei" in taxon_lower:

        group = "H. besckei"
        biological_note = (
            "Band-associated introgression candidate"
        )

    elif "contigua" in taxon_lower:

        group = "H. t. contigua"
        biological_note = (
            "Ray present; Dennis absent"
        )

    elif "meriana" in taxon_lower:

        group = "H. m. meriana"
        biological_note = (
            "Dennis present; Ray absent"
        )

    if group is None:
        continue

    B = phenotype_text(
        row[
            "B"
        ]
    )

    R = phenotype_text(
        row[
            "R"
        ]
    )

    D = phenotype_text(
        row[
            "D"
        ]
    )

    data_type = clean_text(
        row[
            "Data Type"
        ]
    )

    if sample in focal_samples:

        phenotype_label = (
            "focal expectation"
        )

    else:

        phenotype_label = (
            f"B={B} R={R} D={D}"
        )

    sample_label = (
        f"{sample} | "
        f"{taxon} | "
        f"{phenotype_label} | "
        f"{data_type}"
    )

    sample_rows.append(
        {
            "sample":
                sample,

            "group":
                group,

            "taxon":
                taxon,

            "data_type":
                data_type,

            "B":
                B,

            "R":
                R,

            "D":
                D,

            "biological_note":
                biological_note,

            "sample_label":
                sample_label,
        }
    )

samples = pd.DataFrame(
    sample_rows
).drop_duplicates(
    subset=[
        "sample"
    ]
)

group_order = {
    "Focal gynandromorph": 0,
    "H. besckei": 1,
    "H. t. contigua": 2,
    "H. m. meriana": 3,
}

samples["group_order"] = (
    samples[
        "group"
    ]
    .map(
        group_order
    )
)

samples = samples.sort_values(
    [
        "group_order",
        "taxon",
        "sample",
    ]
).reset_index(
    drop=True
)

samples["plot_order"] = np.arange(
    1,
    len(samples) + 1,
)

site_lookup = {}

for _, row in sites.iterrows():

    key = (
        str(
            row[
                "chrom"
            ]
        ),
        int(
            row[
                "pos"
            ]
        ),
        str(
            row[
                "ref"
            ]
        ),
        str(
            row[
                "alt"
            ]
        ),
    )

    site_lookup.setdefault(
        key,
        [],
    ).append(
        row.to_dict()
    )

requested_samples = set(
    samples[
        "sample"
    ]
)

long_rows = []
found_site_keys = set()
vcf_samples = []
sample_indices = {}

with VCF.open() as handle:

    for line in handle:

        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):

            header = line.rstrip(
                "\n"
            ).split(
                "\t"
            )

            vcf_samples = header[9:]

            sample_indices = {
                sample:
                    vcf_samples.index(
                        sample
                    )
                for sample in requested_samples
                if sample in vcf_samples
            }

            continue

        fields = line.rstrip(
            "\n"
        ).split(
            "\t"
        )

        if len(fields) < 10:
            continue

        key = (
            fields[0],
            int(
                fields[1]
            ),
            fields[3],
            fields[4],
        )

        selected_records = site_lookup.get(
            key
        )

        if selected_records is None:
            continue

        found_site_keys.add(
            key
        )

        format_fields = fields[
            8
        ].split(
            ":"
        )

        format_index = {
            name: index
            for index, name
            in enumerate(
                format_fields
            )
        }

        for selected_record in selected_records:

            module = selected_record[
                "module"
            ]

            for _, sample_row in samples.iterrows():

                sample = sample_row[
                    "sample"
                ]

                if sample not in sample_indices:
                    continue

                sample_field = fields[
                    9
                    +
                    sample_indices[
                        sample
                    ]
                ].split(
                    ":"
                )

                def get_format(name):

                    index = format_index.get(
                        name
                    )

                    if (
                        index is None
                        or
                        index >= len(
                            sample_field
                        )
                    ):
                        return "."

                    return sample_field[
                        index
                    ]

                raw_gt = get_format(
                    "GT"
                )

                gt = canonical_gt(
                    raw_gt
                )

                ad_value = get_format(
                    "AD"
                )

                ref_depth, alt_depth = parse_ad(
                    ad_value
                )

                dp = numeric_or_nan(
                    get_format(
                        "DP"
                    )
                )

                gq = numeric_or_nan(
                    get_format(
                        "GQ"
                    )
                )

                dosage = pattern_dosage(
                    gt,
                    module,
                )

                qc_pass = (
                    gt is not None
                    and
                    np.isfinite(dp)
                    and
                    np.isfinite(gq)
                    and
                    dp >= 2
                    and
                    gq >= 20
                )

                long_rows.append(
                    {
                        "module":
                            module,

                        "chrom":
                            key[0],

                        "pos":
                            key[1],

                        "ref":
                            key[2],

                        "alt":
                            key[3],

                        "minimum_BA_all_analyses":
                            selected_record[
                                "minimum_BA_all_analyses"
                            ],

                        "sample":
                            sample,

                        "group":
                            sample_row[
                                "group"
                            ],

                        "plot_order":
                            sample_row[
                                "plot_order"
                            ],

                        "sample_label":
                            sample_row[
                                "sample_label"
                            ],

                        "B":
                            sample_row[
                                "B"
                            ],

                        "R":
                            sample_row[
                                "R"
                            ],

                        "D":
                            sample_row[
                                "D"
                            ],

                        "raw_GT":
                            raw_gt,

                        "canonical_GT":
                            (
                                gt
                                if gt is not None
                                else "missing"
                            ),

                        "AD":
                            ad_value,

                        "ref_depth":
                            ref_depth,

                        "alt_depth":
                            alt_depth,

                        "DP":
                            dp,

                        "GQ":
                            gq,

                        "pattern_dosage":
                            dosage,

                        "biological_GT":
                            biological_gt_label(
                                dosage
                            ),

                        "DP2_GQ20_pass":
                            qc_pass,
                    }
                )

missing_sample_names = sorted(
    requested_samples
    -
    set(
        sample_indices
    )
)

missing_sites = []

for key, records in site_lookup.items():

    if key in found_site_keys:
        continue

    for record in records:
        missing_sites.append(
            record
        )

long_table = pd.DataFrame(
    long_rows
)

samples.to_csv(
    OUT_SAMPLES,
    sep="\t",
    index=False,
    na_rep="NA",
)

sites.to_csv(
    OUT_SITES,
    sep="\t",
    index=False,
    na_rep="NA",
)

long_table.to_csv(
    OUT_LONG,
    sep="\t",
    index=False,
    na_rep="NA",
)

pd.DataFrame(
    missing_sites
).to_csv(
    OUT_MISSING,
    sep="\t",
    index=False,
    na_rep="NA",
)

print(
    "CORE EXPLANATORY OPTIX GENOTYPE EXTRACTION"
)

print(
    "\nSAMPLES BY GROUP"
)

print(
    samples.groupby(
        "group"
    )
    .size()
    .to_string()
)

print(
    "\nPHENOTYPE COMBINATIONS"
)

print(
    samples.loc[
        ~samples[
            "sample"
        ].isin(
            focal_samples
        )
    ]
    .groupby(
        [
            "group",
            "B",
            "R",
            "D",
        ]
    )
    .size()
    .to_string()
)

print(
    "\nSELECTED SITES"
)

print(
    sites.groupby(
        "module"
    )
    .size()
    .to_string()
)

print(
    f"\nRequested samples missing from VCF: "
    f"{len(missing_sample_names)}"
)

for sample in missing_sample_names:
    print(
        f"  {sample}"
    )

print(
    f"Selected module-site records missing from VCF: "
    f"{len(missing_sites)}"
)

print(
    f"\nSample manifest:\n{OUT_SAMPLES}"
)

print(
    f"Site manifest:\n{OUT_SITES}"
)

print(
    f"Genotype table:\n{OUT_LONG}"
)

print(
    f"Missing-site table:\n{OUT_MISSING}"
)
