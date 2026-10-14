#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

CHECKPOINTS = (
    PROJECT
    / "analysis"
    / "drop_problematic_reanalysis"
    / "checkpoints"
)

PERFECT_FILE = (
    CHECKPOINTS
    / "perfect_D_R_comparative_SNPs.DP2_GQ20.tsv"
)

MANIFEST_FILE = (
    CHECKPOINTS
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)

METADATA_FILE = (
    PROJECT
    / "data"
    / "derived"
    / "region_panel_metadata.tsv"
)

VCF_FILE = (
    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
    / "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

OUT_LONG = (
    CHECKPOINTS
    / "perfect_D_R_all_comparative_genotypes.long.tsv"
)

OUT_SAMPLES = (
    CHECKPOINTS
    / "perfect_D_R_all_comparative_samples.tsv"
)

OUT_SUMMARY = (
    CHECKPOINTS
    / "perfect_D_R_all_comparative_genotypes.summary.tsv"
)


MODULES = {
    "Dennis": {
        "manifest_module": "Dennis2",
        "trait": "D",
    },
    "Ray": {
        "manifest_module": "Ray2",
        "trait": "R",
    },
}


def clean(value):

    if value is None:
        return ""

    return (
        str(value)
        .replace("\n", " ")
        .strip()
    )


def as_bool(value):

    return clean(value).lower() in {
        "true",
        "t",
        "1",
        "yes",
        "y",
    }


def normalized(value):

    return (
        clean(value)
        .lower()
        .replace(" ", "")
        .replace("_", "")
        .replace("-", "")
    )


def find_column(columns, target):

    target = normalized(target)

    for column in columns:
        if normalized(column) == target:
            return column

    raise RuntimeError(
        f"Could not find column {target}; "
        f"available columns: {list(columns)}"
    )


def canonical_gt(value):

    value = clean(value).replace("|", "/")

    if value in {
        "",
        ".",
        "./.",
    }:
        return None

    parts = value.split("/")

    if (
        len(parts) != 2
        or "." in parts
    ):
        return None

    try:
        alleles = sorted(
            int(part)
            for part in parts
        )
    except ValueError:
        return None

    return "/".join(
        str(allele)
        for allele in alleles
    )


def numeric(value):

    value = clean(value)

    if value in {
        "",
        ".",
    }:
        return None

    try:
        return float(value)
    except ValueError:
        return None


def biological_gt(gt):

    mapping = {
        "0/0": "A/A",
        "0/1": "P/A",
        "1/1": "P/P",
    }

    return mapping.get(
        gt,
        "No call",
    )


for required in (
    PERFECT_FILE,
    MANIFEST_FILE,
    METADATA_FILE,
    VCF_FILE,
):
    if not required.exists():
        raise RuntimeError(
            f"Missing required input: {required}"
        )


perfect = pd.read_csv(
    PERFECT_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

perfect = perfect.loc[
    perfect["module"].isin(
        MODULES
    )
].copy()

perfect["pos"] = pd.to_numeric(
    perfect["pos"],
    errors="raise",
).astype(int)

perfect["site_key"] = (
    perfect["chrom"]
    + ":"
    + perfect["pos"].astype(str)
    + ":"
    + perfect["ref"]
    + ":"
    + perfect["alt"]
)

perfect = perfect.drop_duplicates(
    [
        "module",
        "site_key",
    ]
)


manifest = pd.read_csv(
    MANIFEST_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

allowed = {}

for module, specification in MODULES.items():

    subset = manifest.loc[
        (
            manifest["cohort"]
            == "WGS_plus_CS"
        )
        & (
            manifest["module"]
            == specification["manifest_module"]
        )
        & manifest["keep"].map(as_bool)
    ].copy()

    allowed[module] = list(
        dict.fromkeys(
            subset["sample"].map(clean)
        )
    )


metadata = pd.read_csv(
    METADATA_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

accession_column = find_column(
    metadata.columns,
    "Accession",
)

taxon_column = find_column(
    metadata.columns,
    "Taxon Name",
)

country_column = find_column(
    metadata.columns,
    "Country",
)

data_type_column = find_column(
    metadata.columns,
    "Data Type",
)

trait_columns = {
    trait: find_column(
        metadata.columns,
        trait,
    )
    for trait in (
        "D",
        "R",
    )
}

metadata_records = {}

for _, row in metadata.iterrows():

    accession = clean(
        row[accession_column]
    )

    if not accession:
        continue

    record = {
        "taxon": " ".join(
            clean(
                row[taxon_column]
            ).split()
        ),
        "country": " ".join(
            clean(
                row[country_column]
            ).split()
        ),
        "data_type": clean(
            row[data_type_column]
        ),
    }

    for trait, column in trait_columns.items():

        value = numeric(
            row[column]
        )

        record[trait] = (
            int(value)
            if value in {
                0.0,
                1.0,
            }
            else None
        )

    metadata_records[accession] = record


sample_rows = []

for module, samples in allowed.items():

    trait = MODULES[module]["trait"]

    for sample in samples:

        record = metadata_records.get(
            sample
        )

        if record is None:
            raise RuntimeError(
                f"{module}: no metadata for {sample}"
            )

        phenotype = record[trait]

        if phenotype not in {
            0,
            1,
        }:
            raise RuntimeError(
                f"{module}: missing {trait} phenotype "
                f"for retained sample {sample}"
            )

        sample_rows.append(
            {
                "module": module,
                "trait": trait,
                "sample": sample,
                "phenotype": phenotype,
                "taxon": record["taxon"],
                "country": record["country"],
                "data_type": record["data_type"],
            }
        )

sample_table = pd.DataFrame(
    sample_rows
)

sample_table = sample_table.sort_values(
    [
        "module",
        "phenotype",
        "taxon",
        "data_type",
        "sample",
    ]
)

sample_table["plot_order"] = (
    sample_table
    .groupby("module")
    .cumcount()
    + 1
)

sample_table["sample_label"] = (
    sample_table["sample"]
    + " | "
    + sample_table["taxon"]
    + " | "
    + sample_table["data_type"]
)

sample_table.to_csv(
    OUT_SAMPLES,
    sep="\t",
    index=False,
)


target_to_modules = {}

for _, row in perfect.iterrows():

    key = (
        row["chrom"],
        int(row["pos"]),
        row["ref"],
        row["alt"],
    )

    target_to_modules.setdefault(
        key,
        []
    ).append(
        row["module"]
    )

target_keys = set(
    target_to_modules
)

union_samples = sorted(
    set().union(
        *(
            set(samples)
            for samples in allowed.values()
        )
    )
)

sample_indices = None
recovered = {}
long_rows = []


with VCF_FILE.open(
    "rt",
    encoding="utf-8",
    errors="replace",
) as handle:

    for line in handle:

        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):

            header = line.rstrip("\n").split("\t")
            vcf_samples = header[9:]

            sample_to_index = {
                sample: index
                for index, sample in enumerate(
                    vcf_samples
                )
            }

            missing_samples = [
                sample
                for sample in union_samples
                if sample not in sample_to_index
            ]

            if missing_samples:
                raise RuntimeError(
                    "Retained samples missing from VCF: "
                    + ",".join(
                        missing_samples
                    )
                )

            sample_indices = {
                sample: sample_to_index[sample]
                for sample in union_samples
            }

            continue

        if sample_indices is None:
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

        key = (
            chrom,
            pos,
            ref,
            alt,
        )

        if key not in target_keys:
            continue

        recovered[key] = True

        format_names = fields[8].split(":")

        for module in target_to_modules[key]:

            perfect_row = perfect.loc[
                (
                    perfect["module"]
                    == module
                )
                & (
                    perfect["chrom"]
                    == chrom
                )
                & (
                    perfect["pos"]
                    == pos
                )
                & (
                    perfect["ref"]
                    == ref
                )
                & (
                    perfect["alt"]
                    == alt
                )
            ].iloc[0]

            trait = MODULES[module]["trait"]

            module_samples = sample_table.loc[
                sample_table["module"]
                == module
            ]

            for _, sample_record in module_samples.iterrows():

                sample = sample_record["sample"]
                index = sample_indices[sample]

                values = fields[
                    9 + index
                ].split(":")

                call = {
                    name: (
                        values[field_index]
                        if field_index < len(values)
                        else ""
                    )
                    for field_index, name in enumerate(
                        format_names
                    )
                }

                gt = canonical_gt(
                    call.get(
                        "GT",
                        "",
                    )
                )

                dp = numeric(
                    call.get(
                        "DP",
                        "",
                    )
                )

                gq = numeric(
                    call.get(
                        "GQ",
                        "",
                    )
                )

                ad = clean(
                    call.get(
                        "AD",
                        "",
                    )
                )

                qc_pass = (
                    gt is not None
                    and dp is not None
                    and dp >= 2
                    and gq is not None
                    and gq >= 20
                )

                biological = (
                    biological_gt(gt)
                    if qc_pass
                    else "No call"
                )

                predicted_presence = (
                    int(
                        biological in {
                            "P/A",
                            "P/P",
                        }
                    )
                    if qc_pass
                    else None
                )

                phenotype = int(
                    sample_record["phenotype"]
                )

                discordant = (
                    qc_pass
                    and predicted_presence != phenotype
                )

                long_rows.append(
                    {
                        "module": module,
                        "trait": trait,
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "site_key": perfect_row["site_key"],
                        "perfect_both_schemes": perfect_row.get(
                            "perfect_both_schemes",
                            "",
                        ),
                        "sample": sample,
                        "sample_label": sample_record["sample_label"],
                        "plot_order": sample_record["plot_order"],
                        "phenotype": phenotype,
                        "taxon": sample_record["taxon"],
                        "country": sample_record["country"],
                        "data_type": sample_record["data_type"],
                        "GT": gt if gt is not None else "",
                        "biological_GT": biological,
                        "AD": ad,
                        "DP": dp if dp is not None else "",
                        "GQ": gq if gq is not None else "",
                        "DP2_GQ20_pass": qc_pass,
                        "predicted_presence": (
                            predicted_presence
                            if predicted_presence is not None
                            else ""
                        ),
                        "callable_discordance": discordant,
                    }
                )

        if len(recovered) == len(target_keys):
            break


missing_sites = target_keys.difference(
    recovered
)

if missing_sites:
    raise RuntimeError(
        "Perfect sites missing from broad VCF: "
        + ",".join(
            ":".join(
                map(str, key)
            )
            for key in sorted(
                missing_sites
            )
        )
    )


long_table = pd.DataFrame(
    long_rows
)

long_table = long_table.sort_values(
    [
        "module",
        "pos",
        "plot_order",
    ]
)

long_table.to_csv(
    OUT_LONG,
    sep="\t",
    index=False,
)


summary_rows = []

for module in MODULES:

    subset = long_table.loc[
        long_table["module"]
        == module
    ]

    callable_subset = subset.loc[
        subset["DP2_GQ20_pass"]
    ]

    discordances = int(
        callable_subset[
            "callable_discordance"
        ].sum()
    )

    summary_rows.append(
        {
            "module": module,
            "samples": subset["sample"].nunique(),
            "perfect_SNPs": subset["site_key"].nunique(),
            "genotype_cells": len(subset),
            "callable_cells": len(callable_subset),
            "call_rate": (
                len(callable_subset)
                / len(subset)
            ),
            "callable_discordances": discordances,
            "validation": (
                "perfect_reproduced"
                if discordances == 0
                else "ERROR_discordance_recovered"
            ),
        }
    )

summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False,
)


print("PERFECT D/R ALL-COMPARATIVE GENOTYPE EXTRACTION")
print()
print(summary.to_string(index=False))
print()
print(f"Long genotype table: {OUT_LONG}")
print(f"Sample manifest: {OUT_SAMPLES}")
print(f"Summary: {OUT_SUMMARY}")
