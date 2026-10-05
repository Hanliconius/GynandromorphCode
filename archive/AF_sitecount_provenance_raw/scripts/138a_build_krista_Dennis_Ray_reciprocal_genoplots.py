#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

VCF_DIR = (
    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
)

METADATA = (
    PROJECT
    / "data"
    / "derived"
    / "region_panel_metadata.tsv"
)

RECIPROCAL_CALLS = (
    PROJECT
    / "analysis"
    / "original_502"
    / "reciprocal_L_R_mixture"
    / "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

OUTPUT_DIR = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
    / "krista_DR_regions_reciprocal_genoplots"
)

SUMMARY_FILE = (
    OUTPUT_DIR
    / "krista_DR_regions_reciprocal_genoplots.summary.tsv"
)

FOCAL_SAMPLES = {
    "TN2L_S1",
    "TN3L_S3",
    "TN2R_S2",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
}

PANELS = {
    "Dennis1": {
        "trait": "D",
        "vcf": "dennis1-1215210-1216124.vcf",
        "expected_snps": 204,
    },
    "Dennis2": {
        "trait": "D",
        "vcf": "dennis2-1218557-1221494.vcf",
        "expected_snps": 638,
    },
    "Ray1": {
        "trait": "R",
        "vcf": "ray1-1156746-1160220.vcf",
        "expected_snps": 850,
    },
    "Ray2": {
        "trait": "R",
        "vcf": "ray2-1174547-1187833.vcf",
        "expected_snps": 2459,
    },
}


def clean(value):

    if value is None or pd.isna(value):
        return ""

    return re.sub(
        r"\s+",
        " ",
        str(value),
    ).strip()


def canonical_gt(value):

    gt = clean(value).replace(
        "|",
        "/",
    )

    if gt == "1/0":
        gt = "0/1"

    if gt in {
        "0/0",
        "0/1",
        "1/1",
    }:
        return gt

    return "./."


def number(value):

    try:
        return float(value)

    except (
        TypeError,
        ValueError,
    ):
        return np.nan


def biological_gt(raw_gt):
    """
    Hmel2 is the non-pattern A/A state in the Dennis and
    Ray regions. ALT is therefore the pattern-associated
    allele for these displays.
    """

    return {
        "0/0": "A/A",
        "0/1": "P/A",
        "1/1": "P/P",
    }.get(
        canonical_gt(raw_gt),
        "N",
    )


required_files = [
    METADATA,
    RECIPROCAL_CALLS,
]

required_files.extend(
    VCF_DIR / panel["vcf"]
    for panel in PANELS.values()
)

for path in required_files:

    if not path.is_file():
        raise RuntimeError(
            f"Missing input: {path}"
        )


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

for column in [
    "Accession",
    "Taxon Name",
    "Country",
    "Data Type",
]:

    metadata[column] = metadata[column].map(
        clean
    )


for trait in [
    "D",
    "R",
]:

    metadata[
        f"{trait}_plot"
    ] = pd.to_numeric(
        metadata[trait],
        errors="coerce",
    )


# Plot-only correction requested by the user.
# The source metadata table is not modified.
metadata.loc[
    metadata["Accession"].eq(
        "ERR3652402"
    ),
    [
        "D_plot",
        "R_plot",
    ],
] = [
    1.0,
    1.0,
]


search_text = (
    metadata
    .apply(
        lambda column: column.map(
            clean
        )
    )
    .agg(
        " ".join,
        axis=1,
    )
)

metadata["piano_key"] = (
    search_text
    .str.contains(
        r"piano\s*keys?",
        case=False,
        regex=True,
    )
)

metadata = (
    metadata
    .loc[
        metadata["Accession"].ne("")
    ]
    .drop_duplicates(
        "Accession"
    )
)

metadata_lookup = (
    metadata
    .set_index(
        "Accession"
    )
    .to_dict(
        "index"
    )
)

piano_key_samples = (
    metadata
    .loc[
        metadata["piano_key"],
        "Accession",
    ]
    .tolist()
)

if len(
    piano_key_samples
) != 5:

    raise RuntimeError(
        "Expected 5 piano-key samples; found "
        f"{len(piano_key_samples)}: "
        + ",".join(
            piano_key_samples
        )
    )


# ------------------------------------------------------------
# Read the posterior >= 0.95 reciprocal TN_Left/TN_Right calls.
# These were fitted without phenotype or expected focal state.
# ------------------------------------------------------------

mixture = pd.read_csv(
    RECIPROCAL_CALLS,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

mixture["pos"] = pd.to_numeric(
    mixture["pos"],
    errors="raise",
).astype(
    int
)

pair_column = next(
    (
        column
        for column in [
            "mixture_pair_canonical",
            "mixture_pair_state",
            "pair_state",
        ]
        if column in mixture.columns
    ),
    None,
)

if pair_column is None:

    candidate_columns = [
        column
        for column in mixture.columns
        if (
            "pair" in column.lower()
            and "original" not in column.lower()
        )
    ]

    if len(
        candidate_columns
    ) != 1:

        raise RuntimeError(
            "Could not identify exactly one reciprocal "
            "pair-state column."
        )

    pair_column = candidate_columns[0]


pair_states = (
    mixture[pair_column]
    .str.split(
        "|",
        regex=False,
        n=1,
        expand=True,
    )
)

mixture["TN_Left"] = (
    pair_states[0]
    .map(
        canonical_gt
    )
)

mixture["TN_Right"] = (
    pair_states[1]
    .map(
        canonical_gt
    )
)

mixture_lookup = {
    (
        row.chrom,
        int(row.pos),
        row.ref,
        row.alt,
    ): {
        "TN_Left": row.TN_Left,
        "TN_Right": row.TN_Right,
    }
    for row in mixture[
        [
            "chrom",
            "pos",
            "ref",
            "alt",
            "TN_Left",
            "TN_Right",
        ]
    ].itertuples(
        index=False
    )
}


def read_vcf(path):

    records = []
    samples = None

    with path.open(
        "rt",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        for line in handle:

            if line.startswith(
                "##"
            ):
                continue

            if line.startswith(
                "#CHROM"
            ):

                samples = (
                    line
                    .rstrip("\n")
                    .split("\t")[9:]
                )

                continue

            if line.startswith(
                "#"
            ):
                continue

            if samples is None:
                raise RuntimeError(
                    "VCF sample header was not found: "
                    f"{path}"
                )

            fields = (
                line
                .rstrip("\n")
                .split("\t")
            )

            chrom = fields[0]
            pos = int(
                fields[1]
            )
            ref = fields[3]
            alt = fields[4]

            if (
                len(ref) != 1
                or len(alt) != 1
                or "," in alt
            ):
                continue

            format_fields = (
                fields[8]
                .split(":")
            )

            calls = {}

            for sample, value in zip(
                samples,
                fields[9:],
            ):

                call = dict(
                    zip(
                        format_fields,
                        value.split(":"),
                    )
                )

                gt = canonical_gt(
                    call.get(
                        "GT"
                    )
                )

                dp = number(
                    call.get(
                        "DP"
                    )
                )

                gq = number(
                    call.get(
                        "GQ"
                    )
                )

                qc_pass = (
                    gt != "./."
                    and np.isfinite(dp)
                    and dp >= 2
                    and np.isfinite(gq)
                    and gq >= 20
                )

                calls[sample] = (
                    gt,
                    call.get(
                        "AD",
                        ".",
                    ),
                    dp,
                    gq,
                    qc_pass,
                )

            records.append(
                (
                    chrom,
                    pos,
                    ref,
                    alt,
                    calls,
                )
            )

    return records


def sample_label(
    sample,
    trait=None,
):

    row = metadata_lookup[
        sample
    ]

    label_values = [
        sample,
        row["Taxon Name"],
        row["Country"],
        row["Data Type"],
    ]

    label = " | ".join(
        value
        for value in label_values
        if value
    )

    if trait is not None:

        phenotype = int(
            row[
                f"{trait}_plot"
            ]
        )

        label = (
            f"{label} | "
            f"{trait}={phenotype}"
        )

    return label


summary_rows = []


for panel_name, specification in PANELS.items():

    trait = specification[
        "trait"
    ]

    expected_snps = specification[
        "expected_snps"
    ]

    vcf_path = (
        VCF_DIR
        / specification["vcf"]
    )

    records = read_vcf(
        vcf_path
    )

    if len(
        records
    ) != expected_snps:

        raise RuntimeError(
            f"{panel_name}: expected "
            f"{expected_snps} SNPs; found "
            f"{len(records)}"
        )

    vcf_samples = list(
        records[0][4]
    )

    missing_metadata = sorted(
        set(vcf_samples)
        - set(metadata_lookup)
    )

    if missing_metadata:

        raise RuntimeError(
            f"{panel_name}: metadata missing for "
            + ",".join(
                missing_metadata[:10]
            )
        )

    comparative_samples = [
        sample
        for sample in vcf_samples
        if (
            sample not in FOCAL_SAMPLES
            and sample not in piano_key_samples
            and metadata_lookup[
                sample
            ][
                f"{trait}_plot"
            ] in {
                0.0,
                1.0,
            }
        )
    ]

    if len(
        comparative_samples
    ) != 481:

        raise RuntimeError(
            f"{panel_name}: expected original "
            "comparative cohort n=481; found "
            f"{len(comparative_samples)}"
        )

    comparative_samples.sort(
        key=lambda sample: (
            int(
                metadata_lookup[
                    sample
                ][
                    f"{trait}_plot"
                ]
            ),
            metadata_lookup[
                sample
            ][
                "Taxon Name"
            ].lower(),
            metadata_lookup[
                sample
            ][
                "Country"
            ].lower(),
            sample,
        )
    )

    order_keys = (
        piano_key_samples
        + [
            "TN_Left|mix",
            "TN_Left|raw",
            "TN_Right|mix",
            "TN_Right|raw",
        ]
        + comparative_samples
    )

    display_order = {
        key: index + 1
        for index, key in enumerate(
            order_keys
        )
    }

    output_rows = []


    def add_row(
        record,
        sample,
        label,
        group,
        source,
        gt,
        ad,
        dp,
        gq,
        qc_pass,
        changed,
        order_key,
    ):

        chrom, pos, ref, alt, _ = record

        if qc_pass:

            display_gt = biological_gt(
                gt
            )

        elif source.startswith(
            "reciprocal"
        ):

            display_gt = "Uncertain"

        else:

            display_gt = "N"

        output_rows.append(
            {
                "panel": panel_name,
                "trait": trait,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "site_label": f"{pos:,}",
                "sample": sample,
                "sample_label": label,
                "plot_group": group,
                "call_source": source,
                "raw_GT": gt,
                "AD": ad,
                "DP": dp,
                "GQ": gq,
                "display_GT": display_gt,
                "changed_from_original": changed,
                "display_order": display_order[
                    order_key
                ],
            }
        )


    for record in records:

        chrom, pos, ref, alt, calls = record

        site_key = (
            chrom,
            pos,
            ref,
            alt,
        )

        for sample in piano_key_samples:

            gt, ad, dp, gq, qc_pass = (
                calls[sample]
            )

            add_row(
                record=record,
                sample=sample,
                label=sample_label(
                    sample
                ),
                group="Butterfly World piano keys",
                source="original VCF",
                gt=gt,
                ad=ad,
                dp=dp,
                gq=gq,
                qc_pass=qc_pass,
                changed=False,
                order_key=sample,
            )

        reciprocal = mixture_lookup.get(
            site_key,
            {},
        )

        for sample in [
            "TN_Left",
            "TN_Right",
        ]:

            gt, ad, dp, gq, qc_pass = (
                calls[sample]
            )

            add_row(
                record=record,
                sample=sample,
                label=(
                    f"{sample} | original VCF call"
                ),
                group="Focal pooled samples",
                source="original VCF",
                gt=gt,
                ad=ad,
                dp=dp,
                gq=gq,
                qc_pass=qc_pass,
                changed=False,
                order_key=f"{sample}|raw",
            )

            new_gt = canonical_gt(
                reciprocal.get(
                    sample,
                    "./.",
                )
            )

            confident = (
                new_gt != "./."
            )

            add_row(
                record=record,
                sample=sample,
                label=(
                    f"{sample} | reciprocal-mixture call"
                ),
                group="Focal pooled samples",
                source=(
                    "reciprocal mixture "
                    "(posterior >= 0.95)"
                ),
                gt=new_gt,
                ad=ad,
                dp=dp,
                gq=gq,
                qc_pass=confident,
                changed=(
                    confident
                    and new_gt != gt
                ),
                order_key=f"{sample}|mix",
            )

        for sample in comparative_samples:

            gt, ad, dp, gq, qc_pass = (
                calls[sample]
            )

            phenotype = int(
                metadata_lookup[
                    sample
                ][
                    f"{trait}_plot"
                ]
            )

            if phenotype == 0:

                group = (
                    f"{panel_name} absent "
                    f"({trait}=0)"
                )

            else:

                group = (
                    f"{panel_name} present "
                    f"({trait}=1)"
                )

            add_row(
                record=record,
                sample=sample,
                label=sample_label(
                    sample,
                    trait,
                ),
                group=group,
                source="original VCF",
                gt=gt,
                ad=ad,
                dp=dp,
                gq=gq,
                qc_pass=qc_pass,
                changed=False,
                order_key=sample,
            )


    output = pd.DataFrame(
        output_rows
    ).sort_values(
        [
            "display_order",
            "pos",
        ],
        kind="stable",
    )

    expected_rows = (
        expected_snps
        * (
            5
            + 4
            + 481
        )
    )

    if len(
        output
    ) != expected_rows:

        raise RuntimeError(
            f"{panel_name}: final table has "
            f"{len(output)} rows; expected "
            f"{expected_rows}"
        )

    output_file = (
        OUTPUT_DIR
        / (
            f"optix_{panel_name}_original_region_"
            "reciprocal_genoplot.long.tsv.gz"
        )
    )

    output.to_csv(
        output_file,
        sep="\t",
        index=False,
        compression={
            "method": "gzip",
            "mtime": 0,
        },
        na_rep="NA",
    )

    summary_rows.append(
        {
            "panel": panel_name,
            "trait": trait,
            "biallelic_SNPs": expected_snps,
            "comparative_samples": 481,
            "piano_key_samples": 5,
            "TN_Left_posterior95": int(
                sum(
                    (
                        output["sample_label"]
                        == (
                            "TN_Left | "
                            "reciprocal-mixture call"
                        )
                    )
                    & (
                        output["display_GT"]
                        != "Uncertain"
                    )
                )
            ),
            "TN_Right_posterior95": int(
                sum(
                    (
                        output["sample_label"]
                        == (
                            "TN_Right | "
                            "reciprocal-mixture call"
                        )
                    )
                    & (
                        output["display_GT"]
                        != "Uncertain"
                    )
                )
            ),
            "reciprocal_cells_changed": int(
                output[
                    "changed_from_original"
                ].sum()
            ),
            "output": str(
                output_file
            ),
        }
    )


summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    SUMMARY_FILE,
    sep="\t",
    index=False,
)

print(
    "KRISTA DENNIS/RAY ORIGINAL-REGION "
    "GENOPLOT TABLES"
)

print()

print(
    summary
    .drop(
        columns="output"
    )
    .to_string(
        index=False
    )
)

print()
print(
    "Plot-only correction: ERR3652402 "
    "D=1, R=1; source metadata unchanged"
)

print()
print(
    f"Summary: {SUMMARY_FILE}"
)

for output_file in summary[
    "output"
]:

    print(
        output_file
    )
