#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

VCF = (
    PROJECT
    / "external/krista_genotype_plots_original/Genotype_Plots"
    / "optix_band-1112651-1114314.vcf"
)

METADATA = (
    PROJECT
    / "data/derived/region_panel_metadata.tsv"
)

MIXTURE = (
    PROJECT
    / "analysis/original_502/reciprocal_L_R_mixture"
    / "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

OUTPUT = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "optix_band_original_region_reciprocal_genoplot.long.tsv.gz"
)

SUMMARY = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "optix_band_original_region_reciprocal_genoplot.summary.tsv"
)

FOCAL_SAMPLES = {
    "TN2L_S1",
    "TN3L_S3",
    "TN2R_S2",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
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


def numeric(value):
    try:
        return float(value)
    except (
        TypeError,
        ValueError,
    ):
        return np.nan


def band_biological_gt(raw_gt):
    """
    The Band reference state is P/P:
      REF/REF -> P/P
      REF/ALT -> P/A
      ALT/ALT -> A/A
    """

    return {
        "0/0": "P/P",
        "0/1": "P/A",
        "1/1": "A/A",
    }.get(
        canonical_gt(raw_gt),
        "N",
    )


for required_file in [
    VCF,
    METADATA,
    MIXTURE,
]:
    if not required_file.is_file():
        raise RuntimeError(
            f"Missing input: {required_file}"
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

metadata["B_plot"] = pd.to_numeric(
    metadata["B"],
    errors="coerce",
)

# Plot-only correction supplied by the user.
# The source metadata file is not modified.
metadata.loc[
    metadata["Accession"].eq(
        "ERR3652402"
    ),
    "B_plot",
] = 0.0

searchable = metadata.apply(
    lambda column: column.map(clean)
).agg(
    " ".join,
    axis=1,
)

metadata["piano_key"] = searchable.str.contains(
    r"piano\s*keys?",
    case=False,
    regex=True,
)

metadata = metadata.loc[
    metadata["Accession"].ne("")
].drop_duplicates(
    "Accession"
)

metadata_lookup = metadata.set_index(
    "Accession"
).to_dict(
    "index"
)

piano_samples = metadata.loc[
    metadata["piano_key"],
    "Accession",
].tolist()

if len(piano_samples) != 5:
    raise RuntimeError(
        "Expected five piano-key samples; found "
        f"{len(piano_samples)}: "
        + ",".join(piano_samples)
    )


# Load the posterior >= 0.95 reciprocal calls.
mixture = pd.read_csv(
    MIXTURE,
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

    candidates = [
        column
        for column in mixture.columns
        if (
            "pair"
            in column.lower()
            and "original"
            not in column.lower()
        )
    ]

    if len(candidates) != 1:
        raise RuntimeError(
            "Could not identify one reciprocal pair-state column."
        )

    pair_column = candidates[0]

pair_parts = mixture[
    pair_column
].str.split(
    "|",
    regex=False,
    n=1,
    expand=True,
)

mixture["TN_Left"] = pair_parts[0].map(
    canonical_gt
)

mixture["TN_Right"] = pair_parts[1].map(
    canonical_gt
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


# Read all biallelic SNPs from Krista's original Band VCF.
records = []

with VCF.open(
    "rt",
    encoding="utf-8",
    errors="replace",
) as handle:

    sample_names = None

    for line in handle:

        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            sample_names = line.rstrip(
                "\n"
            ).split(
                "\t"
            )[9:]
            continue

        if line.startswith("#"):
            continue

        if sample_names is None:
            raise RuntimeError(
                "VCF sample header was not found."
            )

        fields = line.rstrip(
            "\n"
        ).split(
            "\t"
        )

        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]

        if (
            len(ref) != 1
            or len(alt) != 1
            or "," in alt
        ):
            continue

        format_fields = fields[8].split(
            ":"
        )

        calls = {}

        for sample, sample_value in zip(
            sample_names,
            fields[9:],
        ):

            call = dict(
                zip(
                    format_fields,
                    sample_value.split(":"),
                )
            )

            raw_gt = canonical_gt(
                call.get("GT")
            )

            dp = numeric(
                call.get("DP")
            )

            gq = numeric(
                call.get("GQ")
            )

            qc_pass = (
                raw_gt != "./."
                and np.isfinite(dp)
                and dp >= 2
                and np.isfinite(gq)
                and gq >= 20
            )

            calls[sample] = (
                raw_gt,
                call.get("AD", "."),
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


if len(records) != 153:
    raise RuntimeError(
        "Expected 153 biallelic Band SNPs; found "
        f"{len(records)}."
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
        "Metadata missing for: "
        + ",".join(
            missing_metadata[:10]
        )
    )


# Reproduce the original 481-sample Band scoring cohort.
comparative_samples = [
    sample
    for sample in vcf_samples
    if (
        sample not in FOCAL_SAMPLES
        and sample not in piano_samples
        and metadata_lookup[sample]["B_plot"]
            in {0.0, 1.0}
    )
]

if len(comparative_samples) != 481:
    raise RuntimeError(
        "Expected the original Band cohort n=481; found "
        f"{len(comparative_samples)}."
    )

comparative_samples.sort(
    key=lambda sample: (
        int(
            metadata_lookup[sample]["B_plot"]
        ),
        metadata_lookup[sample][
            "Taxon Name"
        ].lower(),
        metadata_lookup[sample][
            "Country"
        ].lower(),
        sample,
    )
)


def sample_label(
    sample,
    include_phenotype=False,
):
    row = metadata_lookup[sample]

    values = [
        sample,
        row["Taxon Name"],
        row["Country"],
        row["Data Type"],
    ]

    label = " | ".join(
        value
        for value in values
        if value
    )

    if include_phenotype:
        label += (
            f" | B={int(row['B_plot'])}"
        )

    return label


order_keys = (
    piano_samples
    + [
        "TN_Left|mixture",
        "TN_Left|original",
        "TN_Right|mixture",
        "TN_Right|original",
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
    call_source,
    raw_gt,
    ad,
    dp,
    gq,
    qc_pass,
    changed,
    order_key,
):
    chrom, pos, ref, alt, _ = record

    if qc_pass:
        displayed_gt = band_biological_gt(
            raw_gt
        )
    elif call_source.startswith(
        "reciprocal"
    ):
        displayed_gt = "Uncertain"
    else:
        displayed_gt = "N"

    output_rows.append(
        {
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
            "site_label": f"{pos:,}",
            "sample": sample,
            "sample_label": label,
            "plot_group": group,
            "call_source": call_source,
            "raw_GT": raw_gt,
            "AD": ad,
            "DP": dp,
            "GQ": gq,
            "display_GT": displayed_gt,
            "changed_from_original":
                changed,
            "display_order":
                display_order[order_key],
        }
    )


for record in records:

    chrom, pos, ref, alt, calls = record

    key = (
        chrom,
        pos,
        ref,
        alt,
    )

    for sample in piano_samples:

        raw_gt, ad, dp, gq, qc_pass = calls[
            sample
        ]

        add_row(
            record=record,
            sample=sample,
            label=sample_label(sample),
            group="Butterfly World piano keys",
            call_source="original VCF",
            raw_gt=raw_gt,
            ad=ad,
            dp=dp,
            gq=gq,
            qc_pass=qc_pass,
            changed=False,
            order_key=sample,
        )

    reciprocal_pair = mixture_lookup.get(
        key,
        {},
    )

    for sample in [
        "TN_Left",
        "TN_Right",
    ]:

        raw_gt, ad, dp, gq, qc_pass = calls[
            sample
        ]

        add_row(
            record=record,
            sample=sample,
            label=f"{sample} | original VCF call",
            group="Focal pooled samples",
            call_source="original VCF",
            raw_gt=raw_gt,
            ad=ad,
            dp=dp,
            gq=gq,
            qc_pass=qc_pass,
            changed=False,
            order_key=f"{sample}|original",
        )

        recalled_gt = canonical_gt(
            reciprocal_pair.get(
                sample,
                "./.",
            )
        )

        confident = (
            recalled_gt != "./."
        )

        add_row(
            record=record,
            sample=sample,
            label=(
                f"{sample} | "
                "reciprocal-mixture call"
            ),
            group="Focal pooled samples",
            call_source=(
                "reciprocal mixture "
                "(posterior >= 0.95)"
            ),
            raw_gt=recalled_gt,
            ad=ad,
            dp=dp,
            gq=gq,
            qc_pass=confident,
            changed=(
                confident
                and recalled_gt != raw_gt
            ),
            order_key=f"{sample}|mixture",
        )

    for sample in comparative_samples:

        raw_gt, ad, dp, gq, qc_pass = calls[
            sample
        ]

        phenotype = int(
            metadata_lookup[sample][
                "B_plot"
            ]
        )

        if phenotype == 0:
            group = "Band absent (B=0)"
        else:
            group = "Band present (B=1)"

        add_row(
            record=record,
            sample=sample,
            label=sample_label(
                sample,
                include_phenotype=True,
            ),
            group=group,
            call_source="original VCF",
            raw_gt=raw_gt,
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
    153
    * (
        5
        + 4
        + 481
    )
)

if len(output) != expected_rows:
    raise RuntimeError(
        f"Final table has {len(output)} rows; "
        f"expected {expected_rows}."
    )

output.to_csv(
    OUTPUT,
    sep="\t",
    index=False,
    compression={
        "method": "gzip",
        "mtime": 0,
    },
    na_rep="NA",
)


summary = pd.DataFrame(
    [
        (
            "biallelic_SNPs",
            153,
        ),
        (
            "comparative_samples",
            481,
        ),
        (
            "piano_key_samples",
            5,
        ),
        (
            "TN_Left_posterior95",
            int(
                sum(
                    (
                        output["sample_label"]
                        == "TN_Left | reciprocal-mixture call"
                    )
                    & (
                        output["display_GT"]
                        != "Uncertain"
                    )
                )
            ),
        ),
        (
            "TN_Right_posterior95",
            int(
                sum(
                    (
                        output["sample_label"]
                        == "TN_Right | reciprocal-mixture call"
                    )
                    & (
                        output["display_GT"]
                        != "Uncertain"
                    )
                )
            ),
        ),
        (
            "reciprocal_cells_changed",
            int(
                output[
                    "changed_from_original"
                ].sum()
            ),
        ),
    ],
    columns=[
        "statistic",
        "value",
    ],
)

summary.to_csv(
    SUMMARY,
    sep="\t",
    index=False,
)


print(
    "KRISTA ORIGINAL BAND-REGION GENOPLOT TABLE\n"
)

print(
    summary.to_string(
        index=False
    )
)

print(
    "\nPiano keys: "
    + ",".join(
        piano_samples
    )
)

print(
    "Plot-only correction: "
    "ERR3652402 B=0; "
    "source metadata unchanged"
)

print(
    f"\nLong table: {OUTPUT}"
)

print(
    f"Summary: {SUMMARY}"
)
