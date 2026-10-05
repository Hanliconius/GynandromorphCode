#!/usr/bin/env python3

from pathlib import Path
import hashlib
import re

import numpy as np
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

SOURCE = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "optix_R_all_BA90_comparative_selected_genotypes.long.tsv"
)

MIXTURE = (
    PROJECT
    / "analysis/original_502/reciprocal_L_R_mixture"
    / "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

VCF = (
    PROJECT
    / "external/krista_genotype_plots_original/Genotype_Plots"
    / "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

METADATA = (
    PROJECT
    / "data/derived/region_panel_metadata.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis/original_502/checkpoints"
)

OUTPUT = (
    OUTDIR
    / "optix_R_BA90_reciprocal_genoplot.long.tsv.gz"
)

SUMMARY = (
    OUTDIR
    / "optix_R_BA90_reciprocal_genoplot.summary.tsv"
)

MANIFEST = (
    OUTDIR
    / "optix_R_BA90_reciprocal_genoplot.manifest.tsv"
)

KEY = [
    "chrom",
    "pos",
    "ref",
    "alt",
]


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


def alt_dosage(value):
    return {
        "0/0": 0.0,
        "0/1": 0.5,
        "1/1": 1.0,
    }.get(
        canonical_gt(value),
        np.nan,
    )


def biological_call(
    raw_gt,
    pattern_is_alt,
):
    dosage = alt_dosage(
        raw_gt
    )

    if not np.isfinite(
        dosage
    ):
        return (
            np.nan,
            "N",
        )

    if pattern_is_alt:
        pattern_dosage = dosage
    else:
        pattern_dosage = (
            1.0 - dosage
        )

    biological_gt = {
        0.0: "A/A",
        0.5: "P/A",
        1.0: "P/P",
    }[
        pattern_dosage
    ]

    return (
        pattern_dosage,
        biological_gt,
    )


def numeric(value):
    try:
        return float(
            value
        )
    except (
        TypeError,
        ValueError,
    ):
        return np.nan


def sha256_file(path):
    digest = hashlib.sha256()

    with path.open(
        "rb"
    ) as handle:

        for block in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(
                block
            )

    return digest.hexdigest()


for required_file in [
    SOURCE,
    MIXTURE,
    VCF,
    METADATA,
]:
    if not required_file.is_file():
        raise RuntimeError(
            f"Missing input: {required_file}"
        )


OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


source = pd.read_csv(
    SOURCE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

source["pos"] = pd.to_numeric(
    source["pos"],
    errors="raise",
).astype(
    int
)


sites = source[
    KEY
    + [
        "site_label",
        "minimum_BA_all_analyses",
    ]
].drop_duplicates()

sites = sites.sort_values(
    KEY,
    kind="stable",
).reset_index(
    drop=True
)


if len(
    sites
) != 255:
    raise RuntimeError(
        "Expected 255 frozen BA90 SNPs; "
        f"recovered {len(sites)}."
    )


comparative = source.loc[
    source[
        "sample_type"
    ].str.lower().eq(
        "comparative"
    )
].copy()

comparative_samples = list(
    dict.fromkeys(
        comparative[
            "sample"
        ].tolist()
    )
)


if len(
    comparative_samples
) != 372:
    raise RuntimeError(
        "Expected 372 comparative samples; "
        f"recovered {len(comparative_samples)}."
    )


focal = source.loc[
    source[
        "sample"
    ].isin(
        [
            "TN_Left",
            "TN_Right",
        ]
    )
].copy()


for sample in [
    "TN_Left",
    "TN_Right",
]:
    recovered = int(
        focal[
            "sample"
        ].eq(
            sample
        ).sum()
    )

    if recovered != 255:
        raise RuntimeError(
            f"{sample}: expected 255 original rows; "
            f"recovered {recovered}."
        )


# Infer whether REF or ALT is the pattern-associated allele
# independently at every selected SNP.
pattern_is_alt = {}

for key, rows in source.groupby(
    KEY,
    sort=False,
):

    evidence = []

    for row in rows.itertuples(
        index=False
    ):

        raw_dosage = alt_dosage(
            row.raw_GT
        )

        pattern_dosage = numeric(
            row.pattern_dosage
        )

        if (
            not np.isfinite(
                raw_dosage
            )
            or not np.isfinite(
                pattern_dosage
            )
            or np.isclose(
                raw_dosage,
                0.5,
            )
        ):
            continue

        if np.isclose(
            pattern_dosage,
            raw_dosage,
        ):
            evidence.append(
                True
            )

        elif np.isclose(
            pattern_dosage,
            1.0 - raw_dosage,
        ):
            evidence.append(
                False
            )

    if (
        not evidence
        or len(
            set(
                evidence
            )
        ) != 1
    ):
        raise RuntimeError(
            "Could not infer one pattern-allele orientation "
            f"for {key}."
        )

    pattern_is_alt[
        (
            str(
                key[0]
            ),
            int(
                key[1]
            ),
            str(
                key[2]
            ),
            str(
                key[3]
            ),
        )
    ] = evidence[0]


# The posterior95 table contains only confident reciprocal calls.
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

    if len(
        candidates
    ) != 1:
        raise RuntimeError(
            "Could not identify the mixture pair-state column. "
            "Columns were: "
            + ",".join(
                mixture.columns
            )
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

mixture[
    "TN_Left"
] = pair_parts[
    0
].map(
    canonical_gt
)

mixture[
    "TN_Right"
] = pair_parts[
    1
].map(
    canonical_gt
)

mixture = mixture[
    KEY
    + [
        "TN_Left",
        "TN_Right",
    ]
].drop_duplicates(
    KEY
)


# Identify the five previously excluded Butterfly World piano-key
# samples from metadata.
metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

searchable = metadata.apply(
    lambda column: column.map(
        clean
    )
).agg(
    " ".join,
    axis=1,
)

piano_metadata = metadata.loc[
    searchable.str.contains(
        r"piano\s*keys?",
        case=False,
        regex=True,
    )
].copy()

piano_metadata[
    "Accession"
] = piano_metadata[
    "Accession"
].map(
    clean
)

piano_metadata = piano_metadata.loc[
    piano_metadata[
        "Accession"
    ].ne(
        ""
    )
].drop_duplicates(
    "Accession"
)


if len(
    piano_metadata
) != 5:
    raise RuntimeError(
        "Expected five piano-key samples; recovered "
        f"{len(piano_metadata)}: "
        + ",".join(
            piano_metadata[
                "Accession"
            ]
        )
    )


piano_samples = piano_metadata[
    "Accession"
].tolist()

piano_lookup = piano_metadata.set_index(
    "Accession"
).to_dict(
    "index"
)


site_lookup = {
    (
        str(
            row.chrom
        ),
        int(
            row.pos
        ),
        str(
            row.ref
        ),
        str(
            row.alt
        ),
    ): row
    for row in sites.itertuples(
        index=False
    )
}

selected_positions = {
    (
        key[0],
        key[1],
    )
    for key in site_lookup
}


piano_rows = []
found_sites = set()


# Stream the broad VCF once and retain only the 255 selected
# positions for the five piano-key samples.
with VCF.open(
    "rt",
    encoding="utf-8",
    errors="replace",
) as handle:

    sample_indexes = None

    for line in handle:

        if line.startswith(
            "##"
        ):
            continue

        if line.startswith(
            "#CHROM"
        ):
            header = line.rstrip(
                "\n"
            ).split(
                "\t"
            )

            missing_samples = sorted(
                set(
                    piano_samples
                )
                - set(
                    header[9:]
                )
            )

            if missing_samples:
                raise RuntimeError(
                    "Piano-key samples missing from VCF: "
                    + ",".join(
                        missing_samples
                    )
                )

            sample_indexes = {
                sample: header.index(
                    sample
                )
                for sample in piano_samples
            }

            continue

        if line.startswith(
            "#"
        ):
            continue

        fields = line.rstrip(
            "\n"
        ).split(
            "\t"
        )

        chrom = fields[0]
        pos = int(
            fields[1]
        )

        if (
            chrom,
            pos,
        ) not in selected_positions:
            continue

        key = (
            chrom,
            pos,
            fields[3],
            fields[4],
        )

        if key not in site_lookup:
            continue

        found_sites.add(
            key
        )

        format_fields = fields[
            8
        ].split(
            ":"
        )

        site = site_lookup[
            key
        ]

        for sample in piano_samples:

            values = fields[
                sample_indexes[
                    sample
                ]
            ].split(
                ":"
            )

            call = dict(
                zip(
                    format_fields,
                    values,
                )
            )

            raw_gt = canonical_gt(
                call.get(
                    "GT",
                    "./.",
                )
            )

            dp = numeric(
                call.get(
                    "DP"
                )
            )

            gq = numeric(
                call.get(
                    "GQ"
                )
            )

            qc_pass = (
                raw_gt != "./."
                and np.isfinite(
                    dp
                )
                and dp >= 2
                and np.isfinite(
                    gq
                )
                and gq >= 20
            )

            dosage, biological_gt = biological_call(
                raw_gt,
                pattern_is_alt[
                    key
                ],
            )

            sample_metadata = piano_lookup[
                sample
            ]

            taxon = clean(
                sample_metadata.get(
                    "Taxon Name"
                )
            )

            data_type = clean(
                sample_metadata.get(
                    "Data Type"
                )
            )

            piano_rows.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": fields[3],
                    "alt": fields[4],
                    "site_label": site.site_label,
                    "minimum_BA_all_analyses":
                        site.minimum_BA_all_analyses,
                    "sample": sample,
                    "sample_label":
                        f"{sample} | {taxon} | {data_type}",
                    "plot_group":
                        "Butterfly World piano keys",
                    "call_source":
                        "original VCF",
                    "raw_GT": raw_gt,
                    "AD": call.get(
                        "AD",
                        ".",
                    ),
                    "DP": dp,
                    "GQ": gq,
                    "display_GT":
                        biological_gt
                        if qc_pass
                        else "N",
                    "pattern_dosage":
                        dosage,
                    "changed_from_original":
                        False,
                }
            )


missing_sites = (
    set(
        site_lookup
    )
    - found_sites
)

if missing_sites:
    raise RuntimeError(
        "Broad VCF lacked "
        f"{len(missing_sites)} selected SNP records."
    )


output_columns = [
    "chrom",
    "pos",
    "ref",
    "alt",
    "site_label",
    "minimum_BA_all_analyses",
    "sample",
    "sample_label",
    "plot_group",
    "call_source",
    "raw_GT",
    "AD",
    "DP",
    "GQ",
    "display_GT",
    "pattern_dosage",
    "changed_from_original",
]


# Comparative genotypes remain unchanged.
comparative_output = comparative.copy()

comparative_output[
    "plot_group"
] = comparative_output[
    "row_group"
]

comparative_output[
    "call_source"
] = "original VCF"

comparative_output[
    "changed_from_original"
] = False

comparative_output = comparative_output[
    output_columns
]


# Retain the original focal rows.
original_focal = focal.copy()

original_focal[
    "plot_group"
] = "Focal pooled samples"

original_focal[
    "call_source"
] = "original VCF"

original_focal[
    "changed_from_original"
] = False

original_focal[
    "sample_label"
] = original_focal[
    "sample"
].map(
    {
        "TN_Left":
            "TN_Left | original VCF call",
        "TN_Right":
            "TN_Right | original VCF call",
    }
)

original_focal = original_focal[
    output_columns
]


# Add a second focal row for each side using reciprocal calls.
recalled_focal = focal.merge(
    mixture,
    on=KEY,
    how="left",
    validate="many_to_one",
)

recalled_focal[
    "original_raw_GT"
] = recalled_focal[
    "raw_GT"
].map(
    canonical_gt
)

recalled_focal[
    "raw_GT"
] = np.where(
    recalled_focal[
        "sample"
    ].eq(
        "TN_Left"
    ),
    recalled_focal[
        "TN_Left"
    ],
    recalled_focal[
        "TN_Right"
    ],
)

recalled_focal[
    "raw_GT"
] = pd.Series(
    recalled_focal[
        "raw_GT"
    ],
    index=recalled_focal.index,
).map(
    canonical_gt
)

recalled_focal[
    "confident"
] = recalled_focal[
    "raw_GT"
].ne(
    "./."
)

recalled_focal[
    "changed_from_original"
] = (
    recalled_focal[
        "confident"
    ]
    & recalled_focal[
        "raw_GT"
    ].ne(
        recalled_focal[
            "original_raw_GT"
        ]
    )
)


recalled_calls = []

for row in recalled_focal.itertuples(
    index=False
):

    key = (
        row.chrom,
        int(
            row.pos
        ),
        row.ref,
        row.alt,
    )

    recalled_calls.append(
        biological_call(
            row.raw_GT,
            pattern_is_alt[
                key
            ],
        )
    )


recalled_focal[
    "pattern_dosage"
] = [
    result[0]
    for result in recalled_calls
]

recalled_focal[
    "display_GT"
] = [
    result[1]
    if confident
    else "Uncertain"
    for result, confident in zip(
        recalled_calls,
        recalled_focal[
            "confident"
        ],
    )
]

recalled_focal[
    "plot_group"
] = "Focal pooled samples"

recalled_focal[
    "call_source"
] = (
    "reciprocal mixture "
    "(posterior >= 0.95)"
)

recalled_focal[
    "sample_label"
] = recalled_focal[
    "sample"
].map(
    {
        "TN_Left":
            "TN_Left | reciprocal-mixture call",
        "TN_Right":
            "TN_Right | reciprocal-mixture call",
    }
)

recalled_focal = recalled_focal[
    output_columns
]


piano_output = pd.DataFrame(
    piano_rows
)[
    output_columns
]


# Explicit top-to-bottom ordering:
# piano keys, focal calls, then the frozen comparative order.
row_order = (
    piano_output.drop_duplicates(
        "sample"
    )[
        "sample_label"
    ].tolist()
    + [
        "TN_Left | reciprocal-mixture call",
        "TN_Left | original VCF call",
        "TN_Right | reciprocal-mixture call",
        "TN_Right | original VCF call",
    ]
    + comparative_output.drop_duplicates(
        "sample"
    )[
        "sample_label"
    ].tolist()
)


if len(
    row_order
) != len(
    set(
        row_order
    )
):
    raise RuntimeError(
        "Sample labels are not unique."
    )


order_lookup = {
    label: index + 1
    for index, label in enumerate(
        row_order
    )
}


combined = pd.concat(
    [
        piano_output,
        recalled_focal,
        original_focal,
        comparative_output,
    ],
    ignore_index=True,
)

combined[
    "display_order"
] = combined[
    "sample_label"
].map(
    order_lookup
)

combined = combined.sort_values(
    [
        "display_order",
        "pos",
        "ref",
        "alt",
    ],
    kind="stable",
)


expected_rows = (
    255
    * (
        5
        + 4
        + 372
    )
)

if len(
    combined
) != expected_rows:
    raise RuntimeError(
        f"Final table has {len(combined)} rows; "
        f"expected {expected_rows}."
    )


combined.to_csv(
    OUTPUT,
    sep="\t",
    index=False,
    compression={
        "method": "gzip",
        "mtime": 0,
    },
    na_rep="NA",
)


focal_summary = recalled_focal.groupby(
    "sample",
    sort=False,
).agg(
    selected_SNPs=(
        "pos",
        "size",
    ),
    posterior95_calls=(
        "display_GT",
        lambda values:
            int(
                sum(
                    values
                    != "Uncertain"
                )
            ),
    ),
    changed_from_original=(
        "changed_from_original",
        "sum",
    ),
).reset_index()


piano_summary = piano_output.groupby(
    "sample",
    sort=False,
).agg(
    selected_SNPs=(
        "pos",
        "size",
    ),
    QC_calls=(
        "display_GT",
        lambda values:
            int(
                sum(
                    values
                    != "N"
                )
            ),
    ),
).reset_index()


summary = pd.concat(
    [
        focal_summary.assign(
            section="focal"
        ),
        piano_summary.assign(
            section="piano_key"
        ),
    ],
    ignore_index=True,
)

summary.to_csv(
    SUMMARY,
    sep="\t",
    index=False,
    na_rep="NA",
)


manifest = pd.DataFrame(
    [
        (
            "analysis_version",
            "1.0.0",
        ),
        (
            "frozen_selection",
            "255 Ray SNPs with minimum comparative BA >= 0.90",
        ),
        (
            "comparative_samples",
            "372; unchanged",
        ),
        (
            "focal_recall",
            "joint reciprocal-mixture posterior >= 0.95",
        ),
        (
            "focal_below_threshold",
            "displayed as Uncertain",
        ),
        (
            "piano_key_QC",
            "DP >= 2 and GQ >= 20",
        ),
        (
            "focal_or_piano_used_for_selection",
            "False",
        ),
        (
            "source_SHA256",
            sha256_file(
                SOURCE
            ),
        ),
        (
            "mixture_posterior95_SHA256",
            sha256_file(
                MIXTURE
            ),
        ),
        (
            "VCF_SHA256",
            sha256_file(
                VCF
            ),
        ),
        (
            "output_SHA256",
            sha256_file(
                OUTPUT
            ),
        ),
    ],
    columns=[
        "field",
        "value",
    ],
)

manifest.to_csv(
    MANIFEST,
    sep="\t",
    index=False,
)


print(
    "R BA90 RECIPROCAL GENOPLOT TABLE"
)

print(
    f"\nFrozen SNPs: {len(sites)}"
)

print(
    "Comparative samples retained: "
    f"{len(comparative_samples)}"
)

print(
    "Piano keys added: "
    + ",".join(
        piano_samples
    )
)

print(
    "\nFOCAL CALL COMPARISON"
)

print(
    focal_summary.to_string(
        index=False
    )
)

print(
    "\nPIANO-KEY CALLABILITY"
)

print(
    piano_summary.to_string(
        index=False
    )
)

print(
    f"\nLong table: {OUTPUT}"
)

print(
    f"Summary: {SUMMARY}"
)

print(
    f"Manifest: {MANIFEST}"
)
