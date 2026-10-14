#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import math
import re

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
    / "fixed_SNP_replication"
)

MANIFEST_OUT = (
    OUTDIR
    / "Wallbank2016_Ray_contigua_replication_samples.tsv"
)

SITE_OUT = (
    OUTDIR
    / "Wallbank2016_Ray_contigua_replication_sites.tsv"
)

WINDOW_OUT = (
    OUTDIR
    / "Wallbank2016_Ray_contigua_replication_5kb_windows.tsv"
)

SUMMARY_OUT = (
    OUTDIR
    / "Wallbank2016_Ray_contigua_replication_summary.tsv"
)


PUBLISHED_APPROXIMATE_PEAKS = {
    "Dennis": 49,
    "Ray": 36,
}

WINDOW_SIZE = 5000
WINDOW_STEP = 1000


def clean_column(value):

    return re.sub(
        r"\s+",
        " ",
        str(value),
    ).strip()


def clean_text(value):

    if value is None:
        return ""

    return re.sub(
        r"\s+",
        " ",
        str(value).replace(
            "'",
            "",
        ),
    ).strip().lower()


def as_trait(value):

    try:
        number = float(value)
    except (
        TypeError,
        ValueError,
    ):
        return None

    if number == 0:
        return 0

    if number == 1:
        return 1

    return None


def taxon_matches(
    taxon,
    terms,
):

    return any(
        term in taxon
        for term in terms
    )


def parse_gt(
    sample_field,
    format_lookup,
    scheme,
):

    values = sample_field.split(":")

    def value_for(name):

        index = format_lookup.get(name)

        if (
            index is None
            or index >= len(values)
        ):
            return ""

        return values[index]

    gt = value_for("GT")

    if (
        not gt
        or gt in {
            ".",
            "./.",
            ".|.",
        }
    ):
        return None

    alleles = re.split(
        r"[/|]",
        gt,
    )

    if (
        len(alleles) != 2
        or any(
            allele not in {
                "0",
                "1",
            }
            for allele in alleles
        )
    ):
        return None

    if scheme == "paper_QC":

        try:
            dp = float(
                value_for("DP")
            )
            gq = float(
                value_for("GQ")
            )
        except (
            TypeError,
            ValueError,
        ):
            return None

        if (
            dp < 4
            or dp > 300
            or gq < 30
        ):
            return None

    if alleles[0] != alleles[1]:
        return "heterozygous"

    return int(
        alleles[0]
    )


def evaluate_fixed_difference(
    calls_a,
    calls_b,
):

    total_a = len(calls_a)
    total_b = len(calls_b)

    callable_a = [
        call
        for call in calls_a
        if call is not None
    ]

    callable_b = [
        call
        for call in calls_b
        if call is not None
    ]

    hom_a = [
        call
        for call in callable_a
        if call in {
            0,
            1,
        }
    ]

    hom_b = [
        call
        for call in callable_b
        if call in {
            0,
            1,
        }
    ]

    heterozygote_present = (
        "heterozygous" in callable_a
        or "heterozygous" in callable_b
    )

    call_rate_a = (
        len(callable_a) / total_a
        if total_a
        else 0.0
    )

    call_rate_b = (
        len(callable_b) / total_b
        if total_b
        else 0.0
    )

    opposing_fixed_alleles = (
        not heterozygote_present
        and len(hom_a) >= 2
        and len(hom_b) >= 2
        and len(set(hom_a)) == 1
        and len(set(hom_b)) == 1
        and hom_a[0] != hom_b[0]
    )

    fixed_among_called = (
        opposing_fixed_alleles
    )

    fixed_80pct_callable = (
        opposing_fixed_alleles
        and call_rate_a >= 0.80
        and call_rate_b >= 0.80
    )

    fixed_complete = (
        opposing_fixed_alleles
        and len(callable_a) == total_a
        and len(callable_b) == total_b
    )

    return {
        "group_A_n": total_a,
        "group_B_n": total_b,
        "group_A_called": len(
            callable_a
        ),
        "group_B_called": len(
            callable_b
        ),
        "group_A_call_rate": call_rate_a,
        "group_B_call_rate": call_rate_b,
        "group_A_fixed_allele": (
            hom_a[0]
            if (
                hom_a
                and len(set(hom_a)) == 1
                and not heterozygote_present
            )
            else None
        ),
        "group_B_fixed_allele": (
            hom_b[0]
            if (
                hom_b
                and len(set(hom_b)) == 1
                and not heterozygote_present
            )
            else None
        ),
        "fixed_among_called_min2":
            fixed_among_called,
        "fixed_80pct_callable":
            fixed_80pct_callable,
        "fixed_complete":
            fixed_complete,
    }


if not VCF.is_file():
    raise RuntimeError(
        f"Missing broad Optix VCF: {VCF}"
    )

if not METADATA.is_file():
    raise RuntimeError(
        f"Missing metadata: {METADATA}"
    )

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)

metadata = pd.read_csv(
    METADATA,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

metadata.columns = [
    clean_column(column)
    for column in metadata.columns
]

required_columns = {
    "Accession",
    "Taxon Name",
    "D",
    "R",
}

missing_columns = (
    required_columns
    - set(metadata.columns)
)

if missing_columns:
    raise RuntimeError(
        "Metadata is missing columns: "
        + ", ".join(
            sorted(
                missing_columns
            )
        )
    )

metadata["taxon_clean"] = (
    metadata[
        "Taxon Name"
    ].map(
        clean_text
    )
)

metadata["D_value"] = (
    metadata["D"].map(
        as_trait
    )
)

metadata["R_value"] = (
    metadata["R"].map(
        as_trait
    )
)


# ============================================================
# Original Figure 2C combined-comparison taxa.
#
# In S1 Table the generic "with phenotype" and "without
# phenotype" headings are biologically reversed for dennis:
# contigua, cythera, plesseni, amaryllis and ordinary
# melpomene are no-D; malleti, ecuadoriensis, meriana and
# thelxiopeia are D-present.
# ============================================================

RAY_PRESENT_TAXA = (
    "timareta f. contigua",
    "ssp. malleti",
    "ssp. ecuadoriensis",
    "ssp. thelxiopeia",
)

RAY_ABSENT_TAXA = (
    "timareta f. timareta",
    "ssp. plesseni",
    "ssp. meriana",
)

DENNIS_ABSENT_TAXA = (
    "ssp. linaresi",
    "timareta f. contigua",
    "ssp. thelxinoe",
    "ssp. cythera",
    "ssp. plesseni",
    "ssp. amaryllis",
    "ssp. melpomene",
)

DENNIS_PRESENT_TAXA = (
    "ssp. florencia",
    "ssp. malleti",
    "ssp. ecuadoriensis",
    "ssp. aglaope",
    "ssp. meriana",
    "ssp. thelxiopeia",
)


def select_taxon_analogue(
    terms,
    trait_column,
    expected_trait,
):

    selected = metadata.loc[
        metadata[
            "taxon_clean"
        ].map(
            lambda value:
            taxon_matches(
                value,
                terms,
            )
        )
        & (
            metadata[
                trait_column
            ]
            == expected_trait
        )
    ].copy()

    return sorted(
        set(
            selected[
                "Accession"
            ]
        )
        - {""}
    )


# Expanded current-panel analogues of the original taxa.
dennis_absent_analogue = (
    select_taxon_analogue(
        DENNIS_ABSENT_TAXA,
        "D_value",
        0,
    )
)

dennis_present_analogue = (
    select_taxon_analogue(
        DENNIS_PRESENT_TAXA,
        "D_value",
        1,
    )
)

ray_absent_analogue = (
    select_taxon_analogue(
        RAY_ABSENT_TAXA,
        "R_value",
        0,
    )
)

ray_present_analogue = (
    select_taxon_analogue(
        RAY_PRESENT_TAXA,
        "R_value",
        1,
    )
)


# Exact retained biological specimens from Wallbank et al.
# after matching paper Sample Code to current EarthCapeID.
dennis_absent_overlap = [
    "ERR260302",  # JM-09-313, thelxinoe
    "ERR1143591", # CAM002856, cythera
    "ERR1143605", # CAM009156, plesseni
    "ERR260283",  # CAM009317, melpomene
]

dennis_present_overlap = [
    "ERR2298235", # CS002403, florencia
    "ERR3653533", # CAM013819, meriana
]


# ============================================================
# Hybrid Wallbank-ray replication.
#
# Presence group:
#   Current H. timareta f. contigua specimens with the
#   diagnostic D=0, R=1 phenotype.
#
# Absence group:
#   Original Wallbank et al. ray-absent biological specimens
#   that are represented in the current 502-sample VCF.
#
# The D=0, R=0 contigua record is excluded automatically.
# ============================================================

ray_present_current_contigua = sorted(
    set(
        metadata.loc[
            metadata[
                "taxon_clean"
            ].map(
                lambda value:
                taxon_matches(
                    value,
                    (
                        "timareta f. contigua",
                    ),
                )
            )
            & (
                metadata[
                    "D_value"
                ]
                == 0
            )
            & (
                metadata[
                    "R_value"
                ]
                == 1
            ),
            "Accession",
        ]
    )
    - {""}
)


ray_absent_original_overlap = [
    "ERR3653537", # CAM008533, H. t. timareta
    "ERR3653539", # CAM009178, H. t. timareta
    "ERR1143606", # CAM016293, H. m. plesseni
    "ERR1143605", # CAM009156, H. m. plesseni
    "ERR3653533", # CAM013819, H. m. meriana
]


analyses = {
    "Ray_current_contigua_vs_original_absent": {
        "trait": "Ray",
        "group_A_label": "original_no_R",
        "group_B_label": "current_contigua_D0R1",
        "group_A": ray_absent_original_overlap,
        "group_B": ray_present_current_contigua,
    },
}


# Read only the VCF header to establish column positions.
vcf_samples = []

with VCF.open(
    "rt",
    errors="replace",
) as handle:

    for line in handle:

        if line.startswith(
            "#CHROM"
        ):
            fields = (
                line.rstrip(
                    "\n"
                ).split(
                    "\t"
                )
            )
            vcf_samples = fields[9:]
            break

if not vcf_samples:
    raise RuntimeError(
        "Could not recover the VCF sample header."
    )

sample_to_column = {
    sample: index
    for index, sample
    in enumerate(
        vcf_samples
    )
}


manifest_rows = []

for analysis_name, specification in analyses.items():

    for side in (
        "group_A",
        "group_B",
    ):

        label = specification[
            f"{side}_label"
        ]

        retained = []

        for sample in specification[
            side
        ]:

            present = (
                sample
                in sample_to_column
            )

            manifest_rows.append(
                {
                    "analysis": analysis_name,
                    "trait":
                        specification["trait"],
                    "group": label,
                    "sample": sample,
                    "present_in_VCF":
                        present,
                }
            )

            if present:
                retained.append(
                    sample
                )

        specification[
            side
        ] = retained

    if (
        len(
            specification[
                "group_A"
            ]
        ) < 2
        or len(
            specification[
                "group_B"
            ]
        ) < 2
    ):
        raise RuntimeError(
            f"{analysis_name} has fewer than "
            "two VCF samples in one group."
        )

pd.DataFrame(
    manifest_rows
).to_csv(
    MANIFEST_OUT,
    sep="\t",
    index=False,
)


site_rows = []

schemes = (
    "called_GT",
    "paper_QC",
)

print(
    "Scanning existing broad VCF once."
)

with VCF.open(
    "rt",
    errors="replace",
) as handle:

    for line_number, line in enumerate(
        handle,
        start=1,
    ):

        if (
            not line
            or line.startswith(
                "#"
            )
        ):
            continue

        fields = (
            line.rstrip(
                "\n"
            ).split(
                "\t"
            )
        )

        if len(fields) < 10:
            continue

        chrom = fields[0]

        try:
            pos = int(
                fields[1]
            )
        except ValueError:
            continue

        ref = fields[3]
        alt = fields[4]

        if (
            len(ref) != 1
            or len(alt) != 1
            or "," in alt
        ):
            continue

        format_names = fields[
            8
        ].split(
            ":"
        )

        format_lookup = {
            name: index
            for index, name
            in enumerate(
                format_names
            )
        }

        for analysis_name, specification in analyses.items():

            group_a_columns = [
                sample_to_column[
                    sample
                ]
                for sample
                in specification[
                    "group_A"
                ]
            ]

            group_b_columns = [
                sample_to_column[
                    sample
                ]
                for sample
                in specification[
                    "group_B"
                ]
            ]

            for scheme in schemes:

                calls_a = [
                    parse_gt(
                        fields[
                            9 + column
                        ],
                        format_lookup,
                        scheme,
                    )
                    for column
                    in group_a_columns
                ]

                calls_b = [
                    parse_gt(
                        fields[
                            9 + column
                        ],
                        format_lookup,
                        scheme,
                    )
                    for column
                    in group_b_columns
                ]

                result = (
                    evaluate_fixed_difference(
                        calls_a,
                        calls_b,
                    )
                )

                site_rows.append(
                    {
                        "analysis":
                            analysis_name,
                        "trait":
                            specification[
                                "trait"
                            ],
                        "scheme": scheme,
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        **result,
                    }
                )

        if (
            line_number % 10000
            == 0
        ):
            print(
                f"Processed VCF line "
                f"{line_number:,}"
            )


sites = pd.DataFrame(
    site_rows
)

sites.to_csv(
    SITE_OUT,
    sep="\t",
    index=False,
)


# ============================================================
# Original publication used 5-kb windows advanced by 1 kb.
# ============================================================

window_rows = []

minimum_position = int(
    sites["pos"].min()
)

maximum_position = int(
    sites["pos"].max()
)

first_window = (
    minimum_position
    // WINDOW_STEP
    * WINDOW_STEP
)

for analysis_name in analyses:

    for scheme in schemes:

        subset = sites.loc[
            (
                sites["analysis"]
                == analysis_name
            )
            & (
                sites["scheme"]
                == scheme
            )
        ]

        for window_start in range(
            first_window,
            maximum_position + 1,
            WINDOW_STEP,
        ):

            window_end = (
                window_start
                + WINDOW_SIZE
                - 1
            )

            in_window = subset.loc[
                (
                    subset["pos"]
                    >= window_start
                )
                & (
                    subset["pos"]
                    <= window_end
                )
            ]

            window_rows.append(
                {
                    "analysis":
                        analysis_name,
                    "trait":
                        analyses[
                            analysis_name
                        ][
                            "trait"
                        ],
                    "scheme": scheme,
                    "window_start":
                        window_start,
                    "window_end":
                        window_end,
                    "window_midpoint":
                        (
                            window_start
                            + WINDOW_SIZE / 2
                        ),
                    "fixed_among_called_min2":
                        int(
                            in_window[
                                "fixed_among_called_min2"
                            ].sum()
                        ),
                    "fixed_80pct_callable":
                        int(
                            in_window[
                                "fixed_80pct_callable"
                            ].sum()
                        ),
                    "fixed_complete":
                        int(
                            in_window[
                                "fixed_complete"
                            ].sum()
                        ),
                }
            )


windows = pd.DataFrame(
    window_rows
)

windows.to_csv(
    WINDOW_OUT,
    sep="\t",
    index=False,
)


summary_rows = []

for analysis_name, specification in analyses.items():

    for scheme in schemes:

        site_subset = sites.loc[
            (
                sites["analysis"]
                == analysis_name
            )
            & (
                sites["scheme"]
                == scheme
            )
        ]

        window_subset = windows.loc[
            (
                windows["analysis"]
                == analysis_name
            )
            & (
                windows["scheme"]
                == scheme
            )
        ]

        row = {
            "analysis": analysis_name,
            "trait":
                specification["trait"],
            "scheme": scheme,
            "group_A_label":
                specification[
                    "group_A_label"
                ],
            "group_A_n":
                len(
                    specification[
                        "group_A"
                    ]
                ),
            "group_B_label":
                specification[
                    "group_B_label"
                ],
            "group_B_n":
                len(
                    specification[
                        "group_B"
                    ]
                ),
            "biallelic_SNPs":
                len(
                    site_subset
                ),
            "total_fixed_among_called_min2":
                int(
                    site_subset[
                        "fixed_among_called_min2"
                    ].sum()
                ),
            "total_fixed_80pct_callable":
                int(
                    site_subset[
                        "fixed_80pct_callable"
                    ].sum()
                ),
            "total_fixed_complete":
                int(
                    site_subset[
                        "fixed_complete"
                    ].sum()
                ),
            "maximum_5kb_fixed_among_called_min2":
                int(
                    window_subset[
                        "fixed_among_called_min2"
                    ].max()
                ),
            "maximum_5kb_fixed_80pct_callable":
                int(
                    window_subset[
                        "fixed_80pct_callable"
                    ].max()
                ),
            "maximum_5kb_fixed_complete":
                int(
                    window_subset[
                        "fixed_complete"
                    ].max()
                ),
            "published_approximate_5kb_peak":
                PUBLISHED_APPROXIMATE_PEAKS[
                    specification[
                        "trait"
                    ]
                ],
        }

        best_index = (
            window_subset[
                "fixed_80pct_callable"
            ].idxmax()
        )

        best_window = windows.loc[
            best_index
        ]

        row[
            "best_80pct_window_start"
        ] = int(
            best_window[
                "window_start"
            ]
        )

        row[
            "best_80pct_window_end"
        ] = int(
            best_window[
                "window_end"
            ]
        )

        summary_rows.append(
            row
        )


summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    SUMMARY_OUT,
    sep="\t",
    index=False,
)


print()
print(
    "WALLBANK ET AL. 2016 FIXED-SNP REPLICATION"
)
print()
print(
    "Original analysis: 5-kb windows, 1-kb step."
)
print(
    "Paper-QC calls require DP 4-300 and GQ >=30."
)
print(
    "Mapping quality cannot be reproduced per genotype "
    "from the present VCF FORMAT fields."
)
print()
print(
    summary.to_string(
        index=False,
    )
)
print()
print(
    f"Sample manifest: {MANIFEST_OUT}"
)
print(
    f"Per-site results: {SITE_OUT}"
)
print(
    f"Sliding windows: {WINDOW_OUT}"
)
print(
    f"Summary: {SUMMARY_OUT}"
)
