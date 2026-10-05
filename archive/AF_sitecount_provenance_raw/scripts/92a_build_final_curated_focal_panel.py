#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

CK = PROJECT / "analysis" / "checkpoints"

CK502 = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
)

EXPECTATIONS = (
    CK502
    / "focal_module_expectations.tsv"
)

# ============================================================
# OPTIX
#
# Site selection: comparative 51-sample candidate set only.
# Focal data: previously constructed pseudopools from the
# four focal legs. Focal data do NOT determine site selection.
# ============================================================

OPTIX_SELECTED = (
    CK
    / "optix_51sample_candidate_pattern_groups.tsv"
)

OPTIX_FOCAL = (
    CK
    / "optix_focal_pseudopool_evidence.long.tsv"
)

# ============================================================
# Yb
#
# Site selection: comparative-only BA>=0.80 table.
# Focal data: pooled TN_Left/TN_Right regional VCF.
# ============================================================

YB_SELECTED = (
    CK502
    / "ivory_comparative_informative_sites_BA80.tsv"
)

YB_VCF = (
    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
    / "Ivory_Hmr15-1627000-1637000_biallelic.vcf"
)

OUT_SITES = (
    CK502
    / "final_curated_comparative_selected_sites.tsv"
)

OUT_LONG = (
    CK502
    / "final_curated_focal_biological_evidence.long.tsv"
)

OUT_WIDE = (
    CK502
    / "final_curated_focal_biological_evidence.wide.tsv"
)

OUT_SUMMARY = (
    CK502
    / "final_curated_focal_biological_evidence.summary.tsv"
)


MODULE_ORDER = [
    "Band",
    "Ray1",
    "Ray2",
    "Dennis1",
    "Dennis2",
    "Yb",
]


DISPLAY = {
    "Band": "Band",
    "Ray1": "Ray 1",
    "Ray2": "Ray 2",
    "Dennis1": "Dennis 1",
    "Dennis2": "Dennis 2",
    "Yb": "Yb",
}


REFERENCE_STATE = {
    "Band": "presence",
    "Ray1": "absence",
    "Ray2": "absence",
    "Dennis1": "absence",
    "Dennis2": "absence",
    "Yb": "presence",
}


PANEL_TO_MODULE = {
    "band": "Band",
    "ray1": "Ray1",
    "ray2": "Ray2",
    "dennis1": "Dennis1",
    "dennis2": "Dennis2",
}


def fail(message):
    raise RuntimeError(message)


def read_tsv(path):
    if not path.exists():
        fail(
            f"Missing required file: {path}"
        )

    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )


def emptyish(value):
    if value is None:
        return True

    return str(value).strip() in {
        "",
        ".",
        "NA",
        "NaN",
        "nan",
        "None",
    }


def fnum(value):
    if emptyish(value):
        return np.nan

    try:
        return float(value)
    except Exception:
        return np.nan


def truthy(value):
    return (
        str(value)
        .strip()
        .lower()
        in {
            "1",
            "true",
            "t",
            "yes",
            "y",
        }
    )


def norm_allele(value):
    if emptyish(value):
        return ""

    return (
        str(value)
        .strip()
        .upper()
    )


def parse_gt(gt):
    if emptyish(gt):
        return None

    parts = re.split(
        r"[/|]",
        str(gt).strip(),
    )

    if len(parts) != 2:
        return None

    if any(
        allele not in {"0", "1"}
        for allele in parts
    ):
        return None

    return sum(
        int(allele)
        for allele in parts
    )


def genotype_from_raw(
    module,
    raw_gt,
):
    """
    Convert VCF GT directly into biological
    phenotype-presence space.

    Reference biological states are locked:

      Band:     reference = presence
      Yb:       reference = presence
      Dennis/R: reference = absence
    """

    alt_copies = parse_gt(
        raw_gt
    )

    if alt_copies is None:
        return None

    if (
        REFERENCE_STATE[module]
        == "presence"
    ):
        presence_copies = (
            2 - alt_copies
        )
    else:
        presence_copies = (
            alt_copies
        )

    return {
        0: "AA",
        1: "PA",
        2: "PP",
    }[
        presence_copies
    ]


def biological_dosage(gt):
    return {
        "AA": 0.0,
        "PA": 0.5,
        "PP": 1.0,
    }.get(
        gt,
        np.nan,
    )


def presence_fraction(
    module,
    ref_depth,
    alt_depth,
):
    ref_depth = fnum(
        ref_depth
    )

    alt_depth = fnum(
        alt_depth
    )

    if (
        not np.isfinite(ref_depth)
        or not np.isfinite(alt_depth)
        or (
            ref_depth
            + alt_depth
        ) <= 0
    ):
        return np.nan

    if (
        REFERENCE_STATE[module]
        == "presence"
    ):
        return (
            ref_depth
            / (
                ref_depth
                + alt_depth
            )
        )

    return (
        alt_depth
        / (
            ref_depth
            + alt_depth
        )
    )


def choose_column(
    row,
    possibilities,
    default="",
):
    for column in possibilities:
        if column in row.index:
            value = row[column]

            if not emptyish(value):
                return value

    return default


def module_from_panel(
    value,
):
    text = (
        str(value)
        .strip()
        .lower()
    )

    compact = re.sub(
        r"[^a-z0-9]",
        "",
        text,
    )

    for key, module in (
        PANEL_TO_MODULE.items()
    ):
        if key in compact:
            return module

    return None


def evaluate_call(
    module,
    gt,
    ref_depth,
    alt_depth,
    dp,
    gq,
    expected_gt,
):
    biological_gt = (
        genotype_from_raw(
            module,
            gt,
        )
    )

    dp_num = fnum(
        dp
    )

    gq_num = fnum(
        gq
    )

    qc = bool(
        biological_gt is not None
        and np.isfinite(
            dp_num
        )
        and dp_num >= 2
        and np.isfinite(
            gq_num
        )
        and gq_num >= 20
    )

    if qc:
        dosage = (
            biological_dosage(
                biological_gt
            )
        )

        p_fraction = (
            presence_fraction(
                module,
                ref_depth,
                alt_depth,
            )
        )

        matches = (
            biological_gt
            == expected_gt
        )

    else:
        dosage = np.nan
        p_fraction = np.nan
        matches = np.nan

    return {
        "GT_raw":
            gt,

        "ref_depth":
            fnum(
                ref_depth
            ),

        "alt_depth":
            fnum(
                alt_depth
            ),

        "DP":
            dp_num,

        "GQ":
            gq_num,

        "QC_pass_DP2_GQ20":
            qc,

        "biological_genotype":
            (
                biological_gt
                if qc
                else "missing"
            ),

        "presence_dosage":
            dosage,

        "presence_allele_fraction":
            p_fraction,

        "expected_genotype":
            expected_gt,

        "matches_expectation":
            matches,
    }


# ============================================================
# Yb VCF parser
# ============================================================

def read_yb_vcf(
    selected_positions,
):
    samples = None
    records = {}

    with YB_VCF.open(
        "rt",
        errors="replace",
    ) as handle:

        for line in handle:

            if line.startswith("##"):
                continue

            if line.startswith(
                "#CHROM\t"
            ):
                fields = (
                    line.rstrip("\n")
                    .split("\t")
                )

                samples = fields[9:]
                continue

            if line.startswith("#"):
                continue

            fields = (
                line.rstrip("\n")
                .split("\t")
            )

            if len(fields) < 10:
                continue

            chrom = fields[0]

            try:
                pos = int(
                    fields[1]
                )
            except Exception:
                continue

            if pos not in selected_positions:
                continue

            ref = norm_allele(
                fields[3]
            )

            alt = norm_allele(
                fields[4]
            )

            fmt_names = (
                fields[8]
                .split(":")
            )

            calls = {}

            for sample, text in zip(
                samples,
                fields[9:],
            ):
                values = (
                    text.split(":")
                )

                calls[sample] = {
                    name:
                        (
                            values[i]
                            if i
                            < len(values)
                            else ""
                        )
                    for i, name
                    in enumerate(
                        fmt_names
                    )
                }

            records[
                (
                    chrom,
                    pos,
                    ref,
                    alt,
                )
            ] = calls

    if samples is None:
        fail(
            f"Could not parse VCF header: {YB_VCF}"
        )

    required = {
        "TN_Left",
        "TN_Right",
    }

    if not required.issubset(
        samples
    ):
        fail(
            "Yb VCF does not contain "
            "TN_Left and TN_Right."
        )

    return records


def main():

    expectations = read_tsv(
        EXPECTATIONS
    )

    exp = (
        expectations
        .set_index("module")
        .to_dict(
            orient="index"
        )
    )

    long_rows = []
    site_rows = []

    # ========================================================
    # OPTIX SELECTED SITES
    # ========================================================

    selected = read_tsv(
        OPTIX_SELECTED
    )

    focal = read_tsv(
        OPTIX_FOCAL
    )

    selected["module"] = [
        module_from_panel(
            row.get(
                "panel",
                "",
            )
        )
        for _, row
        in selected.iterrows()
    ]

    selected = selected.loc[
        selected[
            "module"
        ].isin(
            MODULE_ORDER
        )
    ].copy()

    focal["module"] = [
        module_from_panel(
            row.get(
                "panel",
                "",
            )
        )
        for _, row
        in focal.iterrows()
    ]

    selected[
        "pos_num"
    ] = pd.to_numeric(
        selected["pos"],
        errors="raise",
    ).astype(int)

    focal[
        "pos_num"
    ] = pd.to_numeric(
        focal["pos"],
        errors="raise",
    ).astype(int)

    print(
        "\nOPTIX CURATED COMPARATIVE SITE COUNTS"
    )

    print(
        selected.groupby(
            "module"
        )
        .size()
        .reindex(
            MODULE_ORDER[:-1],
            fill_value=0,
        )
        .to_string()
    )

    print(
        "\nTotal Optix selected:",
        len(selected),
    )

    # Each selected site stays in the analysis even if
    # focal evidence is missing.
    for _, site in (
        selected.iterrows()
    ):

        module = site[
            "module"
        ]

        chrom = site[
            "chrom"
        ]

        pos = int(
            site[
                "pos_num"
            ]
        )

        ref = norm_allele(
            site["ref"]
        )

        alt = norm_allele(
            site["alt"]
        )

        site_id = (
            f"{module}:"
            f"{chrom}:"
            f"{pos}:"
            f"{ref}>{alt}"
        )

        minimum_ba = fnum(
            choose_column(
                site,
                [
                    "minimum_BA",
                    "minimum_BA_51sample",
                    "balanced_accuracy",
                ],
            )
        )

        minimum_mcc = fnum(
            choose_column(
                site,
                [
                    "minimum_MCC",
                    "minimum_MCC_51sample",
                    "MCC",
                ],
            )
        )

        tier = choose_column(
            site,
            [
                "tier",
                "candidate_tier",
            ],
        )

        group_id = (
            choose_column(
                site,
                [
                    "group_id",
                    "pattern_group",
                ],
            )
        )

        is_representative = (
            choose_column(
                site,
                [
                    "is_representative",
                    "representative",
                ],
            )
        )

        site_base = {
            "site_id":
                site_id,

            "module":
                module,

            "module_display":
                DISPLAY[module],

            "chrom":
                chrom,

            "pos":
                pos,

            "ref":
                ref,

            "alt":
                alt,

            "minimum_comparative_BA":
                minimum_ba,

            "minimum_comparative_MCC":
                minimum_mcc,

            "tier":
                tier,

            "group_id":
                group_id,

            "is_representative":
                is_representative,

            "ascertainment_panel":
                "curated_51sample_Optix",

            "ascertainment_rule":
                "predefined comparative candidate-pattern set",

            "reference_state":
                REFERENCE_STATE[
                    module
                ],
        }

        site_rows.append(
            site_base
        )

        focal_match = focal.loc[
            (
                focal[
                    "module"
                ]
                == module
            )
            &
            (
                focal[
                    "chrom"
                ]
                == chrom
            )
            &
            (
                focal[
                    "pos_num"
                ]
                == pos
            )
            &
            (
                focal[
                    "ref"
                ]
                .str.upper()
                == ref
            )
            &
            (
                focal[
                    "alt"
                ]
                .str.upper()
                == alt
            )
        ].copy()

        for side in [
            "Left",
            "Right",
        ]:

            expected_gt = (
                exp[module][
                    (
                        "left_expected_genotype"
                        if side
                        == "Left"
                        else
                        "right_expected_genotype"
                    )
                ]
            )

            side_match = (
                focal_match.loc[
                    focal_match[
                        "side"
                    ]
                    .astype(str)
                    .str.lower()
                    == side.lower()
                ]
            )

            if side_match.empty:

                result = {
                    "GT_raw":
                        "missing",
                    "ref_depth":
                        np.nan,
                    "alt_depth":
                        np.nan,
                    "DP":
                        np.nan,
                    "GQ":
                        np.nan,
                    "QC_pass_DP2_GQ20":
                        False,
                    "biological_genotype":
                        "missing",
                    "presence_dosage":
                        np.nan,
                    "presence_allele_fraction":
                        np.nan,
                    "expected_genotype":
                        expected_gt,
                    "matches_expectation":
                        np.nan,
                }

                focal_status = (
                    "no_pseudopool_row"
                )

            else:

                row = (
                    side_match
                    .iloc[0]
                )

                result = evaluate_call(
                    module=module,

                    gt=row.get(
                        "pseudo_GT",
                        "",
                    ),

                    ref_depth=row.get(
                        "combined_ref_depth",
                        "",
                    ),

                    alt_depth=row.get(
                        "combined_alt_depth",
                        "",
                    ),

                    dp=row.get(
                        "combined_DP",
                        "",
                    ),

                    gq=row.get(
                        "pseudo_GQ",
                        "",
                    ),

                    expected_gt=
                        expected_gt,
                )

                focal_status = (
                    "pseudopool"
                )

            long_rows.append(
                {
                    **site_base,

                    "side":
                        side,

                    "focal_source":
                        "four-leg pseudopool",

                    "focal_status":
                        focal_status,

                    **result,
                }
            )

    # ========================================================
    # Yb SELECTED SITES
    # ========================================================

    yb = read_tsv(
        YB_SELECTED
    )

    yb[
        "pos_num"
    ] = pd.to_numeric(
        yb["pos"],
        errors="raise",
    ).astype(int)

    # Unique physical SNPs.
    yb = (
        yb.sort_values(
            "pos_num"
        )
        .drop_duplicates(
            [
                "chrom",
                "pos_num",
                "ref",
                "alt",
            ]
        )
        .copy()
    )

    yb_positions = set(
        yb[
            "pos_num"
        ].astype(int)
    )

    yb_records = (
        read_yb_vcf(
            yb_positions
        )
    )

    for _, site in yb.iterrows():

        module = "Yb"

        chrom = site[
            "chrom"
        ]

        pos = int(
            site[
                "pos_num"
            ]
        )

        ref = norm_allele(
            site["ref"]
        )

        alt = norm_allele(
            site["alt"]
        )

        site_id = (
            f"Yb:{chrom}:"
            f"{pos}:{ref}>{alt}"
        )

        site_base = {
            "site_id":
                site_id,

            "module":
                "Yb",

            "module_display":
                "Yb",

            "chrom":
                chrom,

            "pos":
                pos,

            "ref":
                ref,

            "alt":
                alt,

            "minimum_comparative_BA":
                fnum(
                    site.get(
                        "minimum_BA",
                        "",
                    )
                ),

            "minimum_comparative_MCC":
                fnum(
                    site.get(
                        "minimum_MCC",
                        "",
                    )
                ),

            "tier":
                "",

            "group_id":
                "",

            "is_representative":
                "",

            "ascertainment_panel":
                "comparative_Yb",

            "ascertainment_rule":
                "comparative minimum BA >= 0.80",

            "reference_state":
                "presence",
        }

        site_rows.append(
            site_base
        )

        calls = (
            yb_records.get(
                (
                    chrom,
                    pos,
                    ref,
                    alt,
                ),
                {},
            )
        )

        for side in [
            "Left",
            "Right",
        ]:

            sample = (
                "TN_Left"
                if side == "Left"
                else "TN_Right"
            )

            expected_gt = (
                exp["Yb"][
                    (
                        "left_expected_genotype"
                        if side
                        == "Left"
                        else
                        "right_expected_genotype"
                    )
                ]
            )

            call = calls.get(
                sample
            )

            if call is None:

                result = {
                    "GT_raw":
                        "missing",
                    "ref_depth":
                        np.nan,
                    "alt_depth":
                        np.nan,
                    "DP":
                        np.nan,
                    "GQ":
                        np.nan,
                    "QC_pass_DP2_GQ20":
                        False,
                    "biological_genotype":
                        "missing",
                    "presence_dosage":
                        np.nan,
                    "presence_allele_fraction":
                        np.nan,
                    "expected_genotype":
                        expected_gt,
                    "matches_expectation":
                        np.nan,
                }

                focal_status = (
                    "site_missing_from_Yb_VCF"
                )

            else:

                ad = (
                    str(
                        call.get(
                            "AD",
                            "",
                        )
                    )
                    .split(",")
                )

                ref_depth = (
                    ad[0]
                    if len(ad) >= 1
                    else ""
                )

                alt_depth = (
                    ad[1]
                    if len(ad) >= 2
                    else ""
                )

                result = evaluate_call(
                    module="Yb",

                    gt=call.get(
                        "GT",
                        "",
                    ),

                    ref_depth=
                        ref_depth,

                    alt_depth=
                        alt_depth,

                    dp=call.get(
                        "DP",
                        "",
                    ),

                    gq=call.get(
                        "GQ",
                        "",
                    ),

                    expected_gt=
                        expected_gt,
                )

                focal_status = (
                    "pooled_side_VCF"
                )

            long_rows.append(
                {
                    **site_base,

                    "side":
                        side,

                    "focal_source":
                        "TN_Left/TN_Right Yb VCF",

                    "focal_status":
                        focal_status,

                    **result,
                }
            )

    # ========================================================
    # FINAL TABLES
    # ========================================================

    sites = pd.DataFrame(
        site_rows
    )

    long = pd.DataFrame(
        long_rows
    )

    module_rank = {
        module: i
        for i, module
        in enumerate(
            MODULE_ORDER
        )
    }

    sites[
        "_module_rank"
    ] = sites[
        "module"
    ].map(
        module_rank
    )

    sites = (
        sites.sort_values(
            [
                "_module_rank",
                "pos",
            ]
        )
        .drop(
            columns="_module_rank"
        )
        .reset_index(
            drop=True
        )
    )

    sites[
        "site_order_within_module"
    ] = (
        sites.groupby(
            "module"
        )
        .cumcount()
        + 1
    )

    order_lookup = (
        sites.set_index(
            "site_id"
        )[
            "site_order_within_module"
        ]
        .to_dict()
    )

    long[
        "site_order_within_module"
    ] = (
        long[
            "site_id"
        ].map(
            order_lookup
        )
    )

    long[
        "_module_rank"
    ] = (
        long[
            "module"
        ].map(
            module_rank
        )
    )

    long = (
        long.sort_values(
            [
                "_module_rank",
                "pos",
                "side",
            ]
        )
        .drop(
            columns="_module_rank"
        )
        .reset_index(
            drop=True
        )
    )

    base_columns = [
        "site_id",
        "module",
        "module_display",
        "chrom",
        "pos",
        "ref",
        "alt",
        "minimum_comparative_BA",
        "minimum_comparative_MCC",
        "tier",
        "group_id",
        "is_representative",
        "ascertainment_panel",
        "ascertainment_rule",
        "reference_state",
        "site_order_within_module",
    ]

    value_columns = [
        "focal_source",
        "focal_status",
        "GT_raw",
        "ref_depth",
        "alt_depth",
        "DP",
        "GQ",
        "QC_pass_DP2_GQ20",
        "biological_genotype",
        "presence_dosage",
        "presence_allele_fraction",
        "expected_genotype",
        "matches_expectation",
    ]

    left = (
        long.loc[
            long["side"]
            == "Left",
            base_columns
            + value_columns,
        ]
        .copy()
        .rename(
            columns={
                column:
                    f"left_{column}"
                for column
                in value_columns
            }
        )
    )

    right = (
        long.loc[
            long["side"]
            == "Right",
            base_columns
            + value_columns,
        ]
        .copy()
        .rename(
            columns={
                column:
                    f"right_{column}"
                for column
                in value_columns
            }
        )
    )

    wide = left.merge(
        right,
        on=base_columns,
        how="outer",
        validate="one_to_one",
    )

    wide[
        "both_callable"
    ] = (
        wide[
            "left_QC_pass_DP2_GQ20"
        ].astype(bool)
        &
        wide[
            "right_QC_pass_DP2_GQ20"
        ].astype(bool)
    )

    wide[
        "pair_matches_expectation"
    ] = np.where(
        wide[
            "both_callable"
        ],
        (
            wide[
                "left_matches_expectation"
            ].astype(bool)
            &
            wide[
                "right_matches_expectation"
            ].astype(bool)
        ),
        np.nan,
    )

    # ========================================================
    # SUMMARY
    # ========================================================

    summary_rows = []

    for module in MODULE_ORDER:

        d = wide.loc[
            wide[
                "module"
            ]
            == module
        ].copy()

        if d.empty:
            continue

        left_callable = (
            d[
                "left_QC_pass_DP2_GQ20"
            ].astype(bool)
        )

        right_callable = (
            d[
                "right_QC_pass_DP2_GQ20"
            ].astype(bool)
        )

        both = (
            d[
                "both_callable"
            ].astype(bool)
        )

        paired = d.loc[
            both
        ]

        summary_rows.append(
            {
                "module":
                    module,

                "selected_sites":
                    len(d),

                "left_callable":
                    int(
                        left_callable.sum()
                    ),

                "right_callable":
                    int(
                        right_callable.sum()
                    ),

                "both_callable":
                    int(
                        both.sum()
                    ),

                "left_exact_expected_n":
                    int(
                        d.loc[
                            left_callable,
                            "left_matches_expectation",
                        ]
                        .fillna(False)
                        .astype(bool)
                        .sum()
                    ),

                "right_exact_expected_n":
                    int(
                        d.loc[
                            right_callable,
                            "right_matches_expectation",
                        ]
                        .fillna(False)
                        .astype(bool)
                        .sum()
                    ),

                "pair_exact_expected_n":
                    int(
                        paired[
                            "pair_matches_expectation"
                        ]
                        .fillna(False)
                        .astype(bool)
                        .sum()
                    ),

                "pair_exact_expected_fraction":
                    (
                        paired[
                            "pair_matches_expectation"
                        ]
                        .astype(float)
                        .mean()
                        if len(
                            paired
                        )
                        else np.nan
                    ),

                "left_expected_genotype":
                    exp[module][
                        "left_expected_genotype"
                    ],

                "right_expected_genotype":
                    exp[module][
                        "right_expected_genotype"
                    ],

                "left_median_presence_fraction":
                    pd.to_numeric(
                        d[
                            "left_presence_allele_fraction"
                        ],
                        errors="coerce",
                    ).median(),

                "right_median_presence_fraction":
                    pd.to_numeric(
                        d[
                            "right_presence_allele_fraction"
                        ],
                        errors="coerce",
                    ).median(),
            }
        )

    summary = pd.DataFrame(
        summary_rows
    )

    sites.to_csv(
        OUT_SITES,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    long.to_csv(
        OUT_LONG,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    wide.to_csv(
        OUT_WIDE,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary.to_csv(
        OUT_SUMMARY,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print(
        "\nFINAL SELECTED SITE COUNTS"
    )

    print(
        sites.groupby(
            "module"
        )
        .size()
        .reindex(
            MODULE_ORDER,
            fill_value=0,
        )
        .to_string()
    )

    print(
        "\nFINAL FOCAL SUMMARY"
    )

    print(
        summary.to_string(
            index=False
        )
    )

    print(
        "\nFOCAL STATUS"
    )

    print(
        long.groupby(
            [
                "module",
                "focal_status",
            ]
        )
        .size()
        .to_string()
    )

    print(
        "\nOutputs:"
    )

    for path in [
        OUT_SITES,
        OUT_LONG,
        OUT_WIDE,
        OUT_SUMMARY,
    ]:
        print(
            f"  {path}"
        )


if __name__ == "__main__":
    main()
