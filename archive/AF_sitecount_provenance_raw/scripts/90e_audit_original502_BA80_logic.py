#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import numpy as np


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

CK = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
)

PAIRED = (
    CK
    / "optix_original_502_paired_site_evidence.tsv"
)

EXISTING_BA80 = (
    CK
    / "optix_comparative_informative_sites_BA80.tsv"
)

OUT = (
    CK
    / "audit_original502_BA80_selection_logic.tsv"
)

OUT_BAND = (
    CK
    / "audit_original502_Band_candidates.tsv"
)


def truthy(x):
    return str(x).strip().lower() in {
        "true", "t", "1", "yes", "y"
    }


def as_num(x):
    return pd.to_numeric(
        x,
        errors="coerce",
    )


paired = pd.read_csv(
    PAIRED,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

existing = pd.read_csv(
    EXISTING_BA80,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

print("\nPAIRED TABLE")
print(f"rows: {len(paired)}")
print("\ncolumns:")
for i, col in enumerate(
    paired.columns,
    start=1,
):
    print(f"{i:3d}  {col}")


# ------------------------------------------------------------
# Trait / panel naming
# ------------------------------------------------------------

print("\nTRAIT COUNTS IN PAIRED TABLE")
print(
    paired.groupby(
        ["panel", "trait"],
        dropna=False,
    )
    .size()
    .to_string()
)

print("\nTRAIT COUNTS IN EXISTING BA80 TABLE")
print(
    existing.groupby(
        ["panel", "trait"],
        dropna=False,
    )
    .size()
    .to_string()
)


# ------------------------------------------------------------
# Convert the existing comparative flags.
# ------------------------------------------------------------

flag_columns = [
    "both_score_eligible",
    "both_q_significant",
    "both_missingness_ok",
    "comparative_stable",
    "comparative_high_BA80",
    "comparative_very_high_BA90",
    "matches_reference_expectation",
]

for col in flag_columns:
    if col in paired.columns:
        paired[col + "_bool"] = (
            paired[col].map(truthy)
        )


for col in [
    "raw_BA",
    "QC_BA",
    "minimum_BA",
    "raw_MCC",
    "QC_MCC",
    "minimum_MCC",
]:
    if col in paired.columns:
        paired[col + "_num"] = (
            as_num(
                paired[col]
            )
        )


# ------------------------------------------------------------
# Candidate filtering combinations.
#
# These are DIAGNOSTIC ONLY. Nothing gets selected here.
# ------------------------------------------------------------

conditions = {
    "minimum_BA_ge_0.80":
        paired[
            "minimum_BA_num"
        ] >= 0.80,
}

for col in flag_columns:
    bool_col = col + "_bool"

    if bool_col in paired.columns:
        conditions[col] = (
            paired[bool_col]
        )


if (
    "comparative_stable_bool"
    in paired.columns
):
    conditions[
        "stable_AND_BA80"
    ] = (
        paired[
            "comparative_stable_bool"
        ]
        &
        (
            paired[
                "minimum_BA_num"
            ] >= 0.80
        )
    )


if (
    "both_score_eligible_bool"
    in paired.columns
):
    conditions[
        "eligible_AND_BA80"
    ] = (
        paired[
            "both_score_eligible_bool"
        ]
        &
        (
            paired[
                "minimum_BA_num"
            ] >= 0.80
        )
    )


if all(
    x in paired.columns
    for x in [
        "both_score_eligible_bool",
        "comparative_stable_bool",
    ]
):
    conditions[
        "eligible_stable_BA80"
    ] = (
        paired[
            "both_score_eligible_bool"
        ]
        &
        paired[
            "comparative_stable_bool"
        ]
        &
        (
            paired[
                "minimum_BA_num"
            ] >= 0.80
        )
    )


if all(
    x in paired.columns
    for x in [
        "both_score_eligible_bool",
        "comparative_stable_bool",
        "both_missingness_ok_bool",
    ]
):
    conditions[
        "eligible_stable_missingness_BA80"
    ] = (
        paired[
            "both_score_eligible_bool"
        ]
        &
        paired[
            "comparative_stable_bool"
        ]
        &
        paired[
            "both_missingness_ok_bool"
        ]
        &
        (
            paired[
                "minimum_BA_num"
            ] >= 0.80
        )
    )


# ------------------------------------------------------------
# Summarise every rule by trait.
# ------------------------------------------------------------

summary_rows = []

for rule_name, mask in conditions.items():

    sub = paired.loc[
        mask.fillna(False)
    ].copy()

    grouped = (
        sub.groupby(
            ["panel", "trait"],
            dropna=False,
        )
        .size()
        .reset_index(
            name="n"
        )
    )

    for _, row in grouped.iterrows():
        summary_rows.append(
            {
                "rule":
                    rule_name,
                "panel":
                    row["panel"],
                "trait":
                    row["trait"],
                "n":
                    int(row["n"]),
            }
        )


summary = pd.DataFrame(
    summary_rows
)

summary.to_csv(
    OUT,
    sep="\t",
    index=False,
)


print("\nSELECTION-RULE COUNTS BY TRAIT")
if not summary.empty:
    print(
        summary.pivot_table(
            index=[
                "panel",
                "trait",
            ],
            columns="rule",
            values="n",
            fill_value=0,
            aggfunc="sum",
        ).to_string()
    )


# ------------------------------------------------------------
# Band audit.
#
# Print every Band row with any respectable comparative BA,
# including focal fields ONLY FOR LATER INSPECTION.
# They are not involved in any filter below.
# ------------------------------------------------------------

trait_text = (
    paired["trait"]
    .astype(str)
    .str.lower()
)

panel_text = (
    paired["panel"]
    .astype(str)
    .str.lower()
)

band = paired.loc[
    trait_text.str.contains(
        "band",
        regex=False,
    )
    |
    panel_text.str.contains(
        "band",
        regex=False,
    )
].copy()


band = band.loc[
    band[
        "minimum_BA_num"
    ].fillna(-np.inf) >= 0.70
].copy()


wanted = [
    "panel",
    "trait",
    "chrom",
    "pos",
    "ref",
    "alt",
    "selected_orientation",
    "stable_orientation",
    "matches_reference_expectation",
    "both_score_eligible",
    "both_q_significant",
    "both_missingness_ok",
    "comparative_stable",
    "comparative_high_BA80",
    "comparative_very_high_BA90",
    "raw_BA",
    "QC_BA",
    "minimum_BA",
    "raw_MCC",
    "QC_MCC",
    "minimum_MCC",
    "raw_errors",
    "QC_errors",
    "TN_Left_GT",
    "TN_Left_AD",
    "TN_Left_DP",
    "TN_Left_GQ",
    "TN_Right_GT",
    "TN_Right_AD",
    "TN_Right_DP",
    "TN_Right_GQ",
]

wanted = [
    col
    for col in wanted
    if col in band.columns
]


band = band[
    wanted
].copy()


if "minimum_BA" in band.columns:
    band["_sort_BA"] = (
        pd.to_numeric(
            band[
                "minimum_BA"
            ],
            errors="coerce",
        )
    )

    band = (
        band.sort_values(
            [
                "_sort_BA",
                "pos",
            ],
            ascending=[
                False,
                True,
            ],
        )
        .drop(
            columns="_sort_BA"
        )
    )


band.to_csv(
    OUT_BAND,
    sep="\t",
    index=False,
)


print("\nBAND ROWS WITH minimum_BA >= 0.70")
if band.empty:
    print("NONE")
else:
    print(
        band.to_string(
            index=False
        )
    )


# ------------------------------------------------------------
# Compare exact keys in the existing BA80 table against paired.
# ------------------------------------------------------------

keys = [
    "panel",
    "trait",
    "chrom",
    "pos",
    "ref",
    "alt",
]

keys = [
    x
    for x in keys
    if (
        x in paired.columns
        and x in existing.columns
    )
]


existing_keys = set(
    tuple(row)
    for row in existing[
        keys
    ].astype(str).itertuples(
        index=False,
        name=None,
    )
)


paired["_already_in_existing_BA80"] = [
    tuple(
        str(row[x])
        for x in keys
    )
    in existing_keys
    for _, row in paired.iterrows()
]


print(
    "\nEXISTING BA80 TABLE REPRODUCES "
    f"{paired['_already_in_existing_BA80'].sum()} "
    f"PAIRED ROWS (expected {len(existing)})"
)


print("\nOutputs:")
print(f"  {OUT}")
print(f"  {OUT_BAND}")
