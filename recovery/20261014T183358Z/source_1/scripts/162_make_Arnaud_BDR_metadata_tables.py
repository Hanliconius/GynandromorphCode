#!/usr/bin/env python3

from pathlib import Path
import re

import numpy as np
import pandas as pd


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

MANIFEST = (
    PROJECT
    / "analysis/figures/key_regions_population_context_genoplots"
    / "exact_BA_comparative_sample_manifest.tsv"
)

METADATA = (
    PROJECT
    / "data/derived/region_panel_metadata.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis/shareable_metadata_for_Arnaud"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)

SAMPLE_OUT = (
    OUTDIR
    / "Arnaud_BDR_exact_BA_comparative_metadata.tsv"
)

TAXON_OUT = (
    OUTDIR
    / "Arnaud_BDR_exact_BA_taxon_summary.tsv"
)

WIDE_OUT = (
    OUTDIR
    / "Arnaud_BDR_exact_BA_taxon_summary.wide.tsv"
)

OVERVIEW_OUT = (
    OUTDIR
    / "Arnaud_BDR_exact_BA_overview.tsv"
)


# ============================================================
# HELPERS
# ============================================================

def clean(x):
    if pd.isna(x):
        return ""

    return " ".join(
        str(x)
        .replace("\n", " ")
        .replace("\r", " ")
        .replace("\t", " ")
        .split()
    )


def canonicalize_taxon(x):
    """
    Match the intentional canonicalization used in analysis 155:

      ssp. 'meriana' -> ssp. meriana
      spp. melpomene -> ssp. melpomene
      spp. tristero  -> ssp. tristero
      nov. 'vicina'  -> nov. vicina

    plus whitespace cleanup.
    """

    x = clean(x)

    x = re.sub(
        r"\bspp\.",
        "ssp.",
        x,
        flags=re.IGNORECASE,
    )

    x = x.replace(
        "'",
        "",
    )

    return clean(x)


def norm_col(x):
    return re.sub(
        r"[^a-z0-9]+",
        "",
        clean(x).lower(),
    )


def find_col(df, names):
    lookup = {
        norm_col(col): col
        for col in df.columns
    }

    for name in names:
        key = norm_col(name)

        if key in lookup:
            return lookup[key]

    raise RuntimeError(
        "Could not find any of "
        + ", ".join(names)
        + "\nAvailable columns:\n"
        + "\n".join(df.columns)
    )


def binary(x):
    x = pd.to_numeric(
        pd.Series([x]),
        errors="coerce",
    ).iloc[0]

    if pd.isna(x):
        return np.nan

    if abs(float(x)) < 1e-8:
        return 0

    if abs(float(x) - 1.0) < 1e-8:
        return 1

    return np.nan


def count_state(x, state):
    return int(
        pd.to_numeric(
            x,
            errors="coerce",
        )
        .eq(state)
        .sum()
    )


def observed_states(frame):
    states = []

    for _, row in frame.iterrows():
        if (
            pd.notna(row["B"])
            and pd.notna(row["D"])
            and pd.notna(row["R"])
        ):
            states.append(
                "B{}D{}R{}".format(
                    int(row["B"]),
                    int(row["D"]),
                    int(row["R"]),
                )
            )

    if not states:
        return ""

    return ",".join(
        sorted(
            set(states)
        )
    )


# ============================================================
# READ EXACT BA MANIFEST
# ============================================================

for path in (
    MANIFEST,
    METADATA,
):
    if not path.is_file():
        raise RuntimeError(
            f"Missing required input: {path}"
        )


m = pd.read_csv(
    MANIFEST,
    sep="\t",
    low_memory=False,
)


required_manifest = {
    "BA_pattern",
    "sample",
    "phenotype",
    "taxon",
    "data_type",
}

missing = (
    required_manifest
    - set(m.columns)
)

if missing:
    raise RuntimeError(
        "Exact BA manifest missing: "
        + ", ".join(
            sorted(missing)
        )
    )


m["BA_pattern"] = m[
    "BA_pattern"
].map(clean)

m["sample"] = m[
    "sample"
].map(clean)

m["taxon_original"] = m[
    "taxon"
].map(clean)

m["taxon"] = m[
    "taxon"
].map(
    canonicalize_taxon
)

m["data_type"] = m[
    "data_type"
].map(clean)

m["phenotype"] = m[
    "phenotype"
].map(binary)


if m[
    [
        "BA_pattern",
        "sample",
    ]
].duplicated().any():

    bad = m.loc[
        m[
            [
                "BA_pattern",
                "sample",
            ]
        ].duplicated(
            keep=False
        )
    ]

    raise RuntimeError(
        "Duplicate pattern/sample rows:\n"
        + bad[
            [
                "BA_pattern",
                "sample",
                "taxon",
            ]
        ].to_string(
            index=False
        )
    )


# ============================================================
# LOCK EXACT COHORT COUNTS
# ============================================================

expected_n = {
    "Band": 163,
    "Dennis": 369,
    "Ray": 368,
}

observed_n = (
    m.groupby(
        "BA_pattern"
    )[
        "sample"
    ]
    .nunique()
    .to_dict()
)


for pattern, expected in expected_n.items():

    observed = observed_n.get(
        pattern,
        0,
    )

    if observed != expected:
        raise RuntimeError(
            f"{pattern}: expected "
            f"{expected} samples, "
            f"found {observed}."
        )


unexpected = (
    set(
        observed_n
    )
    - set(
        expected_n
    )
)

if unexpected:
    raise RuntimeError(
        "Unexpected BA patterns: "
        + ", ".join(
            sorted(unexpected)
        )
    )


# ============================================================
# READ MASTER B/D/R METADATA
# ============================================================

meta0 = pd.read_csv(
    METADATA,
    sep="\t",
    low_memory=False,
)


accession_col = find_col(
    meta0,
    ["Accession"],
)

taxon_col = find_col(
    meta0,
    ["Taxon Name"],
)

B_col = find_col(
    meta0,
    ["B"],
)

D_col = find_col(
    meta0,
    ["D"],
)

R_col = find_col(
    meta0,
    ["R"],
)


meta = pd.DataFrame(
    {
        "sample": (
            meta0[
                accession_col
            ]
            .map(clean)
        ),

        "metadata_taxon": (
            meta0[
                taxon_col
            ]
            .map(
                canonicalize_taxon
            )
        ),

        "B": (
            meta0[
                B_col
            ]
            .map(binary)
        ),

        "D": (
            meta0[
                D_col
            ]
            .map(binary)
        ),

        "R": (
            meta0[
                R_col
            ]
            .map(binary)
        ),
    }
)


meta = meta.loc[
    ~meta[
        "sample"
    ].eq("")
].copy()


# Duplicate accession rows are okay only if their B/D/R metadata agree.
duplicate_accessions = meta.loc[
    meta[
        "sample"
    ].duplicated(
        keep=False
    )
].copy()


if len(
    duplicate_accessions
):

    conflicts = []

    for sample, group in duplicate_accessions.groupby(
        "sample",
        sort=False,
    ):

        state_rows = (
            group[
                [
                    "B",
                    "D",
                    "R",
                ]
            ]
            .drop_duplicates()
        )

        if len(
            state_rows
        ) > 1:
            conflicts.append(
                sample
            )

    if conflicts:
        raise RuntimeError(
            "Conflicting duplicate phenotype metadata "
            "for accessions:\n"
            + "\n".join(
                conflicts
            )
        )


meta = meta.drop_duplicates(
    "sample"
)


# ============================================================
# JOIN: EXACT MEMBERSHIP + FULL B/D/R
# ============================================================

x = m.merge(
    meta,
    on="sample",
    how="left",
    validate="many_to_one",
)


missing_metadata = x[
    [
        "B",
        "D",
        "R",
    ]
].isna().all(
    axis=1
)


if missing_metadata.any():

    bad = x.loc[
        missing_metadata,
        [
            "BA_pattern",
            "sample",
            "taxon",
            "data_type",
        ],
    ]

    raise RuntimeError(
        "Exact-BA samples missing from "
        "region_panel_metadata.tsv:\n"
        + bad.to_string(
            index=False
        )
    )


# ============================================================
# VALIDATE THAT MANIFEST PHENOTYPE MATCHES THE CORRECT TRAIT
# ============================================================

trait_for_pattern = {
    "Band": "B",
    "Dennis": "D",
    "Ray": "R",
}

x["trait"] = x[
    "BA_pattern"
].map(
    trait_for_pattern
)


x["metadata_trait_phenotype"] = [
    row[
        trait_for_pattern[
            row["BA_pattern"]
        ]
    ]
    for _, row
    in x.iterrows()
]


comparable = (
    x[
        "phenotype"
    ].notna()
    & x[
        "metadata_trait_phenotype"
    ].notna()
)


discordant = x.loc[
    comparable
    & ~np.isclose(
        x[
            "phenotype"
        ].astype(float),
        x[
            "metadata_trait_phenotype"
        ].astype(float),
    )
]


if len(
    discordant
):

    raise RuntimeError(
        "Exact BA phenotype disagrees with "
        "master B/D/R metadata:\n"
        + discordant[
            [
                "BA_pattern",
                "sample",
                "taxon",
                "phenotype",
                "B",
                "D",
                "R",
            ]
        ].to_string(
            index=False
        )
    )


x["BDR_state"] = [
    (
        "B{}D{}R{}".format(
            int(row["B"]),
            int(row["D"]),
            int(row["R"]),
        )
        if (
            pd.notna(row["B"])
            and pd.notna(row["D"])
            and pd.notna(row["R"])
        )
        else ""
    )
    for _, row
    in x.iterrows()
]


x["phenotype_label"] = [
    (
        f"{trait}=0"
        if phenotype == 0
        else f"{trait}=1"
        if phenotype == 1
        else f"{trait}=?"
    )
    for trait, phenotype
    in zip(
        x["trait"],
        x["phenotype"],
    )
]


# ============================================================
# SHAREABLE SAMPLE-LEVEL TABLE
# ============================================================

sample_table = x[
    [
        "BA_pattern",
        "trait",
        "sample",
        "taxon",
        "data_type",
        "phenotype",
        "phenotype_label",
        "B",
        "D",
        "R",
        "BDR_state",
    ]
].copy()


sample_table = sample_table.sort_values(
    [
        "BA_pattern",
        "phenotype",
        "taxon",
        "sample",
    ]
)


sample_table.to_csv(
    SAMPLE_OUT,
    sep="\t",
    index=False,
)


# ============================================================
# TAXON-LEVEL LONG SUMMARY
# ============================================================

rows = []


for (
    pattern,
    taxon,
), group in x.groupby(
    [
        "BA_pattern",
        "taxon",
    ],
    sort=True,
    dropna=False,
):

    trait = trait_for_pattern[
        pattern
    ]

    phenotype0_n = (
        group[
            "phenotype"
        ]
        .eq(0)
        .sum()
    )

    phenotype1_n = (
        group[
            "phenotype"
        ]
        .eq(1)
        .sum()
    )

    n = len(
        group
    )

    rows.append(
        {
            "BA_pattern": pattern,
            "trait": trait,
            "taxon": taxon,
            "samples": n,

            "phenotype_0": int(
                phenotype0_n
            ),

            "phenotype_1": int(
                phenotype1_n
            ),

            "phenotype_1_fraction": (
                phenotype1_n / n
                if n
                else np.nan
            ),

            "WGS": int(
                group[
                    "data_type"
                ]
                .str.upper()
                .eq(
                    "WGS"
                )
                .sum()
            ),

            "CS": int(
                group[
                    "data_type"
                ]
                .str.upper()
                .eq(
                    "CS"
                )
                .sum()
            ),

            "BDR_states": (
                observed_states(
                    group
                )
            ),
        }
    )


taxon_summary = pd.DataFrame(
    rows
)


taxon_summary[
    "phenotype_1_fraction"
] = (
    taxon_summary[
        "phenotype_1_fraction"
    ]
    .round(4)
)


pattern_order = {
    "Band": 0,
    "Dennis": 1,
    "Ray": 2,
}


taxon_summary[
    "_pattern_order"
] = taxon_summary[
    "BA_pattern"
].map(
    pattern_order
)


taxon_summary = (
    taxon_summary
    .sort_values(
        [
            "_pattern_order",
            "phenotype_1_fraction",
            "taxon",
        ],
        ascending=[
            True,
            False,
            True,
        ],
    )
    .drop(
        columns=[
            "_pattern_order"
        ]
    )
)


taxon_summary.to_csv(
    TAXON_OUT,
    sep="\t",
    index=False,
)


# ============================================================
# WIDE TAXON TABLE FOR QUICK SHARING/SCANNING
# ============================================================

wide_rows = []


for taxon in sorted(
    taxon_summary[
        "taxon"
    ].unique()
):

    row = {
        "taxon": taxon,
    }

    subset = taxon_summary.loc[
        taxon_summary[
            "taxon"
        ].eq(
            taxon
        )
    ]

    for pattern in (
        "Band",
        "Dennis",
        "Ray",
    ):

        z = subset.loc[
            subset[
                "BA_pattern"
            ].eq(
                pattern
            )
        ]

        if len(
            z
        ) == 1:

            z = z.iloc[
                0
            ]

            row[
                f"{pattern}_samples"
            ] = int(
                z[
                    "samples"
                ]
            )

            row[
                f"{pattern}_phenotype0"
            ] = int(
                z[
                    "phenotype_0"
                ]
            )

            row[
                f"{pattern}_phenotype1"
            ] = int(
                z[
                    "phenotype_1"
                ]
            )

            row[
                f"{pattern}_presence_fraction"
            ] = z[
                "phenotype_1_fraction"
            ]

            row[
                f"{pattern}_BDR_states"
            ] = z[
                "BDR_states"
            ]

        else:

            row[
                f"{pattern}_samples"
            ] = ""

            row[
                f"{pattern}_phenotype0"
            ] = ""

            row[
                f"{pattern}_phenotype1"
            ] = ""

            row[
                f"{pattern}_presence_fraction"
            ] = ""

            row[
                f"{pattern}_BDR_states"
            ] = ""

    wide_rows.append(
        row
    )


wide = pd.DataFrame(
    wide_rows
)


wide.to_csv(
    WIDE_OUT,
    sep="\t",
    index=False,
)


# ============================================================
# OVERVIEW
# ============================================================

overview_rows = []


for pattern in (
    "Band",
    "Dennis",
    "Ray",
):

    group = x.loc[
        x[
            "BA_pattern"
        ].eq(
            pattern
        )
    ]

    trait = trait_for_pattern[
        pattern
    ]

    overview_rows.append(
        {
            "BA_pattern": pattern,
            "trait": trait,
            "samples": (
                group[
                    "sample"
                ].nunique()
            ),
            "taxa": (
                group[
                    "taxon"
                ].nunique()
            ),
            "phenotype_0": int(
                group[
                    "phenotype"
                ].eq(0)
                .sum()
            ),
            "phenotype_1": int(
                group[
                    "phenotype"
                ].eq(1)
                .sum()
            ),
        }
    )


overview = pd.DataFrame(
    overview_rows
)


overview.to_csv(
    OVERVIEW_OUT,
    sep="\t",
    index=False,
)


# ============================================================
# TERMINAL REPORT
# ============================================================

print(
    "\n========================================"
)

print(
    "ANALYSIS 162: SHAREABLE B/D/R METADATA"
)

print(
    "========================================\n"
)

print(
    "Exact cohort membership is taken only "
    "from exact_BA_comparative_sample_manifest.tsv."
)

print(
    "B/D/R phenotype columns are joined from "
    "region_panel_metadata.tsv by Accession.\n"
)


print(
    "OVERVIEW\n"
)

print(
    overview.to_string(
        index=False
    )
)


print(
    "\nTAXON SUMMARY\n"
)

print(
    taxon_summary.to_string(
        index=False
    )
)


print(
    "\nOUTPUTS"
)

for path in (
    SAMPLE_OUT,
    TAXON_OUT,
    WIDE_OUT,
    OVERVIEW_OUT,
):
    print(
        path
    )
