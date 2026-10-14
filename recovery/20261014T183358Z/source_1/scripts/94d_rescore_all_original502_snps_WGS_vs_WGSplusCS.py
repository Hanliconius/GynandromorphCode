#!/usr/bin/env python3

from pathlib import Path
import csv
import importlib.util
import re

import numpy as np
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

ORIGINAL_SCRIPT = (
    PROJECT
    / "scripts"
    / "41_score_original_502_optix_callsets.py"
)

META = (
    PROJECT
    / "data"
    / "derived"
    / "region_panel_metadata.tsv"
)

FILTER_MANIFEST = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
)

OUT_RAW = (
    OUTDIR
    / "original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.tsv"
)

OUT_PRIMARY = (
    OUTDIR
    / "original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.DP2_GQ20.tsv"
)

OUT_PAIRED = (
    OUTDIR
    / "original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.minimum_across_schemes.tsv"
)

OUT_SUMMARY = (
    OUTDIR
    / "original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.summary.tsv"
)

OUT_COHORT_SUMMARY = (
    OUTDIR
    / "original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.cohorts.tsv"
)

OUT_SCORE_MANIFEST = (
    OUTDIR
    / "original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.score_manifest.tsv"
)


MODULES = [
    "Band",
    "Ray2",
    "Dennis2",
]

COHORTS = [
    "WGS_only",
    "WGS_plus_CS",
]

TRAIT_FROM_MODULE = {
    "Band": "B",
    "Ray2": "R",
    "Dennis2": "D",
}

PRIMARY_SCHEME = "DP2_GQ20"

EXCLUDED_TAXON_TERMS = (
    "cydno",
    "numata",
    "pachinus",
    "hecale",
    "ismenius",
    "contigua",
    "linaresi",
)


def fail(message):
    raise RuntimeError(message)


def as_bool(value):

    return str(value).strip().lower() in {
        "true",
        "t",
        "1",
        "yes",
        "y",
    }


def clean(value):

    if value is None:
        return ""

    return str(value).strip()


def load_original_module():

    if not ORIGINAL_SCRIPT.exists():
        fail(
            f"Missing original scorer: {ORIGINAL_SCRIPT}"
        )

    spec = importlib.util.spec_from_file_location(
        "original502_scorer",
        ORIGINAL_SCRIPT,
    )

    if spec is None or spec.loader is None:
        fail(
            "Could not construct import specification "
            f"for {ORIGINAL_SCRIPT}"
        )

    module = importlib.util.module_from_spec(
        spec
    )

    spec.loader.exec_module(
        module
    )

    return module


def discover_panel_specifications(module):
    """
    Locate the original script's panel specification dictionary
    without assuming its variable name.

    We require each entry to contain:
        vcf
        trait
        expected_orientation
    """

    candidates = []

    for name, value in vars(module).items():

        if not isinstance(
            value,
            dict,
        ):
            continue

        subset = {}

        for key, item in value.items():

            if not isinstance(
                item,
                dict,
            ):
                continue

            required = {
                "vcf",
                "trait",
                "expected_orientation",
            }

            if required.issubset(
                set(
                    item.keys()
                )
            ):

                subset[
                    key
                ] = item

        if subset:

            candidates.append(
                (
                    name,
                    subset,
                )
            )

    if not candidates:
        fail(
            "Could not discover the original panel "
            "specification dictionary."
        )

    # Prefer the dictionary containing the most panel specs.
    candidates.sort(
        key=lambda x: len(x[1]),
        reverse=True,
    )

    name, specifications = candidates[0]

    print(
        "\nOriginal panel specification dictionary:",
        name,
    )

    print(
        "Available panel keys:",
        ", ".join(
            str(x)
            for x in specifications.keys()
        ),
    )

    return specifications


def normalize_panel_name(value):

    return re.sub(
        r"[^a-z0-9]",
        "",
        str(value).lower(),
    )


def find_panel_key(
    specifications,
    module,
):

    target = normalize_panel_name(
        module
    )

    exact = []

    partial = []

    for key in specifications:

        normalized = normalize_panel_name(
            key
        )

        if normalized == target:
            exact.append(
                key
            )

        elif target in normalized:
            partial.append(
                key
            )

    if len(exact) == 1:
        return exact[0]

    if len(partial) == 1:
        return partial[0]

    # Trait alone is insufficient for Ray1/Ray2 or Dennis1/Dennis2,
    # so do not guess here.
    fail(
        f"Could not uniquely map module {module!r} "
        f"to panel keys {list(specifications.keys())}"
    )


def load_metadata(original):
    """
    Reconstruct the metadata dictionary expected by score_panel(),
    using the original parse_phenotype() implementation.

    This gives us the original phenotype parsing while allowing the
    known Band metadata corrections to be applied explicitly.
    """

    if not META.exists():
        fail(
            f"Missing metadata: {META}"
        )

    metadata = {}

    with META.open(
        newline="",
        errors="replace",
    ) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:

            accession = clean(
                row.get(
                    "Accession"
                )
            )

            if not accession:
                continue

            taxon = clean(
                row.get(
                    "Taxon Name"
                )
            )

            metadata[
                accession
            ] = {
                "taxon":
                    taxon,

                "data_type":
                    clean(
                        row.get(
                            "Data Type"
                        )
                    ),

                "B":
                    original.parse_phenotype(
                        row.get(
                            "B"
                        )
                    ),

                "R":
                    original.parse_phenotype(
                        row.get(
                            "R"
                        )
                    ),

                "D":
                    original.parse_phenotype(
                        row.get(
                            "D"
                        )
                    ),
            }

    # --------------------------------------------------------
    # Locked Band phenotype corrections.
    # --------------------------------------------------------

    corrected_total = 0
    changed_total = 0

    for accession, entry in metadata.items():

        taxon_lower = (
            entry[
                "taxon"
            ].lower()
        )

        if (
            "malleti" in taxon_lower
            or
            "ecuadoriensis" in taxon_lower
        ):

            corrected_total += 1

            if (
                entry["B"]
                != 0
            ):
                changed_total += 1

            entry["B"] = 0

    print(
        "\nBand phenotype correction:"
    )

    print(
        "  malleti/ecuadoriensis rows forced to B=0:",
        corrected_total,
    )

    print(
        "  values actually changed:",
        changed_total,
    )

    return metadata


def load_final_filtered_manifest():

    if not FILTER_MANIFEST.exists():
        fail(
            f"Missing final cohort manifest: {FILTER_MANIFEST}"
        )

    d = pd.read_csv(
        FILTER_MANIFEST,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "cohort",
        "module",
        "sample",
        "taxon_name",
        "data_type",
        "keep",
    }

    missing = required - set(
        d.columns
    )

    if missing:
        fail(
            "Filtered cohort manifest is missing columns: "
            + ", ".join(
                sorted(
                    missing
                )
            )
        )

    d[
        "keep_bool"
    ] = d[
        "keep"
    ].map(
        as_bool
    )

    retained = d.loc[
        d[
            "keep_bool"
        ]
    ].copy()

    # --------------------------------------------------------
    # Safety check: the final locked taxon exclusions really
    # must not survive into either cohort.
    # --------------------------------------------------------

    taxon_lower = (
        retained[
            "taxon_name"
        ]
        .astype(str)
        .str.lower()
    )

    excluded = pd.Series(
        False,
        index=retained.index,
    )

    for term in EXCLUDED_TAXON_TERMS:

        excluded = (
            excluded
            |
            taxon_lower.str.contains(
                term,
                regex=False,
                na=False,
            )
        )

    if excluded.any():

        bad = retained.loc[
            excluded,
            [
                "cohort",
                "module",
                "sample",
                "taxon_name",
            ],
        ]

        fail(
            "Final manifest still retains excluded taxa:\n"
            + bad.to_string(
                index=False
            )
        )

    return d


def make_filtered_cohort_function(
    allowed_samples,
):
    """
    Return a make_cohort() function with the same interface expected
    by the original script.

    Everything about SNP scoring remains original-script logic.
    Only which comparative individuals enter the cohort is changed.
    """

    allowed_samples = set(
        allowed_samples
    )

    def filtered_make_cohort(
        samples,
        metadata,
        trait,
    ):

        included = []
        manifest = []

        for sample in samples:

            entry = metadata.get(
                sample
            )

            if entry is None:

                phenotype = None
                reason = (
                    "no_metadata_match"
                )

            else:

                phenotype = entry.get(
                    trait
                )

                if (
                    sample
                    not in allowed_samples
                ):

                    reason = (
                        "filtered_out"
                    )

                elif phenotype is None:

                    reason = (
                        "missing_phenotype"
                    )

                else:

                    reason = (
                        "included"
                    )

                    included.append(
                        {
                            "sample":
                                sample,

                            "phenotype":
                                phenotype,

                            "data_type":
                                entry[
                                    "data_type"
                                ],

                            "taxon":
                                entry[
                                    "taxon"
                                ],
                        }
                    )

            manifest.append(
                {
                    "sample":
                        sample,

                    "trait":
                        trait,

                    "phenotype":
                        (
                            ""
                            if phenotype
                            is None
                            else phenotype
                        ),

                    "reason":
                        reason,

                    "taxon":
                        (
                            ""
                            if entry
                            is None
                            else entry[
                                "taxon"
                            ]
                        ),

                    "data_type":
                        (
                            ""
                            if entry
                            is None
                            else entry[
                                "data_type"
                            ]
                        ),
                }
            )

        return (
            included,
            manifest,
        )

    return filtered_make_cohort


def unpack_score_panel_result(
    result,
):
    """
    Identify score_panel() return objects by their contents rather
    than assuming a particular tuple order.

    Expected objects:
      - score rows: list of dicts containing scheme/chrom/pos
      - manifest rows: list of dicts containing sample/reason
      - cohort summary: dict
    """

    if not isinstance(
        result,
        tuple,
    ):
        fail(
            "score_panel() returned an unexpected "
            f"object type: {type(result)}"
        )

    rows = None
    manifest = None
    cohort_summary = None

    unidentified_lists = []

    for part in result:

        if isinstance(
            part,
            dict,
        ):
            if cohort_summary is not None:
                fail(
                    "score_panel() returned more than one "
                    "dictionary-like summary object."
                )

            cohort_summary = part
            continue

        if isinstance(
            part,
            list,
        ):

            if len(part) == 0:
                unidentified_lists.append(
                    part
                )
                continue

            first = part[0]

            if not isinstance(
                first,
                dict,
            ):
                fail(
                    "score_panel() returned a list whose "
                    "elements are not dictionaries."
                )

            keys = set(
                first.keys()
            )

            if {
                "scheme",
                "chrom",
                "pos",
            }.issubset(keys):

                if rows is not None:
                    fail(
                        "More than one score-row list detected."
                    )

                rows = part

            elif {
                "sample",
                "reason",
            }.issubset(keys):

                if manifest is not None:
                    fail(
                        "More than one manifest-row list detected."
                    )

                manifest = part

            else:
                unidentified_lists.append(
                    part
                )

            continue

        fail(
            "Unexpected object returned by score_panel(): "
            f"{type(part)}"
        )

    if rows is None:
        fail(
            "Could not identify score rows returned by "
            "score_panel()."
        )

    if manifest is None:

        if len(
            unidentified_lists
        ) == 1:
            manifest = unidentified_lists[0]

        else:
            fail(
                "Could not uniquely identify manifest rows "
                "returned by score_panel()."
            )

    if cohort_summary is None:
        fail(
            "Could not identify cohort summary returned by "
            "score_panel()."
        )

    # Strong validation before anything is written.
    if not all(
        isinstance(
            row,
            dict,
        )
        for row in rows
    ):
        fail(
            "Score-row object contains non-dictionary entries."
        )

    if not all(
        isinstance(
            row,
            dict,
        )
        for row in manifest
    ):
        fail(
            "Manifest object contains non-dictionary entries."
        )

    return (
        rows,
        manifest,
        cohort_summary,
    )


def numeric_series(
    series,
):

    return pd.to_numeric(
        series,
        errors="coerce",
    )


def prepare_raw_scores(
    df,
):

    for column in [
        "score_eligible",
        "callability_eligible",
        "matches_reference_expectation",
    ]:

        if column in df.columns:

            df[
                column + "_bool"
            ] = df[
                column
            ].map(
                as_bool
            )

    for column in [
        "balanced_accuracy",
        "MCC",
        "tp",
        "tn",
        "fp",
        "fn",
        "misclassified_n",
        "call_rate",
        "phenotype_0_call_rate",
        "phenotype_1_call_rate",
    ]:

        if column in df.columns:

            df[
                column
            ] = numeric_series(
                df[
                    column
                ]
            )

    denominator = (
        df[
            "tp"
        ].fillna(0)
        +
        df[
            "tn"
        ].fillna(0)
        +
        df[
            "fp"
        ].fillna(0)
        +
        df[
            "fn"
        ].fillna(0)
    )

    correct = (
        df[
            "tp"
        ].fillna(0)
        +
        df[
            "tn"
        ].fillna(0)
    )

    df[
        "match_fraction"
    ] = np.where(
        denominator > 0,
        correct / denominator,
        np.nan,
    )

    df[
        "expected_pattern_eligible"
    ] = (
        df[
            "score_eligible_bool"
        ]
        &
        df[
            "matches_reference_expectation_bool"
        ]
    )

    return df


def make_primary_table(
    raw,
):

    available_schemes = sorted(
        raw[
            "scheme"
        ]
        .dropna()
        .astype(str)
        .unique()
    )

    print(
        "\nQC schemes produced by original scorer:",
        ", ".join(
            available_schemes
        ),
    )

    if (
        PRIMARY_SCHEME
        in available_schemes
    ):

        scheme = PRIMARY_SCHEME

    else:

        fail(
            f"Expected primary QC scheme {PRIMARY_SCHEME!r} "
            f"not found. Available: {available_schemes}"
        )

    primary = raw.loc[
        raw[
            "scheme"
        ]
        == scheme
    ].copy()

    return primary


def make_across_scheme_table(
    raw,
):
    """
    Conservative site-level summary across all original QC schemes.

    This gives us a quantity analogous to the previous
    'minimum BA' idea:

      - minimum BA across schemes
      - whether every scheme was score eligible
      - whether every scheme selected the expected biological
        orientation

    No focal genotype fields enter this calculation.
    """

    group_cols = [
        "analysis_cohort",
        "module",
        "panel",
        "trait",
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    rows = []

    for keys, d in raw.groupby(
        group_cols,
        sort=False,
        dropna=False,
    ):

        record = dict(
            zip(
                group_cols,
                keys,
            )
        )

        ba = numeric_series(
            d[
                "balanced_accuracy"
            ]
        )

        eligible = d[
            "score_eligible_bool"
        ].astype(bool)

        expected = d[
            "matches_reference_expectation_bool"
        ].astype(bool)

        expected_eligible = (
            eligible
            &
            expected
        )

        record.update(
            {
                "scheme_n":
                    len(d),

                "schemes":
                    ",".join(
                        sorted(
                            d[
                                "scheme"
                            ]
                            .astype(str)
                            .unique()
                        )
                    ),

                "all_schemes_score_eligible":
                    bool(
                        eligible.all()
                    ),

                "all_schemes_match_expected_orientation":
                    bool(
                        expected.all()
                    ),

                "all_schemes_expected_pattern_eligible":
                    bool(
                        expected_eligible.all()
                    ),

                "minimum_BA_across_schemes":
                    (
                        ba.min()
                        if ba.notna().any()
                        else np.nan
                    ),

                "maximum_BA_across_schemes":
                    (
                        ba.max()
                        if ba.notna().any()
                        else np.nan
                    ),

                "minimum_match_fraction_across_schemes":
                    (
                        numeric_series(
                            d[
                                "match_fraction"
                            ]
                        ).min()
                    ),

                "maximum_misclassified_n_across_schemes":
                    (
                        numeric_series(
                            d[
                                "misclassified_n"
                            ]
                        ).max()
                    ),
            }
        )

        rows.append(
            record
        )

    return pd.DataFrame(
        rows
    )


def build_threshold_summary(
    primary,
    paired,
):

    output = []

    for cohort in COHORTS:

        for module in MODULES:

            p = primary.loc[
                (
                    primary[
                        "analysis_cohort"
                    ]
                    == cohort
                )
                &
                (
                    primary[
                        "module"
                    ]
                    == module
                )
            ].copy()

            q = paired.loc[
                (
                    paired[
                        "analysis_cohort"
                    ]
                    == cohort
                )
                &
                (
                    paired[
                        "module"
                    ]
                    == module
                )
            ].copy()

            p_ba = numeric_series(
                p[
                    "balanced_accuracy"
                ]
            )

            p_pattern = p[
                "expected_pattern_eligible"
            ].astype(bool)

            q_ba = numeric_series(
                q[
                    "minimum_BA_across_schemes"
                ]
            )

            q_pattern = q[
                "all_schemes_expected_pattern_eligible"
            ].astype(bool)

            row = {
                "cohort":
                    cohort,

                "module":
                    module,

                "candidate_SNPs_primary_scheme":
                    len(p),

                "score_eligible_primary":
                    int(
                        p[
                            "score_eligible_bool"
                        ].sum()
                    ),

                "expected_orientation_primary":
                    int(
                        p[
                            "matches_reference_expectation_bool"
                        ].sum()
                    ),

                "expected_pattern_eligible_primary":
                    int(
                        p_pattern.sum()
                    ),

                "BA_median_expected_pattern_primary":
                    (
                        p_ba[
                            p_pattern
                        ].median()
                    ),

                "BA_max_expected_pattern_primary":
                    (
                        p_ba[
                            p_pattern
                        ].max()
                    ),

                "candidate_SNPs_across_schemes":
                    len(q),

                "expected_pattern_eligible_all_schemes":
                    int(
                        q_pattern.sum()
                    ),

                "minimum_BA_median_all_schemes":
                    (
                        q_ba[
                            q_pattern
                        ].median()
                    ),

                "minimum_BA_max_all_schemes":
                    (
                        q_ba[
                            q_pattern
                        ].max()
                    ),
            }

            for threshold in [
                0.80,
                0.85,
                0.90,
                0.95,
            ]:

                suffix = (
                    f"{threshold:.2f}"
                )

                row[
                    "primary_expected_pattern_BA_ge_"
                    + suffix
                ] = int(
                    (
                        p_pattern
                        &
                        (
                            p_ba
                            >= threshold
                        )
                    ).sum()
                )

                row[
                    "all_schemes_expected_pattern_minBA_ge_"
                    + suffix
                ] = int(
                    (
                        q_pattern
                        &
                        (
                            q_ba
                            >= threshold
                        )
                    ).sum()
                )

            output.append(
                row
            )

    return pd.DataFrame(
        output
    )


def main():

    original = load_original_module()

    specifications = (
        discover_panel_specifications(
            original
        )
    )

    panel_from_module = {}

    for module in MODULES:

        panel_from_module[
            module
        ] = find_panel_key(
            specifications,
            module,
        )

    print(
        "\nModule -> original panel mapping:"
    )

    for module in MODULES:

        print(
            f"  {module:8s} -> "
            f"{panel_from_module[module]}"
        )

    metadata = load_metadata(
        original
    )

    filtered_manifest = (
        load_final_filtered_manifest()
    )

    retained = filtered_manifest.loc[
        filtered_manifest[
            "keep_bool"
        ]
    ].copy()

    print(
        "\nFinal retained cohort sizes:"
    )

    print(
        retained.groupby(
            [
                "cohort",
                "module",
            ]
        )
        .size()
        .to_string()
    )

    raw_rows = []
    score_manifest_rows = []
    cohort_rows = []

    original_make_cohort = (
        original.make_cohort
    )

    try:

        for cohort_name in COHORTS:

            for module in MODULES:

                panel = panel_from_module[
                    module
                ]

                specification = (
                    specifications[
                        panel
                    ]
                )

                allowed = set(
                    retained.loc[
                        (
                            retained[
                                "cohort"
                            ]
                            == cohort_name
                        )
                        &
                        (
                            retained[
                                "module"
                            ]
                            == module
                        ),
                        "sample",
                    ]
                )

                print(
                    "\n=============================================="
                )

                print(
                    cohort_name,
                    "/",
                    module,
                    "/",
                    panel,
                )

                print(
                    "Allowed comparative samples:",
                    len(
                        allowed
                    ),
                )

                # Replace only cohort construction.
                original.make_cohort = (
                    make_filtered_cohort_function(
                        allowed
                    )
                )

                result = original.score_panel(
                    panel,
                    specification,
                    metadata,
                )

                (
                    rows,
                    manifest,
                    cohort_summary,
                ) = unpack_score_panel_result(
                    result
                )

                print(
                    "Scored rows:",
                    len(
                        rows
                    ),
                )

                for row in rows:

                    raw_rows.append(
                        {
                            "analysis_cohort":
                                cohort_name,

                            "module":
                                module,

                            **row,
                        }
                    )

                for row in manifest:

                    score_manifest_rows.append(
                        {
                            "analysis_cohort":
                                cohort_name,

                            "module":
                                module,

                            **row,
                        }
                    )

                # Our own cohort summary uses the locked manifest,
                # rather than relying on reason counts internal to 41.
                subset = retained.loc[
                    (
                        retained[
                            "cohort"
                        ]
                        == cohort_name
                    )
                    &
                    (
                        retained[
                            "module"
                        ]
                        == module
                    )
                ].copy()

                trait = (
                    TRAIT_FROM_MODULE[
                        module
                    ]
                )

                phenotypes = []

                for sample in subset[
                    "sample"
                ]:

                    entry = metadata.get(
                        sample
                    )

                    phenotype = (
                        None
                        if entry
                        is None
                        else entry.get(
                            trait
                        )
                    )

                    phenotypes.append(
                        phenotype
                    )

                cohort_rows.append(
                    {
                        "cohort":
                            cohort_name,

                        "module":
                            module,

                        "original_panel":
                            panel,

                        "trait":
                            trait,

                        "retained_n":
                            len(
                                subset
                            ),

                        "WGS_n":
                            int(
                                (
                                    subset[
                                        "data_type"
                                    ]
                                    == "WGS"
                                ).sum()
                            ),

                        "CS_n":
                            int(
                                (
                                    subset[
                                        "data_type"
                                    ]
                                    == "CS"
                                ).sum()
                            ),

                        "phenotype_0_n":
                            sum(
                                x == 0
                                for x
                                in phenotypes
                            ),

                        "phenotype_1_n":
                            sum(
                                x == 1
                                for x
                                in phenotypes
                            ),
                    }
                )

    finally:

        # Do not leave imported module modified.
        original.make_cohort = (
            original_make_cohort
        )

    raw = pd.DataFrame(
        raw_rows
    )

    if raw.empty:
        fail(
            "No all-SNP score rows produced."
        )

    raw = prepare_raw_scores(
        raw
    )

    raw.to_csv(
        OUT_RAW,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    primary = make_primary_table(
        raw
    )

    primary.to_csv(
        OUT_PRIMARY,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    paired = make_across_scheme_table(
        raw
    )

    paired.to_csv(
        OUT_PAIRED,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary = build_threshold_summary(
        primary,
        paired,
    )

    summary.to_csv(
        OUT_SUMMARY,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    pd.DataFrame(
        cohort_rows
    ).to_csv(
        OUT_COHORT_SUMMARY,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    pd.DataFrame(
        score_manifest_rows
    ).to_csv(
        OUT_SCORE_MANIFEST,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print(
        "\n=============================================="
    )

    print(
        "ALL-SNP RESCORING SUMMARY"
    )

    print(
        "=============================================="
    )

    display_cols = [
        "cohort",
        "module",
        "candidate_SNPs_primary_scheme",
        "score_eligible_primary",
        "expected_pattern_eligible_primary",
        "primary_expected_pattern_BA_ge_0.80",
        "primary_expected_pattern_BA_ge_0.85",
        "primary_expected_pattern_BA_ge_0.90",
        "primary_expected_pattern_BA_ge_0.95",
        "BA_max_expected_pattern_primary",
        "all_schemes_expected_pattern_minBA_ge_0.80",
        "all_schemes_expected_pattern_minBA_ge_0.90",
        "minimum_BA_max_all_schemes",
    ]

    print(
        summary[
            display_cols
        ].to_string(
            index=False
        )
    )

    print(
        "\nOutputs:"
    )

    for path in [
        OUT_RAW,
        OUT_PRIMARY,
        OUT_PAIRED,
        OUT_SUMMARY,
        OUT_COHORT_SUMMARY,
        OUT_SCORE_MANIFEST,
    ]:

        print(
            f"  {path}"
        )


if __name__ == "__main__":
    main()
