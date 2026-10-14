#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

SOURCE_ANALYSIS = (
    "B_heurippa_besckei_no_DR_positive_reanalysis"
)

NEW_ANALYSIS = (
    "B_west_of_Andes_no_DR_positive_reanalysis"
)

SOURCE_CHECKPOINT = (
    PROJECT
    / "analysis"
    / SOURCE_ANALYSIS
    / "checkpoints"
)

OUTDIR = (
    PROJECT
    / "analysis"
    / NEW_ANALYSIS
    / "checkpoints"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)

SOURCE_MANIFEST = (
    SOURCE_CHECKPOINT
    / (
        "original502_filtered_cohort_manifest."
        "B_heurippa_besckei_no_DR_positive."
        "WGS_vs_WGSplusCS.tsv"
    )
)

NEW_MANIFEST = (
    OUTDIR
    / (
        "original502_filtered_cohort_manifest."
        "B_west_of_Andes_no_DR_positive."
        "WGS_vs_WGSplusCS.tsv"
    )
)


SOURCE_BASE = (
    PROJECT
    / "scripts"
    / "110b_base_scorer_B_no_DR_positive.py"
)

SOURCE_SCAN = (
    PROJECT
    / "scripts"
    / "110c_scan_B_no_DR_positive.py"
)

SOURCE_PLOT = (
    PROJECT
    / "scripts"
    / "110d_plot_B_no_DR_positive_BA_and_focal.R"
)


BASE_OUTPUT = (
    PROJECT
    / "scripts"
    / "119b_base_scorer_B_west_of_Andes.py"
)

SCAN_OUTPUT = (
    PROJECT
    / "scripts"
    / "119c_scan_B_west_of_Andes.py"
)

PLOT_OUTPUT = (
    PROJECT
    / "scripts"
    / "119d_plot_B_west_of_Andes_BA_and_focal.R"
)


# ------------------------------------------------------------
# These are the additional exclusions agreed for the
# operational west-of-Andes Band sensitivity analysis.
#
# Everything removed by the previous no-D/no-R analysis stays
# removed automatically because we start from its manifest.
# ------------------------------------------------------------

WEST_EXCLUSION_TERMS = (
    "amaryllis",
    "burchelli",
    "nanna",
    "xenoclea",
    "plesseni",
    "besckei",
)


def as_bool(value):

    return str(value).strip().lower() in {
        "true",
        "t",
        "1",
        "yes",
        "y",
    }


def fail(message):

    raise RuntimeError(
        message
    )


def main():

    # ========================================================
    # Validate inputs
    # ========================================================

    for path in (
        SOURCE_MANIFEST,
        SOURCE_BASE,
        SOURCE_SCAN,
        SOURCE_PLOT,
    ):

        if not path.exists():

            fail(
                f"Missing required input: {path}"
            )


    # ========================================================
    # Start from yesterday's already-filtered cohort.
    # ========================================================

    d = pd.read_csv(
        SOURCE_MANIFEST,
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
            "Source manifest missing columns: "
            + ", ".join(
                sorted(
                    missing
                )
            )
        )


    d[
        "prior_keep"
    ] = d[
        "keep"
    ].map(
        as_bool
    )


    taxon_lower = (
        d[
            "taxon_name"
        ]
        .astype(str)
        .str.lower()
    )


    # ========================================================
    # Explicit west-of-Andes operational exclusions.
    # ========================================================

    d[
        "west_named_exclusion"
    ] = False

    d[
        "west_named_exclusion_term"
    ] = ""


    for term in WEST_EXCLUSION_TERMS:

        hit = taxon_lower.str.contains(
            term,
            regex=False,
            na=False,
        )

        d.loc[
            hit,
            "west_named_exclusion",
        ] = True

        empty = (
            hit
            &
            (
                d[
                    "west_named_exclusion_term"
                ]
                == ""
            )
        )

        d.loc[
            empty,
            "west_named_exclusion_term",
        ] = term


    # New cohort = exact previous cohort minus named taxa.
    d[
        "keep"
    ] = (
        d[
            "prior_keep"
        ]
        &
        ~d[
            "west_named_exclusion"
        ]
    )


    # ========================================================
    # Safety checks
    # ========================================================

    retained = d.loc[
        d[
            "keep"
        ]
    ].copy()


    retained_taxon_lower = (
        retained[
            "taxon_name"
        ]
        .astype(str)
        .str.lower()
    )


    bad = pd.Series(
        False,
        index=retained.index,
    )

    for term in WEST_EXCLUSION_TERMS:

        bad = (
            bad
            |
            retained_taxon_lower.str.contains(
                term,
                regex=False,
                na=False,
            )
        )


    if bad.any():

        fail(
            "Named west-of-Andes exclusions survived:\n"
            + retained.loc[
                bad,
                [
                    "cohort",
                    "module",
                    "sample",
                    "taxon_name",
                ],
            ].to_string(
                index=False
            )
        )


    # Source-manifest filtering must only become stricter.
    resurrected = (
        d[
            "keep"
        ]
        &
        ~d[
            "prior_keep"
        ]
    )

    if resurrected.any():

        fail(
            "A previously excluded sample was accidentally "
            "reintroduced."
        )


    # ========================================================
    # Write manifest
    # ========================================================

    d.to_csv(
        NEW_MANIFEST,
        sep="\t",
        index=False,
    )


    # ========================================================
    # Print exact cohort changes.
    # ========================================================

    print(
        "\n=============================================="
    )
    print(
        "PRIOR NO-D/NO-R COHORT"
    )
    print(
        "=============================================="
    )

    print(
        d.loc[
            d[
                "prior_keep"
            ]
        ]
        .groupby(
            [
                "cohort",
                "module",
            ]
        )
        .size()
        .to_string()
    )


    print(
        "\n=============================================="
    )
    print(
        "WEST-OF-ANDES COHORT"
    )
    print(
        "=============================================="
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


    newly_removed = d.loc[
        d[
            "prior_keep"
        ]
        &
        ~d[
            "keep"
        ]
    ].copy()


    print(
        "\n=============================================="
    )
    print(
        "NEWLY REMOVED ROWS"
    )
    print(
        "=============================================="
    )

    print(
        "Manifest rows:",
        len(
            newly_removed
        ),
    )

    print(
        "Unique samples:",
        newly_removed[
            "sample"
        ].nunique(),
    )


    if len(
        newly_removed
    ):

        print(
            "\nBy exclusion term / data type:"
        )

        print(
            newly_removed.groupby(
                [
                    "west_named_exclusion_term",
                    "data_type",
                ]
            )
            .size()
            .to_string()
        )


        print(
            "\nUnique samples by taxon:"
        )

        sample_taxa = (
            newly_removed[
                [
                    "sample",
                    "taxon_name",
                    "data_type",
                ]
            ]
            .drop_duplicates()
        )

        print(
            sample_taxa.groupby(
                [
                    "taxon_name",
                    "data_type",
                ]
            )
            .size()
            .to_string()
        )


    # ========================================================
    # Clone the exact established scorer.
    # ========================================================

    base_text = SOURCE_BASE.read_text()

    base_text = base_text.replace(
        SOURCE_ANALYSIS,
        NEW_ANALYSIS,
    )

    base_text = base_text.replace(
        (
            "original502_filtered_cohort_manifest."
            "B_heurippa_besckei_no_DR_positive."
            "WGS_vs_WGSplusCS.tsv"
        ),
        NEW_MANIFEST.name,
    )


    # Extend its internal safety assertion with the six newly
    # excluded taxon-name terms.
    old_safety = '''EXCLUDED_TAXON_TERMS = (
    "numata",
    "pachinus",
    "hecale",
    "ismenius",
    "timareta",
)'''

    new_safety = '''EXCLUDED_TAXON_TERMS = (
    "numata",
    "pachinus",
    "hecale",
    "ismenius",
    "timareta",
    "amaryllis",
    "burchelli",
    "nanna",
    "xenoclea",
    "plesseni",
    "besckei",
)'''

    if old_safety not in base_text:

        fail(
            "Could not locate EXCLUDED_TAXON_TERMS "
            "in 110b; no derivative scripts written."
        )

    base_text = base_text.replace(
        old_safety,
        new_safety,
        1,
    )

    BASE_OUTPUT.write_text(
        base_text
    )


    # ========================================================
    # Clone the exact established broad Band scanner.
    # ========================================================

    scan_text = SOURCE_SCAN.read_text()

    scan_text = scan_text.replace(
        SOURCE_ANALYSIS,
        NEW_ANALYSIS,
    )

    scan_text = scan_text.replace(
        "110b_base_scorer_B_no_DR_positive.py",
        BASE_OUTPUT.name,
    )

    scan_text = scan_text.replace(
        "Final Band cohort sizes from 94b:",
        "Final west-of-Andes Band cohort sizes:",
    )

    SCAN_OUTPUT.write_text(
        scan_text
    )


    # ========================================================
    # Clone the exact established figure.
    #
    # No aesthetic/layout change: only routing and title.
    # ========================================================

    plot_text = SOURCE_PLOT.read_text()

    plot_text = plot_text.replace(
        SOURCE_ANALYSIS,
        NEW_ANALYSIS,
    )

    plot_text = plot_text.replace(
        "broad_B_heurippa_besckei_no_DR_positive_BA_manhattan",
        "broad_B_west_of_Andes_no_DR_positive_BA_manhattan",
    )

    plot_text = plot_text.replace(
        (
            '"Band association after excluding D-positive '
            'and R-positive samples",'
        ),
        (
            '"Band association in the west-of-Andes subset '
            'after excluding D-positive and R-positive samples",'
        ),
    )

    PLOT_OUTPUT.write_text(
        plot_text
    )


    print(
        "\n=============================================="
    )
    print(
        "OUTPUTS"
    )
    print(
        "=============================================="
    )

    print(
        "Manifest:",
        NEW_MANIFEST,
    )

    print(
        "Base scorer:",
        BASE_OUTPUT,
    )

    print(
        "Broad scanner:",
        SCAN_OUTPUT,
    )

    print(
        "Plotter:",
        PLOT_OUTPUT,
    )


if __name__ == "__main__":
    main()
