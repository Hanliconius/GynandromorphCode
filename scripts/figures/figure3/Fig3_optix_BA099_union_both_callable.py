#!/usr/bin/env python3
"""Run the Figure 3 BA>0.99 genoplots with the actual 170A/170B union logic.

A SNP is admitted when ANY of Band, Dennis, or Ray has comparative BA >0.99.
Then the existing Figure 3 code requires a posterior>=0.95 reciprocal pair call
for both TN sides. The recalled genotype state itself is never a selection rule.

For display orientation only (not selection), retain the established optix
biological coding used by the regional genoplot tables:
  Band:   REF is pattern-associated
  Ray:    ALT is pattern-associated
  Dennis: ALT is pattern-associated
"""

from pathlib import Path
import importlib.util
import sys

import pandas as pd


HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / "Fig3_optix_BA099_both_callable.py"

spec = importlib.util.spec_from_file_location(
    "fig3_ba99_base",
    BASE_PATH,
)

base = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = base
spec.loader.exec_module(base)


def load_ba_union_for_all_panels():
    keys = [
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    band = pd.read_csv(
        base.BAND_BA_PATH,
        sep="\t",
        low_memory=False,
    )

    base.require_columns(
        band,
        keys
        + [
            "analysis_cohort",
            "scheme",
            "balanced_accuracy",
            "score_eligible",
        ],
        "Analysis-142 Band BA table",
    )

    band["pos"] = pd.to_numeric(
        band["pos"],
        errors="raise",
    ).astype(int)

    band["balanced_accuracy"] = pd.to_numeric(
        band["balanced_accuracy"],
        errors="coerce",
    )

    band = band.loc[
        band["analysis_cohort"].eq(
            "WGS_plus_CS"
        )
        & band["scheme"].eq(
            "DP2_GQ20"
        )
        & band["score_eligible"].map(
            base.as_bool
        )
        & band["balanced_accuracy"].gt(
            base.MIN_BA
        )
    ].copy()

    band["support_module"] = "Band"

    dr = pd.read_csv(
        base.DR_BA_PATH,
        sep="\t",
        low_memory=False,
    )

    base.require_columns(
        dr,
        keys
        + [
            "module",
            "balanced_accuracy",
        ],
        "Wallbank Dennis/Ray BA table",
    )

    dr["pos"] = pd.to_numeric(
        dr["pos"],
        errors="raise",
    ).astype(int)

    dr["balanced_accuracy"] = pd.to_numeric(
        dr["balanced_accuracy"],
        errors="coerce",
    )

    dr = dr.loc[
        dr["module"].isin(
            [
                "Dennis",
                "Ray",
            ]
        )
        & dr["balanced_accuracy"].gt(
            base.MIN_BA
        )
    ].copy()

    dr["support_module"] = dr["module"]

    union = pd.concat(
        [
            band[
                keys
                + [
                    "support_module",
                    "balanced_accuracy",
                ]
            ],
            dr[
                keys
                + [
                    "support_module",
                    "balanced_accuracy",
                ]
            ],
        ],
        ignore_index=True,
        sort=False,
    )

    if union.empty:
        raise RuntimeError(
            "No BA >0.99 B/D/R sites were recovered."
        )

    grouped = []

    for key, group in union.groupby(
        keys,
        sort=False,
    ):
        grouped.append(
            {
                **dict(zip(keys, key)),
                "maximum_BA": float(
                    group[
                        "balanced_accuracy"
                    ].max()
                ),
                "support_modules": ",".join(
                    sorted(
                        set(
                            group[
                                "support_module"
                            ].astype(str)
                        )
                    )
                ),
            }
        )

    union = pd.DataFrame(
        grouped
    )

    print(
        "BA_GT_099_UNION_BEFORE_TN_CALLABILITY"
    )
    print(
        f"unique_sites\t{len(union)}"
    )

    print(
        "support_combination_counts"
    )
    print(
        union[
            "support_modules"
        ]
        .value_counts()
        .sort_index()
        .to_string()
    )
    print()

    # The base plotter expects a module-specific BA table.
    # Replicate the UNION into each display trait; this changes
    # only which sites are eligible to enter a panel.
    frames = []

    for module, pattern_alt in (
        ("Band", False),
        ("Ray", True),
        ("Dennis", True),
    ):
        frame = union.copy()
        frame["module"] = module
        frame["balanced_accuracy"] = (
            frame["maximum_BA"]
        )
        frame["selected_orientation"] = (
            "ALT_pattern_display"
            if pattern_alt
            else "REF_pattern_display"
        )
        frame["pattern_is_alt"] = pattern_alt
        frame["BA_source"] = (
            "170A/170B union support="
            + frame["support_modules"]
        )

        frames.append(
            frame[
                [
                    "module",
                    *keys,
                    "balanced_accuracy",
                    "selected_orientation",
                    "pattern_is_alt",
                    "BA_source",
                ]
            ]
        )

    return pd.concat(
        frames,
        ignore_index=True,
        sort=False,
    )


base.load_ba_sites = (
    load_ba_union_for_all_panels
)

base.COMBINED_PDF = (
    base.OUTPUT_DIR
    / "optix_BAgt099_union_both_recalled_callable_five_loci.pdf"
)

base.SITE_OUTPUT = (
    base.OUTPUT_DIR
    / "BAgt099_union_both_recalled_callable_sites.tsv"
)

base.MATRIX_OUTPUT = (
    base.OUTPUT_DIR
    / "BAgt099_union_both_recalled_callable_plot_values.tsv.gz"
)

base.AUDIT_OUTPUT = (
    base.OUTPUT_DIR
    / "BAgt099_union_both_recalled_callable.audit.json"
)


if __name__ == "__main__":
    base.main()
