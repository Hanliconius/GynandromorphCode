#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

CHECKPOINTS = (
    PROJECT
    / "analysis/original_502/checkpoints"
)

OPTIX_SCORES = (
    CHECKPOINTS
    / "optix_original_502_comparative_scores.tsv"
)

YB_SCORES = (
    CHECKPOINTS
    / "ivory_original_502_comparative_scores.tsv"
)

OPTIX_OUTPUT = (
    CHECKPOINTS
    / "optix_comparative_informative_sites_BA80.tsv"
)

YB_OUTPUT = (
    CHECKPOINTS
    / "ivory_comparative_informative_sites_BA80.tsv"
)

SUMMARY_OUTPUT = (
    CHECKPOINTS
    / "comparative_informative_site_selection_BA80_summary.tsv"
)

PLOT_SCRIPTS = [
    PROJECT
    / "scripts/43_make_aligned_optix_genoplot_score_panels.py",
    PROJECT
    / "scripts/45_make_aligned_yb_genoplot_score_panels.py",
]

MINIMUM_BALANCED_ACCURACY = 0.80
MAXIMUM_Q = 0.05

OPTIX_MODELS = {
    "B": {
        "inheritance_model": "dominant",
        "comparative_presence_rule": (
            "heterozygous_or_homozygous_"
            "carrier_of_associated_allele"
        ),
        "focal_expectation_not_used_for_selection": (
            "left_no_band=0_associated_copies;"
            "right_band=1_or_2_associated_copies;"
            "right_heterozygosity_is_plausible_"
            "under_B_vs_DR_haplotype_model"
        ),
    },
    "D": {
        "inheritance_model": "dominant",
        "comparative_presence_rule": (
            "heterozygous_or_homozygous_"
            "carrier_of_associated_allele"
        ),
        "focal_expectation_not_used_for_selection": (
            "left_D_present=1_or_2_associated_copies;"
            "right_D_present=1_or_2_associated_copies;"
            "left_homozygous_and_right_heterozygous_"
            "is_plausible_under_DR_vs_BDR_model"
        ),
    },
    "R": {
        "inheritance_model": "dominant",
        "comparative_presence_rule": (
            "heterozygous_or_homozygous_"
            "carrier_of_associated_allele"
        ),
        "focal_expectation_not_used_for_selection": (
            "left_R_present=1_or_2_associated_copies;"
            "right_R_present=1_or_2_associated_copies;"
            "left_homozygous_and_right_heterozygous_"
            "is_plausible_under_DR_vs_BDR_model"
        ),
    },
}

YB_MODELS = {
    "Yb_dorsal": {
        "inheritance_model": "recessive",
        "comparative_presence_rule": (
            "homozygous_for_associated_allele_only"
        ),
        "focal_expectation_not_used_for_selection": (
            "left_shadow_interpreted_as_heterozygous"
            "=1_associated_copy;"
            "right_Yb_present=2_associated_copies"
        ),
    },
    "Yb_ventral": {
        "inheritance_model": "recessive",
        "comparative_presence_rule": (
            "homozygous_for_associated_allele_only"
        ),
        "focal_expectation_not_used_for_selection": (
            "left_shadow_interpreted_as_heterozygous"
            "=1_associated_copy;"
            "right_Yb_present=2_associated_copies"
        ),
    },
}


def as_float(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None

    if math.isfinite(result):
        return result

    return None


def as_bool(value):
    return (
        str(value)
        .strip()
        .lower()
        in {
            "true",
            "1",
            "yes",
        }
    )


def read_tsv(path):
    with path.open(
        newline=""
    ) as handle:
        return list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )


def write_tsv(
    path,
    rows,
    fieldnames,
):
    with path.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)


def select_sites(
    score_path,
    dataset,
):
    rows = read_tsv(
        score_path
    )

    key_field = (
        "panel"
        if dataset == "optix"
        else "trait"
    )

    models = (
        OPTIX_MODELS
        if dataset == "optix"
        else YB_MODELS
    )

    paired = defaultdict(dict)

    for row in rows:
        key = (
            row[key_field],
            row["trait"],
            row["chrom"],
            int(row["pos"]),
            row["ref"],
            row["alt"],
        )

        paired[key][
            row["scheme"]
        ] = row

    selected = []

    for key, schemes in paired.items():
        raw = schemes.get(
            "called_GT"
        )

        qc = schemes.get(
            "DP2_GQ20"
        )

        if raw is None or qc is None:
            continue

        if not (
            as_bool(
                raw.get(
                    "score_eligible"
                )
            )
            and as_bool(
                qc.get(
                    "score_eligible"
                )
            )
        ):
            continue

        if (
            raw.get(
                "selected_orientation"
            )
            != qc.get(
                "selected_orientation"
            )
        ):
            continue

        raw_q = as_float(
            raw.get("BH_q")
        )

        qc_q = as_float(
            qc.get("BH_q")
        )

        if (
            raw_q is None
            or qc_q is None
            or raw_q > MAXIMUM_Q
            or qc_q > MAXIMUM_Q
        ):
            continue

        raw_ba = as_float(
            raw.get(
                "balanced_accuracy"
            )
        )

        qc_ba = as_float(
            qc.get(
                "balanced_accuracy"
            )
        )

        if (
            raw_ba is None
            or qc_ba is None
        ):
            continue

        minimum_ba = min(
            raw_ba,
            qc_ba,
        )

        if (
            minimum_ba
            < MINIMUM_BALANCED_ACCURACY
        ):
            continue

        (
            group,
            trait,
            chrom,
            position,
            ref,
            alt,
        ) = key

        model = models[trait]

        raw_missingness_q = as_float(
            raw.get(
                "missingness_BH_q"
            )
        )

        qc_missingness_q = as_float(
            qc.get(
                "missingness_BH_q"
            )
        )

        missingness_ok = (
            raw_missingness_q is not None
            and qc_missingness_q is not None
            and raw_missingness_q > 0.05
            and qc_missingness_q > 0.05
        )

        selected.append(
            {
                key_field: group,
                "trait": trait,
                "chrom": chrom,
                "pos": position,
                "ref": ref,
                "alt": alt,
                "selected_orientation": (
                    qc[
                        "selected_orientation"
                    ]
                ),
                "inheritance_model": (
                    model[
                        "inheritance_model"
                    ]
                ),
                "comparative_presence_rule": (
                    model[
                        "comparative_presence_rule"
                    ]
                ),
                "raw_BA": raw_ba,
                "QC_BA": qc_ba,
                "minimum_BA": minimum_ba,
                "raw_MCC": raw.get(
                    "MCC",
                    "",
                ),
                "QC_MCC": qc.get(
                    "MCC",
                    "",
                ),
                "raw_q": raw_q,
                "QC_q": qc_q,
                "raw_missingness_q": (
                    raw_missingness_q
                ),
                "QC_missingness_q": (
                    qc_missingness_q
                ),
                "missingness_ok": (
                    missingness_ok
                ),
                "focal_expectation_not_used_for_selection": (
                    model[
                        "focal_expectation_not_used_for_selection"
                    ]
                ),
            }
        )

    selected.sort(
        key=lambda row: (
            row[key_field],
            -row["minimum_BA"],
            row["pos"],
        )
    )

    return (
        selected,
        key_field,
    )


def patch_plot_scripts():
    replacements = {
        (
            "optix_original_502_"
            "joint_focal_candidates.tsv"
        ): (
            "optix_comparative_"
            "informative_sites_BA80.tsv"
        ),
        (
            "ivory_original_502_"
            "inclusive_focal_candidates.tsv"
        ): (
            "ivory_comparative_"
            "informative_sites_BA80.tsv"
        ),
        (
            "Comparative-supported "
            "+ focal-matching"
        ): (
            "Comparative-informative "
            "(min BA >= 0.80)"
        ),
        (
            "Inclusive focal-compatible site"
        ): (
            "Comparative-informative "
            "(min BA >= 0.80)"
        ),
    }

    for script in PLOT_SCRIPTS:
        if not script.is_file():
            print(
                "WARNING: plotting script "
                f"missing: {script}"
            )

            continue

        original = script.read_text()

        backup = script.with_suffix(
            script.suffix
            + ".before_independent_site_selection"
        )

        if not backup.exists():
            backup.write_text(
                original
            )

        updated = original

        for old, new in replacements.items():
            updated = updated.replace(
                old,
                new,
            )

        if updated != original:
            script.write_text(
                updated
            )

            print(
                "Patched plotting input/legend: "
                f"{script}"
            )

        else:
            print(
                "Plotting script already patched: "
                f"{script}"
            )


def main():
    CHECKPOINTS.mkdir(
        parents=True,
        exist_ok=True,
    )

    required = [
        OPTIX_SCORES,
        YB_SCORES,
    ]

    missing = [
        str(path)
        for path in required
        if not path.is_file()
    ]

    if missing:
        print(
            "ERROR: required score tables "
            "are missing:"
        )

        for path in missing:
            print(path)

        return

    (
        optix_rows,
        optix_key,
    ) = select_sites(
        OPTIX_SCORES,
        "optix",
    )

    (
        yb_rows,
        yb_key,
    ) = select_sites(
        YB_SCORES,
        "yb",
    )

    optix_fields = [
        "panel",
        "trait",
        "chrom",
        "pos",
        "ref",
        "alt",
        "selected_orientation",
        "inheritance_model",
        "comparative_presence_rule",
        "raw_BA",
        "QC_BA",
        "minimum_BA",
        "raw_MCC",
        "QC_MCC",
        "raw_q",
        "QC_q",
        "raw_missingness_q",
        "QC_missingness_q",
        "missingness_ok",
        "focal_expectation_not_used_for_selection",
    ]

    yb_fields = [
        "trait",
        "chrom",
        "pos",
        "ref",
        "alt",
        "selected_orientation",
        "inheritance_model",
        "comparative_presence_rule",
        "raw_BA",
        "QC_BA",
        "minimum_BA",
        "raw_MCC",
        "QC_MCC",
        "raw_q",
        "QC_q",
        "raw_missingness_q",
        "QC_missingness_q",
        "missingness_ok",
        "focal_expectation_not_used_for_selection",
    ]

    write_tsv(
        OPTIX_OUTPUT,
        optix_rows,
        optix_fields,
    )

    write_tsv(
        YB_OUTPUT,
        yb_rows,
        yb_fields,
    )

    summary = []

    expected_groups = [
        (
            "optix",
            optix_rows,
            optix_key,
            [
                "band",
                "dennis1",
                "dennis2",
                "ray1",
                "ray2",
            ],
            "dominant",
        ),
        (
            "ivory",
            yb_rows,
            yb_key,
            [
                "Yb_dorsal",
                "Yb_ventral",
            ],
            "recessive",
        ),
    ]

    for (
        dataset,
        rows,
        key_field,
        groups,
        model,
    ) in expected_groups:
        for group in groups:
            subset = [
                row
                for row in rows
                if row[key_field] == group
            ]

            minimum_ba_values = [
                row["minimum_BA"]
                for row in subset
            ]

            summary.append(
                {
                    "dataset": dataset,
                    "panel_or_trait": group,
                    "inheritance_model": model,
                    "informative_sites": (
                        len(subset)
                    ),
                    "minimum_BA_min": (
                        min(
                            minimum_ba_values
                        )
                        if minimum_ba_values
                        else None
                    ),
                    "minimum_BA_max": (
                        max(
                            minimum_ba_values
                        )
                        if minimum_ba_values
                        else None
                    ),
                    "missingness_flagged": (
                        sum(
                            not row[
                                "missingness_ok"
                            ]
                            for row in subset
                        )
                    ),
                }
            )

    write_tsv(
        SUMMARY_OUTPUT,
        summary,
        [
            "dataset",
            "panel_or_trait",
            "inheritance_model",
            "informative_sites",
            "minimum_BA_min",
            "minimum_BA_max",
            "missingness_flagged",
        ],
    )

    patch_plot_scripts()

    print()
    print(
        "COMPARATIVE-ONLY INFORMATIVE "
        "SITE SELECTION"
    )

    print(
        "No focal genotype, allelic depth, "
        "phenotype, or expected focal match "
        "was read."
    )

    print()
    print("SELECTION RULE")

    print(
        "Both called-GT and DP2/GQ20 "
        "score eligible"
    )

    print(
        "Stable selected allele orientation"
    )

    print(
        "BH q <= 0.05 in both schemes"
    )

    print(
        "Balanced accuracy >= 0.80 "
        "in both schemes"
    )

    print(
        "Missingness recorded but not "
        "used as an exclusion"
    )

    print()
    print(
        "ZYGOSITY MODELS USED IN "
        "COMPARATIVE SCORING"
    )

    print(
        "B, D, R: dominant; heterozygote "
        "or associated homozygote predicts "
        "presence"
    )

    print(
        "Yb: recessive; only associated "
        "homozygote predicts presence"
    )

    print()

    print(
        "dataset\tpanel_or_trait\tmodel\t"
        "informative_sites\tmin_BA_range\t"
        "missingness_flagged"
    )

    for row in summary:
        if (
            row["minimum_BA_min"]
            is None
        ):
            ba_range = "NA"

        else:
            ba_range = (
                f'{row["minimum_BA_min"]:.4f}-'
                f'{row["minimum_BA_max"]:.4f}'
            )

        print(
            f'{row["dataset"]}\t'
            f'{row["panel_or_trait"]}\t'
            f'{row["inheritance_model"]}\t'
            f'{row["informative_sites"]}\t'
            f'{ba_range}\t'
            f'{row["missingness_flagged"]}'
        )

    print()
    print(
        "Optix informative sites: "
        f"{OPTIX_OUTPUT}"
    )

    print(
        "Yb informative sites: "
        f"{YB_OUTPUT}"
    )

    print(
        f"Summary: {SUMMARY_OUTPUT}"
    )


if __name__ == "__main__":
    main()
