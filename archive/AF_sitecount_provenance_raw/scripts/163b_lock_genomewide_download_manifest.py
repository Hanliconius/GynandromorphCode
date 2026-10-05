#!/usr/bin/env python3

import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

OUT = (
    PROJECT
    / "analysis"
    / "genomewide_mixing_input"
    / "metadata"
)

METRICS = (
    PROJECT
    / "analysis"
    / "autosomal_three_vs_four_founder"
    / "eligible_control_callability_and_heterozygosity.tsv"
)

REFERENCE = "GCA_900068215.1"
TARGET_CONTROLS = 18
MIN_CALL_RATE = 0.65

FOCAL_LIBRARIES = [
    {
        "accession": "SRX33804284",
        "library_id": "TN_Left_T2",
        "side": "Left",
        "leg": "T2",
        "merge_group": "TN_Left_merged",
    },
    {
        "accession": "SRX33804285",
        "library_id": "TN_Right_T2",
        "side": "Right",
        "leg": "T2",
        "merge_group": "TN_Right_merged",
    },
    {
        "accession": "SRX33804286",
        "library_id": "TN_Left_T3",
        "side": "Left",
        "leg": "T3",
        "merge_group": "TN_Left_merged",
    },
    {
        "accession": "SRX33804287",
        "library_id": "TN_Right_T3",
        "side": "Right",
        "leg": "T3",
        "merge_group": "TN_Right_merged",
    },
]


def evenly_spaced(frame, number):
    ordered = frame.sort_values(
        [
            "heterozygosity_among_callable",
            "sample",
        ],
        kind="stable",
    ).reset_index(drop=True)

    if len(ordered) <= number:
        return ordered.copy()

    indices = np.rint(
        np.linspace(
            0,
            len(ordered) - 1,
            number,
        )
    ).astype(int)

    if len(np.unique(indices)) != number:
        raise RuntimeError(
            "Evenly spaced control ranks were not unique"
        )

    return ordered.iloc[indices].copy()


def main():
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not METRICS.is_file():
        raise FileNotFoundError(
            f"Missing control metrics: {METRICS}"
        )

    metrics = pd.read_csv(
        METRICS,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "sample",
        "call_rate",
        "heterozygosity_among_callable",
    }

    missing = required.difference(
        metrics.columns
    )

    if missing:
        raise RuntimeError(
            "Control metrics are missing: "
            + ", ".join(sorted(missing))
        )

    metrics = metrics.loc[
        metrics["sample"].str.match(
            r"^(?:SRR|ERR|DRR)\d+$"
        )
    ].copy()

    metrics["call_rate"] = pd.to_numeric(
        metrics["call_rate"],
        errors="coerce",
    )

    metrics[
        "heterozygosity_among_callable"
    ] = pd.to_numeric(
        metrics[
            "heterozygosity_among_callable"
        ],
        errors="coerce",
    )

    metrics = metrics.dropna(
        subset=[
            "call_rate",
            "heterozygosity_among_callable",
        ]
    )

    eligible = metrics.loc[
        metrics["call_rate"].ge(
            MIN_CALL_RATE
        )
    ].copy()

    if len(eligible) < TARGET_CONTROLS:
        raise RuntimeError(
            f"Only {len(eligible)} controls have "
            f"call rate >= {MIN_CALL_RATE}; "
            f"{TARGET_CONTROLS} are required"
        )

    selected = evenly_spaced(
        eligible,
        TARGET_CONTROLS,
    )

    selected[
        "provisional_het_stratum"
    ] = pd.qcut(
        selected[
            "heterozygosity_among_callable"
        ].rank(method="first"),
        q=3,
        labels=[
            "low",
            "middle",
            "high",
        ],
    ).astype(str)

    selected["selection_reason"] = (
        "evenly_spaced_across_optix_heterozygosity;"
        "call_rate_ge_0.65"
    )

    focal = pd.DataFrame(
        FOCAL_LIBRARIES
    )

    focal["role"] = "focal"
    focal["provisional_call_rate"] = ""
    focal["provisional_heterozygosity"] = ""
    focal["provisional_het_stratum"] = "focal"
    focal["selection_reason"] = (
        "locked_gynandromorph_experiment"
    )

    controls = pd.DataFrame(
        {
            "accession": selected["sample"],
            "library_id": selected["sample"],
            "side": "",
            "leg": "",
            "merge_group": selected["sample"],
            "role": "control",
            "provisional_call_rate": (
                selected["call_rate"]
            ),
            "provisional_heterozygosity": (
                selected[
                    "heterozygosity_among_callable"
                ]
            ),
            "provisional_het_stratum": (
                selected[
                    "provisional_het_stratum"
                ]
            ),
            "selection_reason": (
                selected["selection_reason"]
            ),
        }
    )

    manifest = pd.concat(
        [
            focal,
            controls,
        ],
        ignore_index=True,
    )

    manifest.insert(
        0,
        "array_id",
        np.arange(
            len(manifest),
            dtype=int,
        ),
    )

    if manifest[
        "accession"
    ].duplicated().any():
        raise RuntimeError(
            "Download manifest contains "
            "duplicate accessions"
        )

    if manifest[
        "library_id"
    ].duplicated().any():
        raise RuntimeError(
            "Download manifest contains "
            "duplicate library IDs"
        )

    selected.to_csv(
        OUT / "selected_control_run_review.tsv",
        sep="\t",
        index=False,
    )

    focal.to_csv(
        OUT / "focal_library_plan.tsv",
        sep="\t",
        index=False,
    )

    manifest.to_csv(
        OUT / "download_manifest.tsv",
        sep="\t",
        index=False,
    )

    analysis_sets = pd.DataFrame(
        [
            {
                "analysis_set": "merged_T2_T3",
                "left_sample": "TN_Left_merged",
                "right_sample": "TN_Right_merged",
                "purpose": (
                    "primary_genomewide_inference"
                ),
            },
            {
                "analysis_set": "T2_only",
                "left_sample": "TN_Left_T2",
                "right_sample": "TN_Right_T2",
                "purpose": (
                    "within_individual_"
                    "consistency_analysis"
                ),
            },
            {
                "analysis_set": "T3_only",
                "left_sample": "TN_Left_T3",
                "right_sample": "TN_Right_T3",
                "purpose": (
                    "within_individual_"
                    "consistency_analysis"
                ),
            },
        ]
    )

    analysis_sets.to_csv(
        OUT / "focal_analysis_sets.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(
        [
            {
                "accession": REFERENCE,
                "taxon": (
                    "Heliconius melpomene rosina"
                ),
                "assembly_name": (
                    "Heliconius_melpomene_rosina.546"
                ),
            }
        ]
    ).to_csv(
        OUT / "reference_manifest.tsv",
        sep="\t",
        index=False,
    )

    audit = pd.DataFrame(
        [
            (
                "focal_experiments",
                len(focal),
            ),
            (
                "focal_left_libraries",
                int(
                    focal["side"]
                    .eq("Left")
                    .sum()
                ),
            ),
            (
                "focal_right_libraries",
                int(
                    focal["side"]
                    .eq("Right")
                    .sum()
                ),
            ),
            (
                "focal_T2_libraries",
                int(
                    focal["leg"]
                    .eq("T2")
                    .sum()
                ),
            ),
            (
                "focal_T3_libraries",
                int(
                    focal["leg"]
                    .eq("T3")
                    .sum()
                ),
            ),
            (
                "candidate_controls",
                len(metrics),
            ),
            (
                "eligible_controls",
                len(eligible),
            ),
            (
                "selected_controls",
                len(selected),
            ),
            (
                (
                    "selected_min_provisional_"
                    "heterozygosity"
                ),
                selected[
                    "heterozygosity_among_callable"
                ].min(),
            ),
            (
                (
                    "selected_max_provisional_"
                    "heterozygosity"
                ),
                selected[
                    "heterozygosity_among_callable"
                ].max(),
            ),
            (
                "download_accessions",
                len(manifest),
            ),
            (
                "reference_accession",
                REFERENCE,
            ),
            (
                "primary_focal_data",
                "T2_plus_T3_merged_by_side",
            ),
            (
                "sensitivity_focal_data",
                "T2_and_T3_analysed_separately",
            ),
            (
                "control_selection_uses_focal_data",
                False,
            ),
            (
                "optix_heterozygosity_is_provisional",
                True,
            ),
            (
                "NCBI_EUtilities_used",
                False,
            ),
        ],
        columns=[
            "item",
            "value",
        ],
    )

    audit.to_csv(
        OUT / "download_manifest.audit.tsv",
        sep="\t",
        index=False,
    )

    print(
        "GENOME-WIDE DOWNLOAD MANIFEST LOCKED"
    )

    print(
        audit.to_string(
            index=False
        )
    )

    print(
        "\nFOCAL LIBRARIES"
    )

    print(
        focal[
            [
                "accession",
                "library_id",
                "side",
                "leg",
                "merge_group",
            ]
        ].to_string(
            index=False
        )
    )

    print(
        "\nSELECTED CONTROLS"
    )

    print(
        controls[
            [
                "accession",
                "provisional_call_rate",
                "provisional_heterozygosity",
                "provisional_het_stratum",
            ]
        ].to_string(
            index=False
        )
    )

    print(
        "\nFOCAL ANALYSIS SETS"
    )

    print(
        analysis_sets.to_string(
            index=False
        )
    )


if __name__ == "__main__":
    try:
        main()

    except Exception as exc:
        print(
            f"ERROR: {exc}",
            file=sys.stderr,
        )
        raise
