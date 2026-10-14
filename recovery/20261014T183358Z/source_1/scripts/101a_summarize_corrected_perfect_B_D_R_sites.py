#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

project = Path("/scratch/martinlab/joe/GynandroCalling")
checkpoint_dir = (
    project
    / "analysis"
    / "drop_problematic_reanalysis"
    / "checkpoints"
)

modules = {
    "Band": "longBand",
    "Dennis": "longDennis",
    "Ray": "longRay",
}

def as_bool(series):
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "t", "1", "yes"})
    )

rows = []
perfect_rows = []

for module, prefix in modules.items():

    score_file = (
        checkpoint_dir
        / (
            f"{prefix}_original502_expected_pattern_sites."
            "WGS_vs_WGSplusCS.tsv"
        )
    )

    paired_file = (
        checkpoint_dir
        / (
            f"{prefix}_original502_allSNP_scores.filtered."
            "WGS_vs_WGSplusCS.minimum_across_schemes.tsv"
        )
    )

    if not score_file.exists():
        print(f"MISSING: {score_file}")
        continue

    scores = pd.read_csv(score_file, sep="\t", low_memory=False)

    scores["balanced_accuracy"] = pd.to_numeric(
        scores["balanced_accuracy"],
        errors="coerce",
    )

    scores["misclassified_n"] = pd.to_numeric(
        scores["misclassified_n"],
        errors="coerce",
    )

    for cohort in ("WGS_only", "WGS_plus_CS"):
        for scheme in ("called_GT", "DP2_GQ20"):

            subset = scores.loc[
                (scores["analysis_cohort"] == cohort)
                & (scores["scheme"] == scheme)
            ].copy()

            if "expected_pattern_eligible" in subset.columns:
                subset = subset.loc[
                    as_bool(subset["expected_pattern_eligible"])
                ].copy()

            exact = subset.loc[
                subset["balanced_accuracy"].ge(
                    1.0 - 1e-12
                )
            ].copy()

            best_ba = (
                subset["balanced_accuracy"].max()
                if len(subset)
                else float("nan")
            )

            minimum_errors = (
                subset.loc[
                    subset["balanced_accuracy"].eq(best_ba),
                    "misclassified_n",
                ].min()
                if len(subset)
                else float("nan")
            )

            rows.append(
                {
                    "module": module,
                    "cohort": cohort,
                    "scheme": scheme,
                    "eligible_sites": len(subset),
                    "perfect_sites": len(exact),
                    "maximum_BA": best_ba,
                    "minimum_errors_at_maximum_BA": (
                        minimum_errors
                    ),
                }
            )

            for _, row in exact.iterrows():
                perfect_rows.append(
                    {
                        "module": module,
                        "cohort": cohort,
                        "scheme": scheme,
                        "chrom": row["chrom"],
                        "pos": row["pos"],
                        "ref": row["ref"],
                        "alt": row["alt"],
                        "balanced_accuracy": (
                            row["balanced_accuracy"]
                        ),
                        "misclassified_n": (
                            row["misclassified_n"]
                        ),
                    }
                )

    if paired_file.exists():

        paired = pd.read_csv(
            paired_file,
            sep="\t",
            low_memory=False,
        )

        paired["minimum_BA_across_schemes"] = pd.to_numeric(
            paired["minimum_BA_across_schemes"],
            errors="coerce",
        )

        for cohort in ("WGS_only", "WGS_plus_CS"):

            subset = paired.loc[
                paired["analysis_cohort"] == cohort
            ].copy()

            if "all_schemes_expected_pattern_eligible" in subset.columns:
                subset = subset.loc[
                    as_bool(
                        subset[
                            "all_schemes_expected_pattern_eligible"
                        ]
                    )
                ].copy()

            exact = subset.loc[
                subset["minimum_BA_across_schemes"].ge(
                    1.0 - 1e-12
                )
            ]

            rows.append(
                {
                    "module": module,
                    "cohort": cohort,
                    "scheme": "both_schemes",
                    "eligible_sites": len(subset),
                    "perfect_sites": len(exact),
                    "maximum_BA": (
                        subset[
                            "minimum_BA_across_schemes"
                        ].max()
                        if len(subset)
                        else float("nan")
                    ),
                    "minimum_errors_at_maximum_BA": "",
                }
            )

summary = pd.DataFrame(rows)
perfect = pd.DataFrame(perfect_rows)

summary_file = (
    checkpoint_dir
    / "corrected_B_D_R_perfect_site_counts.tsv"
)

perfect_file = (
    checkpoint_dir
    / "corrected_B_D_R_perfect_sites.tsv"
)

summary.to_csv(summary_file, sep="\t", index=False)
perfect.to_csv(perfect_file, sep="\t", index=False)

print("CORRECTED PERFECT-SITE COUNTS")
print(summary.to_string(index=False))

print()
print("PERFECT SITE COORDINATES")

if len(perfect):
    print(perfect.to_string(index=False))
else:
    print("None")

print()
print(f"Summary: {summary_file}")
print(f"Perfect sites: {perfect_file}")

