#!/usr/bin/env python3

from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

MIXTURE = (
    PROJECT
    / "analysis/original_502/reciprocal_L_R_mixture/"
    "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

CONTROLS = (
    PROJECT
    / "data/derived/region_genotypes/"
    "optix_51sample.genotypes.tsv.gz"
)

OUTPUT = (
    PROJECT
    / "analysis/original_502/checkpoints/"
    "three_vs_four_founder_input_audit.tsv"
)


def parse_bool(series):
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "t", "1", "yes"})
    )


def numeric(series):
    return pd.to_numeric(
        series.replace(
            {
                "": np.nan,
                ".": np.nan,
                "NA": np.nan,
                "nan": np.nan,
            }
        ),
        errors="coerce",
    )


def canonical_gt(value):
    value = str(value).strip().replace("|", "/")

    if value == "1/0":
        value = "0/1"

    if value in {"0/0", "0/1", "1/1"}:
        return value

    return "unavailable"


def report_estimate_column(data, column):
    values = numeric(data[column])
    finite = values[np.isfinite(values)]

    print(f"\n{column}")

    if finite.empty:
        print("  no finite values")
        return

    print(f"  finite rows: {len(finite)}")
    print(f"  unique values: {finite.nunique()}")
    print(f"  minimum: {finite.min():.12g}")
    print(f"  median: {finite.median():.12g}")
    print(f"  maximum: {finite.max():.12g}")

    print("  most frequent values:")

    counts = finite.value_counts().head(10)

    for value, count in counts.items():
        print(f"    {value:.12g}\t{count}")


def main():
    print("THREE- VERSUS FOUR-FOUNDER INPUT AUDIT")
    print("Autosomes only")
    print("Two-founder model excluded")
    print("T2/T3 excluded as biological samples")

    if not MIXTURE.exists():
        raise FileNotFoundError(MIXTURE)

    if not CONTROLS.exists():
        raise FileNotFoundError(CONTROLS)

    mixture = pd.read_csv(
        MIXTURE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        low_memory=False,
    )

    print("\nMIXTURE TABLE")
    print(f"path: {MIXTURE}")
    print(f"rows: {len(mixture)}")
    print(f"columns: {len(mixture.columns)}")

    for column in mixture.columns:
        print(f"  {column}")

    required_mixture = {
        "chrom",
        "pos",
        "ref",
        "alt",
        "left_ref_depth",
        "left_alt_depth",
        "left_DP",
        "right_ref_depth",
        "right_alt_depth",
        "right_DP",
        "mixture_left_GT",
        "mixture_right_GT",
        "pair_state",
        "pair_posterior",
        "posterior95",
        "estimated_c_left_from_right",
        "estimated_c_right_from_left",
    }

    missing_mixture = (
        required_mixture
        - set(mixture.columns)
    )

    if missing_mixture:
        raise RuntimeError(
            "Mixture table is missing required columns: "
            f"{sorted(missing_mixture)}"
        )

    posterior = numeric(
        mixture["pair_posterior"]
    )

    confident = (
        parse_bool(mixture["posterior95"])
        & posterior.ge(0.95)
    )

    print("\nMIXTURE CONFIDENCE")
    print(
        f"posterior >=0.95 rows: "
        f"{int(confident.sum())}"
    )
    print(
        f"non-confident rows: "
        f"{int((~confident).sum())}"
    )

    report_estimate_column(
        mixture,
        "estimated_c_left_from_right",
    )

    report_estimate_column(
        mixture,
        "estimated_c_right_from_left",
    )

    mixture["left_gt_clean"] = (
        mixture["mixture_left_GT"]
        .map(canonical_gt)
    )

    mixture["right_gt_clean"] = (
        mixture["mixture_right_GT"]
        .map(canonical_gt)
    )

    confident_pairs = mixture.loc[
        confident,
        [
            "left_gt_clean",
            "right_gt_clean",
        ],
    ].copy()

    genotype_order = [
        "0/0",
        "0/1",
        "1/1",
        "unavailable",
    ]

    pair_matrix = pd.crosstab(
        confident_pairs["left_gt_clean"],
        confident_pairs["right_gt_clean"],
        dropna=False,
    ).reindex(
        index=genotype_order,
        columns=genotype_order,
        fill_value=0,
    )

    print(
        "\nPOSTERIOR>=0.95 "
        "LEFT/RIGHT GENOTYPE-PAIR MATRIX"
    )
    print(pair_matrix.to_string())

    opposing = (
        (
            confident_pairs["left_gt_clean"].eq("0/0")
            & confident_pairs[
                "right_gt_clean"
            ].eq("1/1")
        )
        |
        (
            confident_pairs["left_gt_clean"].eq("1/1")
            & confident_pairs[
                "right_gt_clean"
            ].eq("0/0")
        )
    )

    opposing_count = int(opposing.sum())

    print("\nOPPOSING HOMOZYGOTES")
    print(
        "Posterior>=0.95 0/0|1/1 or "
        f"1/1|0/0 calls: {opposing_count}"
    )
    print(
        "These are descriptive only and will not "
        "be used as fixed truth in model fitting."
    )

    print("\nFOCAL DEPTH SUMMARY")

    for column in [
        "left_DP",
        "right_DP",
    ]:
        values = numeric(mixture[column])
        finite = values[np.isfinite(values)]

        print(
            f"{column}: "
            f"n={len(finite)} "
            f"min={finite.min():.6g} "
            f"median={finite.median():.6g} "
            f"mean={finite.mean():.6g} "
            f"max={finite.max():.6g}"
        )

    control_header = pd.read_csv(
        CONTROLS,
        sep="\t",
        compression="infer",
        nrows=0,
    )

    control_columns = list(
        control_header.columns
    )

    print("\nCONTROL TABLE")
    print(f"path: {CONTROLS}")
    print(
        f"columns: {len(control_columns)}"
    )

    for column in control_columns:
        print(f"  {column}")

    required_control = {
        "sample",
        "chrom",
        "pos",
        "ref",
        "alt",
        "GT",
        "DP",
        "GQ",
    }

    missing_control = (
        required_control
        - set(control_columns)
    )

    if missing_control:
        raise RuntimeError(
            "Control table is missing required columns: "
            f"{sorted(missing_control)}"
        )

    normalized = {
        column: (
            column.lower()
            .replace("_", "")
            .replace("-", "")
        )
        for column in control_columns
    }

    count_like_columns = [
        column
        for column, name in normalized.items()
        if (
            name in {
                "ad",
                "refdepth",
                "altdepth",
                "refcount",
                "altcount",
                "rd",
            }
            or "alleledepth" in name
            or "refdepth" in name
            or "altdepth" in name
            or "refcount" in name
            or "altcount" in name
        )
    ]

    print("\nCONTROL READ-COUNT COLUMNS")
    if count_like_columns:
        for column in count_like_columns:
            print(f"  {column}")
    else:
        print("  none detected")

    selected_columns = list(
        dict.fromkeys(
            [
                "sample",
                "chrom",
                "pos",
                "ref",
                "alt",
                "GT",
                "DP",
                "GQ",
            ]
            + count_like_columns
        )
    )

    sample_names = set()
    sites = set()
    genotype_counts = Counter()
    total_rows = 0
    qc_rows = 0

    for chunk in pd.read_csv(
        CONTROLS,
        sep="\t",
        compression="infer",
        dtype=str,
        keep_default_na=False,
        usecols=selected_columns,
        chunksize=250_000,
        low_memory=False,
    ):
        total_rows += len(chunk)

        sample_names.update(
            chunk["sample"].astype(str)
        )

        site_keys = zip(
            chunk["chrom"].astype(str),
            chunk["pos"].astype(str),
            chunk["ref"].astype(str),
            chunk["alt"].astype(str),
        )

        sites.update(site_keys)

        clean_gt = (
            chunk["GT"]
            .map(canonical_gt)
        )

        genotype_counts.update(
            clean_gt.tolist()
        )

        depth = numeric(chunk["DP"])
        quality = numeric(chunk["GQ"])

        qc_pass = (
            depth.ge(2)
            & quality.ge(20)
        )

        qc_rows += int(qc_pass.sum())

    print("\nCONTROL TABLE SUMMARY")
    print(f"rows: {total_rows}")
    print(f"samples: {len(sample_names)}")
    print(f"unique sites: {len(sites)}")
    print(f"DP>=2/GQ>=20 rows: {qc_rows}")

    print("genotype rows:")
    for genotype in genotype_order:
        print(
            f"  {genotype}: "
            f"{genotype_counts[genotype]}"
        )

    forbidden_samples = {
        "TN2L_S1",
        "TN3L_S3",
        "TN2R_S2",
        "TN3R_S4",
    }

    present_forbidden = sorted(
        sample_names & forbidden_samples
    )

    print("\nT2/T3 CONTROL-TABLE AUDIT")
    print(
        "present before explicit exclusion: "
        + (
            ", ".join(present_forbidden)
            if present_forbidden
            else "none"
        )
    )
    print(
        "Founder-model pipeline will exclude "
        "all four explicitly."
    )

    left_estimates = numeric(
        mixture[
            "estimated_c_left_from_right"
        ]
    )

    right_estimates = numeric(
        mixture[
            "estimated_c_right_from_left"
        ]
    )

    audit_rows = [
        {
            "metric": "mixture_rows",
            "value": len(mixture),
        },
        {
            "metric": "posterior95_rows",
            "value": int(confident.sum()),
        },
        {
            "metric": (
                "posterior95_opposing_"
                "homozygote_calls"
            ),
            "value": opposing_count,
        },
        {
            "metric": (
                "unique_c_left_from_right"
            ),
            "value": int(
                left_estimates.nunique(
                    dropna=True
                )
            ),
        },
        {
            "metric": (
                "median_c_left_from_right"
            ),
            "value": left_estimates.median(),
        },
        {
            "metric": (
                "unique_c_right_from_left"
            ),
            "value": int(
                right_estimates.nunique(
                    dropna=True
                )
            ),
        },
        {
            "metric": (
                "median_c_right_from_left"
            ),
            "value": right_estimates.median(),
        },
        {
            "metric": "control_rows",
            "value": total_rows,
        },
        {
            "metric": "control_samples",
            "value": len(sample_names),
        },
        {
            "metric": "control_sites",
            "value": len(sites),
        },
        {
            "metric": (
                "control_count_columns"
            ),
            "value": (
                ",".join(count_like_columns)
                if count_like_columns
                else "none"
            ),
        },
    ]

    OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    pd.DataFrame(
        audit_rows
    ).to_csv(
        OUTPUT,
        sep="\t",
        index=False,
    )

    print("\nAUDIT COMPLETE")
    print(f"Output: {OUTPUT}")


if __name__ == "__main__":
    main()
