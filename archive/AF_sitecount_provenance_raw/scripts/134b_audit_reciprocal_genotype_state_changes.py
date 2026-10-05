#!/usr/bin/env python3
from pathlib import Path
import hashlib

import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

INPUT = (
    PROJECT
    / "analysis/original_502/reciprocal_L_R_mixture"
    / "optix_reciprocal_L_R_genotype_pairs.tsv.gz"
)

EXPECTED_SHA256 = (
    "8dd2ebecdcd039ca643e9b68c122407bf"
    "a8f98745df2aac1e00a61901e0aa9a8"
)

OUTDIR = INPUT.parent

SUMMARY_OUT = (
    OUTDIR
    / "optix_reciprocal_L_R_genotype_change_summary.tsv"
)

PAIR_TRANSITIONS_OUT = (
    OUTDIR
    / "optix_reciprocal_L_R_high_confidence_pair_transitions.tsv"
)

SIDE_TRANSITIONS_OUT = (
    OUTDIR
    / "optix_reciprocal_L_R_high_confidence_side_transitions.tsv"
)


def sha256(path):

    digest = hashlib.sha256()

    with open(path, "rb") as stream:

        for block in iter(
            lambda: stream.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def canonical_gt(value):

    value = str(value).strip()

    if value in {
        "",
        ".",
        "./.",
        ".|.",
        "nan",
        "None",
    }:
        return "./."

    value = value.replace(
        "|",
        "/",
    )

    alleles = value.split("/")

    if (
        len(alleles) == 2
        and all(
            allele in {
                "0",
                "1",
            }
            for allele in alleles
        )
    ):
        alleles = sorted(
            alleles,
            key=int,
        )

        return "/".join(
            alleles
        )

    return value


def main():

    if not INPUT.is_file():

        raise RuntimeError(
            f"Missing input: {INPUT}"
        )

    observed_hash = sha256(
        INPUT
    )

    if observed_hash != EXPECTED_SHA256:

        raise RuntimeError(
            "Input checksum mismatch. "
            f"Expected {EXPECTED_SHA256}; "
            f"observed {observed_hash}."
        )

    data = pd.read_csv(
        INPUT,
        sep="\t",
        low_memory=False,
    )

    data["left_original_canonical"] = (
        data["left_original_GT"]
        .map(
            canonical_gt
        )
    )

    data["right_original_canonical"] = (
        data["right_original_GT"]
        .map(
            canonical_gt
        )
    )

    data["left_mixture_canonical"] = (
        data["mixture_left_GT"]
        .map(
            canonical_gt
        )
    )

    data["right_mixture_canonical"] = (
        data["mixture_right_GT"]
        .map(
            canonical_gt
        )
    )

    data["original_pair_canonical"] = (
        data["left_original_canonical"]
        + "|"
        + data["right_original_canonical"]
    )

    data["mixture_pair_canonical"] = (
        data["left_mixture_canonical"]
        + "|"
        + data["right_mixture_canonical"]
    )

    data["left_changed_canonical"] = (
        data["left_original_canonical"]
        != data["left_mixture_canonical"]
    )

    data["right_changed_canonical"] = (
        data["right_original_canonical"]
        != data["right_mixture_canonical"]
    )

    data["pair_changed_canonical"] = (
        data["left_changed_canonical"]
        | data["right_changed_canonical"]
    )

    posterior95 = (
        data["posterior95"]
        .astype(str)
        .str.lower()
        .isin(
            {
                "true",
                "t",
                "1",
            }
        )
    )

    confident = data.loc[
        posterior95
    ].copy()

    changed = confident.loc[
        confident[
            "pair_changed_canonical"
        ]
    ].copy()

    left_only = (
        confident[
            "left_changed_canonical"
        ]
        & ~confident[
            "right_changed_canonical"
        ]
    )

    right_only = (
        ~confident[
            "left_changed_canonical"
        ]
        & confident[
            "right_changed_canonical"
        ]
    )

    both_changed = (
        confident[
            "left_changed_canonical"
        ]
        & confident[
            "right_changed_canonical"
        ]
    )

    summary = pd.DataFrame(
        [
            (
                "all_depth_passing_sites",
                len(data),
            ),
            (
                "posterior95_sites",
                len(confident),
            ),
            (
                "below_posterior95_sites",
                len(data) - len(confident),
            ),
            (
                "posterior95_unchanged_pairs",
                int(
                    (
                        ~confident[
                            "pair_changed_canonical"
                        ]
                    ).sum()
                ),
            ),
            (
                "posterior95_changed_pairs",
                len(changed),
            ),
            (
                "posterior95_left_only_changed",
                int(
                    left_only.sum()
                ),
            ),
            (
                "posterior95_right_only_changed",
                int(
                    right_only.sum()
                ),
            ),
            (
                "posterior95_both_sides_changed",
                int(
                    both_changed.sum()
                ),
            ),
            (
                "posterior95_changed_fraction",
                (
                    len(changed)
                    / len(confident)
                    if len(confident)
                    else float("nan")
                ),
            ),
        ],
        columns=[
            "statistic",
            "value",
        ],
    )

    pair_transitions = (
        confident.groupby(
            [
                "original_pair_canonical",
                "mixture_pair_canonical",
            ],
            dropna=False,
        )
        .size()
        .reset_index(
            name="sites"
        )
        .sort_values(
            [
                "sites",
                "original_pair_canonical",
                "mixture_pair_canonical",
            ],
            ascending=[
                False,
                True,
                True,
            ],
        )
    )

    side_rows = []

    for side in [
        "left",
        "right",
    ]:

        side_table = (
            confident.groupby(
                [
                    f"{side}_original_canonical",
                    f"{side}_mixture_canonical",
                ],
                dropna=False,
            )
            .size()
            .reset_index(
                name="sites"
            )
        )

        side_table.columns = [
            "original_GT",
            "mixture_GT",
            "sites",
        ]

        side_table.insert(
            0,
            "side",
            side,
        )

        side_rows.append(
            side_table
        )

    side_transitions = pd.concat(
        side_rows,
        ignore_index=True,
    ).sort_values(
        [
            "side",
            "sites",
            "original_GT",
            "mixture_GT",
        ],
        ascending=[
            True,
            False,
            True,
            True,
        ],
    )

    summary.to_csv(
        SUMMARY_OUT,
        sep="\t",
        index=False,
    )

    pair_transitions.to_csv(
        PAIR_TRANSITIONS_OUT,
        sep="\t",
        index=False,
    )

    side_transitions.to_csv(
        SIDE_TRANSITIONS_OUT,
        sep="\t",
        index=False,
    )

    print(
        "RECIPROCAL GENOTYPE-STATE "
        "CHANGE AUDIT"
    )

    print()

    print(
        summary.to_string(
            index=False
        )
    )

    print(
        "\nCHANGED HIGH-CONFIDENCE "
        "PAIR TRANSITIONS"
    )

    print(
        pair_transitions.loc[
            pair_transitions[
                "original_pair_canonical"
            ]
            != pair_transitions[
                "mixture_pair_canonical"
            ]
        ]
        .head(30)
        .to_string(
            index=False
        )
    )

    print(
        f"\nSummary: {SUMMARY_OUT}"
    )

    print(
        f"Pair transitions: "
        f"{PAIR_TRANSITIONS_OUT}"
    )

    print(
        f"Side transitions: "
        f"{SIDE_TRANSITIONS_OUT}"
    )


if __name__ == "__main__":
    main()
