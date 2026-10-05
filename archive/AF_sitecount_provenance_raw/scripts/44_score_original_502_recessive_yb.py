#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

from scipy.stats import fisher_exact


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

VCF = (
    PROJECT
    / "external/krista_genotype_plots_original/Genotype_Plots"
    / "Ivory_Hmr15-1627000-1637000_biallelic.vcf"
)

METADATA = (
    PROJECT
    / "data/derived/region_panel_metadata.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis/original_502/checkpoints"
)

SCORES = (
    OUTDIR
    / "ivory_original_502_comparative_scores.tsv"
)

CANDIDATES = (
    OUTDIR
    / "ivory_original_502_inclusive_focal_candidates.tsv"
)

COHORTS = (
    OUTDIR
    / "ivory_original_502_scoring_cohorts.tsv"
)

TRAITS = {
    "Yb_dorsal": "Yb-dorsal",
    "Yb_ventral": "Yb-ventral",
}

SCHEMES = {
    "called_GT": (
        lambda call:
        call[0] is not None
    ),
    "DP2_GQ20": (
        lambda call:
        call[0] is not None
        and call[1] is not None
        and call[1] >= 2
        and call[2] is not None
        and call[2] >= 20
    ),
}

FOCAL = {
    "TN2L_S1",
    "TN3L_S3",
    "TN2R_S2",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
}


def number(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None

    if math.isfinite(result):
        return result

    return None


def boolean(value):
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


def canonical_dosage(gt):
    if not gt or "." in gt:
        return None

    alleles = (
        gt.replace("|", "/")
        .split("/")
    )

    if len(alleles) != 2:
        return None

    if any(
        allele not in {"0", "1"}
        for allele in alleles
    ):
        return None

    return sum(
        int(allele)
        for allele in alleles
    )


def canonical_gt(dosage):
    return {
        0: "0/0",
        1: "0/1",
        2: "1/1",
    }.get(
        dosage,
        "./.",
    )


def bh_adjust(values):
    result = [
        None
    ] * len(values)

    valid = [
        (
            index,
            value,
        )
        for index, value in enumerate(values)
        if value is not None
    ]

    valid.sort(
        key=lambda item: item[1]
    )

    number_of_tests = len(valid)
    running = 1.0

    for rank_from_end in range(
        number_of_tests - 1,
        -1,
        -1,
    ):
        original_index, pvalue = (
            valid[rank_from_end]
        )

        rank = (
            rank_from_end + 1
        )

        running = min(
            running,
            (
                pvalue
                * number_of_tests
                / rank
            ),
        )

        result[original_index] = min(
            1.0,
            running,
        )

    return result


def calculate_mcc(
    true_positive,
    true_negative,
    false_positive,
    false_negative,
):
    denominator = (
        (true_positive + false_positive)
        * (true_positive + false_negative)
        * (true_negative + false_positive)
        * (true_negative + false_negative)
    )

    if denominator <= 0:
        return None

    return (
        (
            true_positive
            * true_negative
        )
        - (
            false_positive
            * false_negative
        )
    ) / math.sqrt(denominator)


def read_metadata():
    rows = {}

    with METADATA.open(
        newline=""
    ) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            sample = (
                row.get("Accession")
                or ""
            ).strip()

            if sample:
                rows[sample] = row

    return rows


def read_vcf():
    samples = []
    variants = []
    calls_by_variant = []

    with VCF.open() as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                samples = (
                    line.rstrip("\n")
                    .split("\t")[9:]
                )

                continue

            if (
                line.startswith("#")
                or not line.strip()
            ):
                continue

            fields = (
                line.rstrip("\n")
                .split("\t")
            )

            chrom = fields[0]
            position = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            if len(ref) != 1:
                continue

            if (
                len(alt) != 1
                or "," in alt
            ):
                continue

            format_names = (
                fields[8]
                .split(":")
            )

            format_index = {
                name: index
                for index, name
                in enumerate(format_names)
            }

            if "GT" not in format_index:
                continue

            calls = []

            for sample_field in fields[9:]:
                values = (
                    sample_field
                    .split(":")
                )

                gt_index = (
                    format_index["GT"]
                )

                gt = (
                    values[gt_index]
                    if gt_index < len(values)
                    else "."
                )

                dosage = (
                    canonical_dosage(gt)
                )

                depth = None

                if "DP" in format_index:
                    depth_index = (
                        format_index["DP"]
                    )

                    if depth_index < len(values):
                        depth = number(
                            values[depth_index]
                        )

                genotype_quality = None

                if "GQ" in format_index:
                    quality_index = (
                        format_index["GQ"]
                    )

                    if quality_index < len(values):
                        genotype_quality = number(
                            values[quality_index]
                        )

                allelic_depth = ""

                if "AD" in format_index:
                    ad_index = (
                        format_index["AD"]
                    )

                    if ad_index < len(values):
                        allelic_depth = (
                            values[ad_index]
                        )

                calls.append(
                    (
                        dosage,
                        depth,
                        genotype_quality,
                        allelic_depth,
                    )
                )

            variants.append(
                (
                    chrom,
                    position,
                    ref,
                    alt,
                )
            )

            calls_by_variant.append(
                calls
            )

    if not samples:
        raise RuntimeError(
            "VCF sample header was not found"
        )

    return (
        samples,
        variants,
        calls_by_variant,
    )


def phenotype(
    row,
    column,
):
    value = number(
        row.get(column)
    )

    if value in {
        0.0,
        1.0,
    }:
        return int(value)

    return None


def score_orientation(
    records,
    associated_homozygote,
):
    true_positive = 0
    true_negative = 0
    false_positive = 0
    false_negative = 0

    for observed, dosage in records:
        predicted = int(
            dosage
            == associated_homozygote
        )

        if (
            observed == 1
            and predicted == 1
        ):
            true_positive += 1

        elif (
            observed == 0
            and predicted == 0
        ):
            true_negative += 1

        elif (
            observed == 0
            and predicted == 1
        ):
            false_positive += 1

        else:
            false_negative += 1

    sensitivity = (
        true_positive
        / (
            true_positive
            + false_negative
        )
        if (
            true_positive
            + false_negative
        )
        else None
    )

    specificity = (
        true_negative
        / (
            true_negative
            + false_positive
        )
        if (
            true_negative
            + false_positive
        )
        else None
    )

    balanced_accuracy = (
        (
            sensitivity
            + specificity
        )
        / 2
        if (
            sensitivity is not None
            and specificity is not None
        )
        else None
    )

    contingency_table = [
        [
            true_positive,
            false_negative,
        ],
        [
            false_positive,
            true_negative,
        ],
    ]

    fisher_p = (
        fisher_exact(
            contingency_table
        )[1]
        if all(
            sum(row) > 0
            for row in contingency_table
        )
        else None
    )

    return {
        "tp": true_positive,
        "tn": true_negative,
        "fp": false_positive,
        "fn": false_negative,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "balanced_accuracy": balanced_accuracy,
        "MCC": calculate_mcc(
            true_positive,
            true_negative,
            false_positive,
            false_negative,
        ),
        "fisher_p": fisher_p,
    }


def main():
    OUTDIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    required = [
        VCF,
        METADATA,
    ]

    missing = [
        str(path)
        for path in required
        if not path.is_file()
    ]

    if missing:
        print(
            "ERROR: required inputs are missing:"
        )

        for path in missing:
            print(path)

        return

    metadata = read_metadata()

    (
        samples,
        variants,
        all_calls,
    ) = read_vcf()

    sample_index = {
        sample: index
        for index, sample
        in enumerate(samples)
    }

    cohorts = {}
    cohort_rows = []

    for trait, column in TRAITS.items():
        selected_samples = []

        phenotype_counts = {
            0: 0,
            1: 0,
        }

        excluded = defaultdict(int)

        for sample in samples:
            row = metadata.get(sample)

            if row is None:
                excluded[
                    "no_metadata"
                ] += 1

                continue

            if sample in FOCAL:
                excluded[
                    "focal"
                ] += 1

                continue

            taxon = (
                row.get("Taxon Name")
                or ""
            ).lower()

            if "piano keys" in taxon:
                excluded[
                    "piano_keys"
                ] += 1

                continue

            state = phenotype(
                row,
                column,
            )

            if state is None:
                excluded[
                    "missing_phenotype"
                ] += 1

                continue

            selected_samples.append(
                (
                    sample,
                    state,
                )
            )

            phenotype_counts[
                state
            ] += 1

        cohorts[trait] = (
            selected_samples
        )

        cohort_rows.append(
            {
                "trait": trait,
                "VCF_samples": len(samples),
                "comparative_scoring_samples": (
                    len(selected_samples)
                ),
                "phenotype_0_n": (
                    phenotype_counts[0]
                ),
                "phenotype_1_n": (
                    phenotype_counts[1]
                ),
                "focal_excluded": (
                    excluded["focal"]
                ),
                "piano_keys_excluded": (
                    excluded["piano_keys"]
                ),
                "missing_phenotype": (
                    excluded[
                        "missing_phenotype"
                    ]
                ),
                "no_metadata_match": (
                    excluded["no_metadata"]
                ),
            }
        )

    score_rows = []

    for trait, cohort in cohorts.items():
        total_0 = sum(
            state == 0
            for _, state in cohort
        )

        total_1 = sum(
            state == 1
            for _, state in cohort
        )

        for variant, calls in zip(
            variants,
            all_calls,
        ):
            (
                chrom,
                position,
                ref_allele,
                alt_allele,
            ) = variant

            for (
                scheme,
                passes_filter,
            ) in SCHEMES.items():
                records = []

                called_0 = 0
                called_1 = 0
                missing_0 = 0
                missing_1 = 0

                for sample, state in cohort:
                    call = calls[
                        sample_index[sample]
                    ]

                    if passes_filter(call):
                        records.append(
                            (
                                state,
                                call[0],
                            )
                        )

                        if state == 0:
                            called_0 += 1
                        else:
                            called_1 += 1

                    elif state == 0:
                        missing_0 += 1

                    else:
                        missing_1 += 1

                alt_result = (
                    score_orientation(
                        records,
                        2,
                    )
                )

                ref_result = (
                    score_orientation(
                        records,
                        0,
                    )
                )

                alt_accuracy = (
                    alt_result[
                        "balanced_accuracy"
                    ]
                )

                ref_accuracy = (
                    ref_result[
                        "balanced_accuracy"
                    ]
                )

                if (
                    (
                        alt_accuracy
                        if alt_accuracy is not None
                        else -1
                    )
                    >
                    (
                        ref_accuracy
                        if ref_accuracy is not None
                        else -1
                    )
                ):
                    selected_result = (
                        alt_result
                    )

                    orientation = (
                        "ALT_homozygote_"
                        "predicts_presence"
                    )

                    associated_homozygote = 2

                else:
                    selected_result = (
                        ref_result
                    )

                    orientation = (
                        "REF_homozygote_"
                        "predicts_presence"
                    )

                    associated_homozygote = 0

                predicted_1 = sum(
                    dosage
                    == associated_homozygote
                    for _, dosage in records
                )

                predicted_0 = (
                    len(records)
                    - predicted_1
                )

                call_rate_0 = (
                    called_0 / total_0
                    if total_0
                    else 0
                )

                call_rate_1 = (
                    called_1 / total_1
                    if total_1
                    else 0
                )

                missingness_p = (
                    fisher_exact(
                        [
                            [
                                called_0,
                                missing_0,
                            ],
                            [
                                called_1,
                                missing_1,
                            ],
                        ]
                    )[1]
                )

                callability_eligible = (
                    call_rate_0 >= 0.50
                    and call_rate_1 >= 0.50
                )

                score_eligible = (
                    callability_eligible
                    and called_0 >= 10
                    and called_1 >= 10
                    and predicted_0 >= 5
                    and predicted_1 >= 5
                )

                score_rows.append(
                    {
                        "trait": trait,
                        "scheme": scheme,
                        "chrom": chrom,
                        "pos": position,
                        "ref": ref_allele,
                        "alt": alt_allele,
                        "comparative_cohort_n": (
                            len(cohort)
                        ),
                        "phenotype_0_total": (
                            total_0
                        ),
                        "phenotype_1_total": (
                            total_1
                        ),
                        "phenotype_0_called": (
                            called_0
                        ),
                        "phenotype_1_called": (
                            called_1
                        ),
                        "phenotype_0_call_rate": (
                            call_rate_0
                        ),
                        "phenotype_1_call_rate": (
                            call_rate_1
                        ),
                        "missingness_fisher_p": (
                            missingness_p
                        ),
                        "selected_orientation": (
                            orientation
                        ),
                        "matches_reference_expectation": (
                            orientation.startswith(
                                "REF_"
                            )
                        ),
                        "predicted_0": predicted_0,
                        "predicted_1": predicted_1,
                        **selected_result,
                        "misclassified_n": (
                            selected_result["fp"]
                            + selected_result["fn"]
                        ),
                        "callability_eligible": (
                            callability_eligible
                        ),
                        "score_eligible": (
                            score_eligible
                        ),
                        "BH_q": None,
                        "missingness_BH_q": None,
                    }
                )

    for trait in TRAITS:
        for scheme in SCHEMES:
            indices = [
                index
                for index, row
                in enumerate(score_rows)
                if (
                    row["trait"] == trait
                    and row["scheme"] == scheme
                )
            ]

            qvalues = bh_adjust(
                [
                    score_rows[index][
                        "fisher_p"
                    ]
                    for index in indices
                ]
            )

            missingness_qvalues = (
                bh_adjust(
                    [
                        score_rows[index][
                            "missingness_fisher_p"
                        ]
                        for index in indices
                    ]
                )
            )

            for (
                index,
                qvalue,
                missingness_qvalue,
            ) in zip(
                indices,
                qvalues,
                missingness_qvalues,
            ):
                score_rows[index][
                    "BH_q"
                ] = qvalue

                score_rows[index][
                    "missingness_BH_q"
                ] = missingness_qvalue

    score_lookup = {
        (
            row["trait"],
            row["scheme"],
            row["chrom"],
            row["pos"],
            row["ref"],
            row["alt"],
        ): row
        for row in score_rows
    }

    candidate_rows = []

    for trait in TRAITS:
        for variant, calls in zip(
            variants,
            all_calls,
        ):
            (
                chrom,
                position,
                ref_allele,
                alt_allele,
            ) = variant

            raw = score_lookup[
                (
                    trait,
                    "called_GT",
                    chrom,
                    position,
                    ref_allele,
                    alt_allele,
                )
            ]

            qc = score_lookup[
                (
                    trait,
                    "DP2_GQ20",
                    chrom,
                    position,
                    ref_allele,
                    alt_allele,
                )
            ]

            stable_orientation = (
                raw["selected_orientation"]
                == qc["selected_orientation"]
            )

            comparative_supported = (
                raw["score_eligible"]
                and qc["score_eligible"]
                and stable_orientation
                and raw["BH_q"] is not None
                and raw["BH_q"] <= 0.05
                and qc["BH_q"] is not None
                and qc["BH_q"] <= 0.05
            )

            left_call = calls[
                sample_index["TN_Left"]
            ]

            right_call = calls[
                sample_index["TN_Right"]
            ]

            focal_qc = (
                SCHEMES["DP2_GQ20"](
                    left_call
                )
                and SCHEMES["DP2_GQ20"](
                    right_call
                )
            )

            associated_homozygote = (
                0
                if qc[
                    "selected_orientation"
                ].startswith("REF_")
                else 2
            )

            focal_exact = (
                left_call[0] == 1
                and right_call[0]
                == associated_homozygote
            )

            if not (
                comparative_supported
                and focal_qc
                and focal_exact
            ):
                continue

            missingness_ok = (
                raw[
                    "missingness_BH_q"
                ] is not None
                and qc[
                    "missingness_BH_q"
                ] is not None
                and raw[
                    "missingness_BH_q"
                ] > 0.05
                and qc[
                    "missingness_BH_q"
                ] > 0.05
            )

            candidate_rows.append(
                {
                    "trait": trait,
                    "chrom": chrom,
                    "pos": position,
                    "ref": ref_allele,
                    "alt": alt_allele,
                    "orientation": qc[
                        "selected_orientation"
                    ],
                    "raw_BA": raw[
                        "balanced_accuracy"
                    ],
                    "QC_BA": qc[
                        "balanced_accuracy"
                    ],
                    "minimum_BA": min(
                        raw[
                            "balanced_accuracy"
                        ],
                        qc[
                            "balanced_accuracy"
                        ],
                    ),
                    "raw_MCC": raw["MCC"],
                    "QC_MCC": qc["MCC"],
                    "raw_errors": raw[
                        "misclassified_n"
                    ],
                    "QC_errors": qc[
                        "misclassified_n"
                    ],
                    "raw_q": raw["BH_q"],
                    "QC_q": qc["BH_q"],
                    "raw_missingness_q": raw[
                        "missingness_BH_q"
                    ],
                    "QC_missingness_q": qc[
                        "missingness_BH_q"
                    ],
                    "missingness_ok": (
                        missingness_ok
                    ),
                    "TN_Left_GT": (
                        canonical_gt(
                            left_call[0]
                        )
                    ),
                    "TN_Left_AD": left_call[3],
                    "TN_Left_DP": left_call[1],
                    "TN_Left_GQ": left_call[2],
                    "TN_Right_GT": (
                        canonical_gt(
                            right_call[0]
                        )
                    ),
                    "TN_Right_AD": right_call[3],
                    "TN_Right_DP": right_call[1],
                    "TN_Right_GQ": right_call[2],
                }
            )

    with SCORES.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                score_rows[0]
            ),
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(score_rows)

    with COHORTS.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                cohort_rows[0]
            ),
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(cohort_rows)

    with CANDIDATES.open(
        "w",
        newline="",
    ) as handle:
        fieldnames = (
            list(candidate_rows[0])
            if candidate_rows
            else [
                "trait",
                "chrom",
                "pos",
                "ref",
                "alt",
                "orientation",
            ]
        )

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(candidate_rows)

    print(
        "ORIGINAL 502-SAMPLE "
        "RECESSIVE Yb SCORING"
    )

    print()
    print("SCORING COHORTS")

    print(
        "trait\tVCF_samples\tcomparative_n\t"
        "phenotype_0\tphenotype_1\t"
        "focal_excluded\tpiano_keys_excluded\t"
        "missing_phenotype"
    )

    for row in cohort_rows:
        print(
            f'{row["trait"]}\t'
            f'{row["VCF_samples"]}\t'
            f'{row["comparative_scoring_samples"]}\t'
            f'{row["phenotype_0_n"]}\t'
            f'{row["phenotype_1_n"]}\t'
            f'{row["focal_excluded"]}\t'
            f'{row["piano_keys_excluded"]}\t'
            f'{row["missing_phenotype"]}'
        )

    print()

    print(
        "ELIGIBILITY AND "
        "INCLUSIVE FOCAL SITES"
    )

    print(
        "trait\tscheme\tSNPs\tcallable\t"
        "eligible\tBA_ge_0.80\tBA_ge_0.90\t"
        "inclusive_focal_sites"
    )

    for trait in TRAITS:
        inclusive_n = sum(
            row["trait"] == trait
            for row in candidate_rows
        )

        for scheme in SCHEMES:
            subset = [
                row
                for row in score_rows
                if (
                    row["trait"] == trait
                    and row["scheme"] == scheme
                )
            ]

            print(
                f"{trait}\t"
                f"{scheme}\t"
                f"{len(subset)}\t"
                f"{sum(boolean(row['callability_eligible']) for row in subset)}\t"
                f"{sum(boolean(row['score_eligible']) for row in subset)}\t"
                f"{sum(boolean(row['score_eligible']) and row['balanced_accuracy'] >= 0.80 for row in subset)}\t"
                f"{sum(boolean(row['score_eligible']) and row['balanced_accuracy'] >= 0.90 for row in subset)}\t"
                f"{inclusive_n}"
            )

    print()
    print(
        "INCLUSIVE FOCAL-COMPATIBLE SITES"
    )

    print(
        "trait\tposition\torientation\t"
        "min_BA\tleft_GT\tright_GT\t"
        "missingness_ok"
    )

    for row in sorted(
        candidate_rows,
        key=lambda item: (
            item["trait"],
            -item["minimum_BA"],
            item["pos"],
        ),
    ):
        print(
            f'{row["trait"]}\t'
            f'{row["chrom"]}:{row["pos"]}\t'
            f'{row["orientation"]}\t'
            f'{row["minimum_BA"]:.4f}\t'
            f'{row["TN_Left_GT"]}\t'
            f'{row["TN_Right_GT"]}\t'
            f'{row["missingness_ok"]}'
        )

    print()
    print(f"Scores: {SCORES}")
    print(f"Candidates: {CANDIDATES}")
    print(f"Cohorts: {COHORTS}")

    print(
        "Missingness was recorded but was not "
        "used to exclude inclusive candidates."
    )


if __name__ == "__main__":
    main()
