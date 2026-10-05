#!/usr/bin/env python3

from pathlib import Path
import csv
import math
import re

from scipy.stats import fisher_exact


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

SOURCE_DIR = (
    PROJECT
    / "external/krista_genotype_plots_original/"
    "Genotype_Plots"
)

METADATA_FILE = (
    PROJECT
    / "data/derived/region_panel_metadata.tsv"
)

OUTPUT_DIR = (
    PROJECT
    / "analysis/original_502/checkpoints"
)

SCORE_OUTPUT = (
    OUTPUT_DIR
    / "optix_original_502_comparative_scores.tsv"
)

COHORT_OUTPUT = (
    OUTPUT_DIR
    / "optix_original_502_scoring_cohorts.tsv"
)

COHORT_MANIFEST_OUTPUT = (
    OUTPUT_DIR
    / "optix_original_502_cohort_manifest.tsv"
)

PANELS = {
    "band": {
        "trait": "B",
        "expected_orientation": (
            "REF_carrier_predicts_presence"
        ),
        "vcf": (
            SOURCE_DIR
            / "optix_band-1112651-1114314.vcf"
        ),
    },
    "dennis1": {
        "trait": "D",
        "expected_orientation": (
            "ALT_carrier_predicts_presence"
        ),
        "vcf": (
            SOURCE_DIR
            / "dennis1-1215210-1216124.vcf"
        ),
    },
    "dennis2": {
        "trait": "D",
        "expected_orientation": (
            "ALT_carrier_predicts_presence"
        ),
        "vcf": (
            SOURCE_DIR
            / "dennis2-1218557-1221494.vcf"
        ),
    },
    "ray1": {
        "trait": "R",
        "expected_orientation": (
            "ALT_carrier_predicts_presence"
        ),
        "vcf": (
            SOURCE_DIR
            / "ray1-1156746-1160220.vcf"
        ),
    },
    "ray2": {
        "trait": "R",
        "expected_orientation": (
            "ALT_carrier_predicts_presence"
        ),
        "vcf": (
            SOURCE_DIR
            / "ray2-1174547-1187833.vcf"
        ),
    },
}

SCHEMES = {
    "called_GT": {
        "minimum_DP": None,
        "minimum_GQ": None,
    },
    "DP2_GQ20": {
        "minimum_DP": 2,
        "minimum_GQ": 20,
    },
}

FOCAL_SAMPLES = {
    "TN2L_S1",
    "TN3L_S3",
    "TN2R_S2",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
}

POOLED_SAMPLES = [
    "TN_Left",
    "TN_Right",
]

MINIMUM_OVERALL_CALL_RATE = 0.70
MINIMUM_STATE_CALL_RATE = 0.50
MINIMUM_CALLED_PER_PHENOTYPE = 20
MINIMUM_PREDICTED_CLASS = 10


def clean_text(value):
    return " ".join(str(value or "").split())


def parse_phenotype(value):
    value = clean_text(value)

    if value in {"0", "0.0"}:
        return 0

    if value in {"1", "1.0"}:
        return 1

    return None


def canonical_gt(value):
    if value in {None, "", ".", "./.", ".|."}:
        return None

    alleles = re.split(r"[/|]", value)

    if any(allele in {"", "."} for allele in alleles):
        return None

    try:
        numeric = sorted(int(allele) for allele in alleles)
    except ValueError:
        return None

    return "/".join(str(allele) for allele in numeric)


def parse_number(value):
    if value in {None, "", "."}:
        return None

    try:
        return float(value)
    except ValueError:
        return None


def parse_integer_vector(value):
    if value in {None, "", "."}:
        return None

    try:
        return [int(part) for part in value.split(",")]
    except ValueError:
        return None


def parse_call(format_ids, sample_value):
    values = sample_value.split(":")
    fields = dict(zip(format_ids, values))

    gt = canonical_gt(fields.get("GT"))
    dp = parse_number(fields.get("DP"))
    gq = parse_number(fields.get("GQ"))
    ad = parse_integer_vector(fields.get("AD"))

    alt_fraction = None

    if ad is not None and len(ad) == 2:
        total = ad[0] + ad[1]

        if total > 0:
            alt_fraction = ad[1] / total

    return {
        "GT": gt,
        "DP": dp,
        "GQ": gq,
        "AD": ad,
        "ALT_fraction": alt_fraction,
    }


def call_passes_scheme(call, scheme):
    if call["GT"] is None:
        return False

    minimum_dp = scheme["minimum_DP"]
    minimum_gq = scheme["minimum_GQ"]

    if minimum_dp is not None:
        if call["DP"] is None or call["DP"] < minimum_dp:
            return False

    if minimum_gq is not None:
        if call["GQ"] is None or call["GQ"] < minimum_gq:
            return False

    return True


def carrier_prediction(gt, orientation):
    if gt is None:
        return None

    alleles = [int(value) for value in gt.split("/")]

    if orientation == "ALT_carrier_predicts_presence":
        return 1 if any(value > 0 for value in alleles) else 0

    if orientation == "REF_carrier_predicts_presence":
        return 1 if any(value == 0 for value in alleles) else 0

    return None


def calculate_metrics(phenotypes, predictions):
    tp = sum(
        phenotype == 1 and prediction == 1
        for phenotype, prediction
        in zip(phenotypes, predictions)
    )

    tn = sum(
        phenotype == 0 and prediction == 0
        for phenotype, prediction
        in zip(phenotypes, predictions)
    )

    fp = sum(
        phenotype == 0 and prediction == 1
        for phenotype, prediction
        in zip(phenotypes, predictions)
    )

    fn = sum(
        phenotype == 1 and prediction == 0
        for phenotype, prediction
        in zip(phenotypes, predictions)
    )

    sensitivity = (
        tp / (tp + fn)
        if (tp + fn) > 0
        else None
    )

    specificity = (
        tn / (tn + fp)
        if (tn + fp) > 0
        else None
    )

    balanced_accuracy = (
        (sensitivity + specificity) / 2
        if sensitivity is not None
        and specificity is not None
        else None
    )

    denominator = math.sqrt(
        (tp + fp)
        * (tp + fn)
        * (tn + fp)
        * (tn + fn)
    )

    mcc = (
        ((tp * tn) - (fp * fn)) / denominator
        if denominator > 0
        else None
    )

    fisher_p = fisher_exact(
        [
            [tp, fn],
            [fp, tn],
        ],
        alternative="two-sided",
    ).pvalue

    return {
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "predicted_0": tn + fn,
        "predicted_1": tp + fp,
        "misclassified_n": fp + fn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "balanced_accuracy": balanced_accuracy,
        "MCC": mcc,
        "fisher_p": fisher_p,
    }


def select_orientation(
    phenotypes,
    genotypes,
    expected_orientation,
):
    results = []

    for orientation in [
        "ALT_carrier_predicts_presence",
        "REF_carrier_predicts_presence",
    ]:
        predictions = [
            carrier_prediction(gt, orientation)
            for gt in genotypes
        ]

        metrics = calculate_metrics(
            phenotypes,
            predictions,
        )

        metrics["orientation"] = orientation
        metrics["predictions"] = predictions

        results.append(metrics)

    def rank(result):
        ba = (
            result["balanced_accuracy"]
            if result["balanced_accuracy"] is not None
            else -1
        )

        mcc = (
            result["MCC"]
            if result["MCC"] is not None
            else -2
        )

        expected_tie_break = (
            1
            if result["orientation"]
            == expected_orientation
            else 0
        )

        return (
            ba,
            mcc,
            -result["misclassified_n"],
            expected_tie_break,
        )

    return max(results, key=rank)


def bh_adjust(p_values):
    indexed = [
        (index, value)
        for index, value in enumerate(p_values)
        if value is not None
        and math.isfinite(value)
    ]

    adjusted = [None] * len(p_values)

    if not indexed:
        return adjusted

    ordered = sorted(
        indexed,
        key=lambda item: item[1],
    )

    total = len(ordered)
    running_minimum = 1.0

    for reverse_rank in range(total - 1, -1, -1):
        original_index, p_value = ordered[reverse_rank]
        rank = reverse_rank + 1

        candidate = min(
            1.0,
            p_value * total / rank,
        )

        running_minimum = min(
            running_minimum,
            candidate,
        )

        adjusted[original_index] = running_minimum

    return adjusted


def vector_text(value):
    if value is None:
        return ""
    return ",".join(str(part) for part in value)


def load_metadata():
    metadata = {}

    with METADATA_FILE.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            accession = clean_text(
                row.get("Accession")
            )

            if not accession:
                continue

            taxon = clean_text(
                row.get("Taxon Name")
            )

            metadata[accession] = {
                "taxon": taxon,
                "data_type": clean_text(
                    row.get("Data Type")
                ),
                "B": parse_phenotype(
                    row.get("B")
                ),
                "R": parse_phenotype(
                    row.get("R")
                ),
                "D": parse_phenotype(
                    row.get("D")
                ),
                "piano_key": (
                    "piano keys" in taxon.lower()
                ),
                "focal": accession in FOCAL_SAMPLES,
            }

    return metadata


def make_cohort(samples, metadata, trait):
    included = []
    manifest = []

    for sample in samples:
        entry = metadata.get(sample)

        if entry is None:
            reason = "no_metadata_match"
            phenotype = None
        elif entry["focal"]:
            reason = "focal_excluded"
            phenotype = entry[trait]
        elif entry["piano_key"]:
            reason = "piano_key_excluded"
            phenotype = entry[trait]
        elif entry[trait] is None:
            reason = "missing_phenotype"
            phenotype = None
        else:
            reason = "included"
            phenotype = entry[trait]

            included.append(
                {
                    "sample": sample,
                    "phenotype": phenotype,
                    "data_type": entry["data_type"],
                    "taxon": entry["taxon"],
                }
            )

        manifest.append(
            {
                "sample": sample,
                "trait": trait,
                "phenotype": (
                    ""
                    if phenotype is None
                    else phenotype
                ),
                "reason": reason,
                "taxon": (
                    ""
                    if entry is None
                    else entry["taxon"]
                ),
                "data_type": (
                    ""
                    if entry is None
                    else entry["data_type"]
                ),
            }
        )

    return included, manifest


def score_panel(panel, specification, metadata):
    path = specification["vcf"]
    trait = specification["trait"]
    expected_orientation = specification[
        "expected_orientation"
    ]

    rows = []
    manifest_rows = []
    cohort_summary = None

    with path.open(errors="replace") as handle:
        samples = []
        sample_indices = {}
        cohort = []
        cohort_indices = []

        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]

                sample_indices = {
                    sample: index
                    for index, sample in enumerate(samples)
                }

                cohort, manifest = make_cohort(
                    samples,
                    metadata,
                    trait,
                )

                for manifest_row in manifest:
                    manifest_rows.append(
                        {
                            "panel": panel,
                            **manifest_row,
                        }
                    )

                cohort_indices = [
                    sample_indices[row["sample"]]
                    for row in cohort
                ]

                phenotype_0_n = sum(
                    row["phenotype"] == 0
                    for row in cohort
                )

                phenotype_1_n = sum(
                    row["phenotype"] == 1
                    for row in cohort
                )

                cohort_summary = {
                    "panel": panel,
                    "trait": trait,
                    "VCF_samples": len(samples),
                    "comparative_scoring_samples": (
                        len(cohort)
                    ),
                    "phenotype_0_n": phenotype_0_n,
                    "phenotype_1_n": phenotype_1_n,
                    "focal_excluded": sum(
                        row["reason"] == "focal_excluded"
                        for row in manifest
                    ),
                    "piano_keys_excluded": sum(
                        row["reason"]
                        == "piano_key_excluded"
                        for row in manifest
                    ),
                    "missing_phenotype": sum(
                        row["reason"]
                        == "missing_phenotype"
                        for row in manifest
                    ),
                    "no_metadata_match": sum(
                        row["reason"]
                        == "no_metadata_match"
                        for row in manifest
                    ),
                }

                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            if not (
                "," not in alt
                and len(ref) == 1
                and len(alt) == 1
            ):
                continue

            format_ids = fields[8].split(":")
            sample_values = fields[9:]

            cohort_calls = [
                parse_call(
                    format_ids,
                    sample_values[index],
                )
                for index in cohort_indices
            ]

            pooled_calls = {}

            for pooled_sample in POOLED_SAMPLES:
                index = sample_indices.get(pooled_sample)

                if index is None:
                    pooled_calls[pooled_sample] = {
                        "GT": None,
                        "DP": None,
                        "GQ": None,
                        "AD": None,
                        "ALT_fraction": None,
                    }
                else:
                    pooled_calls[pooled_sample] = (
                        parse_call(
                            format_ids,
                            sample_values[index],
                        )
                    )

            for scheme_name, scheme in SCHEMES.items():
                called_phenotypes = []
                called_genotypes = []

                phenotype_called = {
                    0: 0,
                    1: 0,
                }

                phenotype_missing = {
                    0: 0,
                    1: 0,
                }

                for cohort_row, call in zip(
                    cohort,
                    cohort_calls,
                ):
                    phenotype = cohort_row["phenotype"]

                    if call_passes_scheme(call, scheme):
                        called_phenotypes.append(
                            phenotype
                        )

                        called_genotypes.append(
                            call["GT"]
                        )

                        phenotype_called[
                            phenotype
                        ] += 1
                    else:
                        phenotype_missing[
                            phenotype
                        ] += 1

                called_n = len(called_genotypes)
                cohort_n = len(cohort)

                call_rate = (
                    called_n / cohort_n
                    if cohort_n > 0
                    else None
                )

                phenotype_0_total = (
                    phenotype_called[0]
                    + phenotype_missing[0]
                )

                phenotype_1_total = (
                    phenotype_called[1]
                    + phenotype_missing[1]
                )

                phenotype_0_call_rate = (
                    phenotype_called[0]
                    / phenotype_0_total
                    if phenotype_0_total > 0
                    else None
                )

                phenotype_1_call_rate = (
                    phenotype_called[1]
                    / phenotype_1_total
                    if phenotype_1_total > 0
                    else None
                )

                missingness_p = fisher_exact(
                    [
                        [
                            phenotype_called[0],
                            phenotype_missing[0],
                        ],
                        [
                            phenotype_called[1],
                            phenotype_missing[1],
                        ],
                    ],
                    alternative="two-sided",
                ).pvalue

                selected = select_orientation(
                    called_phenotypes,
                    called_genotypes,
                    expected_orientation,
                )

                callability_eligible = (
                    call_rate is not None
                    and call_rate
                    >= MINIMUM_OVERALL_CALL_RATE
                    and phenotype_0_call_rate is not None
                    and phenotype_0_call_rate
                    >= MINIMUM_STATE_CALL_RATE
                    and phenotype_1_call_rate is not None
                    and phenotype_1_call_rate
                    >= MINIMUM_STATE_CALL_RATE
                    and phenotype_called[0]
                    >= MINIMUM_CALLED_PER_PHENOTYPE
                    and phenotype_called[1]
                    >= MINIMUM_CALLED_PER_PHENOTYPE
                )

                score_eligible = (
                    callability_eligible
                    and selected["MCC"] is not None
                    and selected["predicted_0"]
                    >= MINIMUM_PREDICTED_CLASS
                    and selected["predicted_1"]
                    >= MINIMUM_PREDICTED_CLASS
                )

                left = pooled_calls["TN_Left"]
                right = pooled_calls["TN_Right"]

                left_prediction = carrier_prediction(
                    left["GT"],
                    selected["orientation"],
                )

                right_prediction = carrier_prediction(
                    right["GT"],
                    selected["orientation"],
                )

                expected_left = (
                    0 if trait == "B" else 1
                )

                expected_right = 1

                dominant_focal_match = (
                    left_prediction == expected_left
                    and right_prediction == expected_right
                )

                pooled_both_qc = (
                    call_passes_scheme(
                        left,
                        SCHEMES["DP2_GQ20"],
                    )
                    and call_passes_scheme(
                        right,
                        SCHEMES["DP2_GQ20"],
                    )
                )

                rows.append(
                    {
                        "panel": panel,
                        "trait": trait,
                        "scheme": scheme_name,
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "comparative_cohort_n": cohort_n,
                        "phenotype_0_total": (
                            phenotype_0_total
                        ),
                        "phenotype_1_total": (
                            phenotype_1_total
                        ),
                        "phenotype_0_called": (
                            phenotype_called[0]
                        ),
                        "phenotype_1_called": (
                            phenotype_called[1]
                        ),
                        "call_rate": call_rate,
                        "phenotype_0_call_rate": (
                            phenotype_0_call_rate
                        ),
                        "phenotype_1_call_rate": (
                            phenotype_1_call_rate
                        ),
                        "missingness_fisher_p": (
                            missingness_p
                        ),
                        "selected_orientation": (
                            selected["orientation"]
                        ),
                        "reference_expectation": (
                            expected_orientation
                        ),
                        "matches_reference_expectation": (
                            selected["orientation"]
                            == expected_orientation
                        ),
                        "predicted_0": (
                            selected["predicted_0"]
                        ),
                        "predicted_1": (
                            selected["predicted_1"]
                        ),
                        "tp": selected["tp"],
                        "tn": selected["tn"],
                        "fp": selected["fp"],
                        "fn": selected["fn"],
                        "misclassified_n": (
                            selected["misclassified_n"]
                        ),
                        "sensitivity": (
                            selected["sensitivity"]
                        ),
                        "specificity": (
                            selected["specificity"]
                        ),
                        "balanced_accuracy": (
                            selected["balanced_accuracy"]
                        ),
                        "MCC": selected["MCC"],
                        "fisher_p": selected["fisher_p"],
                        "BH_q": None,
                        "missingness_BH_q": None,
                        "callability_eligible": (
                            callability_eligible
                        ),
                        "score_eligible": score_eligible,
                        "TN_Left_GT": left["GT"] or "",
                        "TN_Left_AD": vector_text(
                            left["AD"]
                        ),
                        "TN_Left_DP": left["DP"],
                        "TN_Left_GQ": left["GQ"],
                        "TN_Left_ALT_fraction": (
                            left["ALT_fraction"]
                        ),
                        "TN_Right_GT": right["GT"] or "",
                        "TN_Right_AD": vector_text(
                            right["AD"]
                        ),
                        "TN_Right_DP": right["DP"],
                        "TN_Right_GQ": right["GQ"],
                        "TN_Right_ALT_fraction": (
                            right["ALT_fraction"]
                        ),
                        "pooled_both_DP2_GQ20": (
                            pooled_both_qc
                        ),
                        "pooled_GT_different": (
                            left["GT"] is not None
                            and right["GT"] is not None
                            and left["GT"] != right["GT"]
                        ),
                        "left_11_right_01": (
                            left["GT"] == "1/1"
                            and right["GT"] == "0/1"
                        ),
                        "dominant_focal_model_match": (
                            dominant_focal_match
                        ),
                    }
                )

    for scheme_name in SCHEMES:
        indices = [
            index for index, row in enumerate(rows)
            if row["scheme"] == scheme_name
        ]

        fisher_q = bh_adjust(
            [rows[index]["fisher_p"] for index in indices]
        )

        missingness_q = bh_adjust(
            [
                rows[index]["missingness_fisher_p"]
                for index in indices
            ]
        )

        for local_index, row_index in enumerate(indices):
            rows[row_index]["BH_q"] = fisher_q[local_index]
            rows[row_index][
                "missingness_BH_q"
            ] = missingness_q[local_index]

    return rows, cohort_summary, manifest_rows


def main():
    required = [
        METADATA_FILE,
        *[
            specification["vcf"]
            for specification in PANELS.values()
        ],
    ]

    missing = [path for path in required if not path.is_file()]

    if missing:
        print("ERROR: required inputs are missing:")
        for path in missing:
            print(path)
        return

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    metadata = load_metadata()

    all_scores = []
    cohort_rows = []
    manifest_rows = []

    for panel, specification in PANELS.items():
        scores, cohort, manifest = score_panel(
            panel,
            specification,
            metadata,
        )

        all_scores.extend(scores)
        cohort_rows.append(cohort)
        manifest_rows.extend(manifest)

    with SCORE_OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(all_scores[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(all_scores)

    with COHORT_OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(cohort_rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(cohort_rows)

    with COHORT_MANIFEST_OUTPUT.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(manifest_rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    print("ORIGINAL 502-SAMPLE OPTIX SCORING")
    print()
    print("SCORING COHORTS")
    print(
        "panel\ttrait\tVCF_samples\tcomparative_n\t"
        "phenotype_0\tphenotype_1\tfocal_excluded\t"
        "piano_keys_excluded\tmissing_phenotype\t"
        "no_metadata"
    )

    for row in cohort_rows:
        print(
            f'{row["panel"]}\t'
            f'{row["trait"]}\t'
            f'{row["VCF_samples"]}\t'
            f'{row["comparative_scoring_samples"]}\t'
            f'{row["phenotype_0_n"]}\t'
            f'{row["phenotype_1_n"]}\t'
            f'{row["focal_excluded"]}\t'
            f'{row["piano_keys_excluded"]}\t'
            f'{row["missing_phenotype"]}\t'
            f'{row["no_metadata_match"]}'
        )

    print()
    print("SITE ELIGIBILITY")
    print(
        "panel\tscheme\tSNPs\tcallable\teligible\t"
        "BA_ge_0.80\tBA_ge_0.90\tq_le_0.05\t"
        "missingness_q_le_0.05"
    )

    for panel in PANELS:
        for scheme_name in SCHEMES:
            rows = [
                row for row in all_scores
                if row["panel"] == panel
                and row["scheme"] == scheme_name
            ]

            print(
                f"{panel}\t"
                f"{scheme_name}\t"
                f"{len(rows)}\t"
                f"{sum(row['callability_eligible'] for row in rows)}\t"
                f"{sum(row['score_eligible'] for row in rows)}\t"
                f"{sum(row['score_eligible'] and row['balanced_accuracy'] is not None and row['balanced_accuracy'] >= 0.80 for row in rows)}\t"
                f"{sum(row['score_eligible'] and row['balanced_accuracy'] is not None and row['balanced_accuracy'] >= 0.90 for row in rows)}\t"
                f"{sum(row['BH_q'] is not None and row['BH_q'] <= 0.05 for row in rows)}\t"
                f"{sum(row['missingness_BH_q'] is not None and row['missingness_BH_q'] <= 0.05 for row in rows)}"
            )

    print()
    print("TOP DP2/GQ20 SITES")
    print(
        "panel\tposition\torientation\tBA\tMCC\terrors\t"
        "q\tmissingness_q\tleft_GT\tright_GT\tL11_R01"
    )

    for panel in PANELS:
        rows = [
            row for row in all_scores
            if row["panel"] == panel
            and row["scheme"] == "DP2_GQ20"
            and row["score_eligible"]
        ]

        rows.sort(
            key=lambda row: (
                -(
                    row["balanced_accuracy"]
                    if row["balanced_accuracy"] is not None
                    else -1
                ),
                -(
                    row["MCC"]
                    if row["MCC"] is not None
                    else -2
                ),
                row["misclassified_n"],
                row["pos"],
            )
        )

        for row in rows[:10]:
            print(
                f'{panel}\t'
                f'{row["chrom"]}:{row["pos"]}\t'
                f'{row["selected_orientation"]}\t'
                f'{row["balanced_accuracy"]:.4f}\t'
                f'{row["MCC"]:.4f}\t'
                f'{row["misclassified_n"]}\t'
                f'{row["BH_q"]:.4g}\t'
                f'{row["missingness_BH_q"]:.4g}\t'
                f'{row["TN_Left_GT"]}\t'
                f'{row["TN_Right_GT"]}\t'
                f'{row["left_11_right_01"]}'
            )

    print()
    print(f"Scores: {SCORE_OUTPUT}")
    print(f"Cohorts: {COHORT_OUTPUT}")
    print(f"Cohort manifest: {COHORT_MANIFEST_OUTPUT}")
    print()
    print(
        "Focal and piano-key samples were not used in "
        "comparative scoring."
    )


if __name__ == "__main__":
    main()
