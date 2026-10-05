#!/usr/bin/env python3

import csv
import math
import re
import subprocess
from collections import defaultdict
from pathlib import Path

from openpyxl import load_workbook

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

VCF_DIR = (
    PROJECT
    / "data"
    / "derived"
    / "reconstructed_plot_inputs"
)

WORKBOOK = (
    PROJECT
    / "metadata"
    / "Wing_Pattern_MetaData.xlsx"
)

CHECKPOINT_DIR = PROJECT / "analysis" / "checkpoints"

SCORE_OUTPUT = (
    CHECKPOINT_DIR
    / "optix_51sample_comparative_scores.tsv"
)

ELIGIBILITY_OUTPUT = (
    CHECKPOINT_DIR
    / "optix_51sample_score_eligibility.tsv"
)

COHORT_OUTPUT = (
    CHECKPOINT_DIR
    / "optix_51sample_scoring_cohorts.tsv"
)

PANELS = {
    "band": {
        "vcf": VCF_DIR / "optix_band-1112651-1114314.vcf",
        "trait": "B",
        "reference_expectation": "REF_carrier_predicts_presence",
    },
    "dennis1": {
        "vcf": VCF_DIR / "dennis1-1215210-1216124.vcf",
        "trait": "D",
        "reference_expectation": "ALT_carrier_predicts_presence",
    },
    "dennis2": {
        "vcf": VCF_DIR / "dennis2-1218557-1221494.vcf",
        "trait": "D",
        "reference_expectation": "ALT_carrier_predicts_presence",
    },
    "ray1": {
        "vcf": VCF_DIR / "ray1-1156746-1160220.vcf",
        "trait": "R",
        "reference_expectation": "ALT_carrier_predicts_presence",
    },
    "ray2": {
        "vcf": VCF_DIR / "ray2-1174547-1187833.vcf",
        "trait": "R",
        "reference_expectation": "ALT_carrier_predicts_presence",
    },
}

SCHEMES = {
    "called_GT": {
        "minimum_dp": None,
        "minimum_gq": None,
    },
    "DP2_GQ20": {
        "minimum_dp": 2,
        "minimum_gq": 20,
    },
}

FOCAL_SAMPLES = {
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
}


def normalized_header(value):
    if value is None:
        return ""

    return re.sub(r"\s+", " ", str(value)).strip()


def parse_binary(value):
    if value is None:
        return None

    text = str(value).strip()

    if text in {"0", "0.0"}:
        return 0

    if text in {"1", "1.0"}:
        return 1

    return None


def load_metadata():
    workbook = load_workbook(
        WORKBOOK,
        read_only=True,
        data_only=True,
    )

    worksheet = workbook["Metadata"]
    rows = worksheet.iter_rows(values_only=True)

    headers = [
        normalized_header(value)
        for value in next(rows)
    ]

    header_index = {
        header: index
        for index, header in enumerate(headers)
    }

    required = ["Accession", "B", "R", "D"]

    missing = [
        header
        for header in required
        if header not in header_index
    ]

    if missing:
        return {}, (
            "Missing workbook columns: "
            + ",".join(missing)
        )

    metadata = {}

    for values in rows:
        accession_value = values[
            header_index["Accession"]
        ]

        if accession_value is None:
            continue

        accession = str(accession_value).strip()

        if not accession:
            continue

        metadata[accession] = {
            "B": parse_binary(
                values[header_index["B"]]
            ),
            "R": parse_binary(
                values[header_index["R"]]
            ),
            "D": parse_binary(
                values[header_index["D"]]
            ),
        }

    return metadata, ""


def run_capture(command):
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    return result.returncode, result.stdout, result.stderr


def get_samples(vcf):
    code, stdout, stderr = run_capture(
        ["bcftools", "query", "-l", str(vcf)]
    )

    if code != 0:
        return [], stderr.strip()

    return [
        line.strip()
        for line in stdout.splitlines()
        if line.strip()
    ], ""


def query_records(vcf):
    code, stdout, stderr = run_capture(
        [
            "bcftools",
            "query",
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT"
            "[\t%GT\t%DP\t%GQ]\n",
            str(vcf),
        ]
    )

    if code != 0:
        return [], stderr.strip()

    records = []

    for line in stdout.splitlines():
        fields = line.split("\t")

        if len(fields) < 4:
            continue

        sample_fields = fields[4:]

        if len(sample_fields) % 3 != 0:
            continue

        calls = []

        for index in range(
            0,
            len(sample_fields),
            3,
        ):
            calls.append(
                {
                    "GT": sample_fields[index],
                    "DP": sample_fields[index + 1],
                    "GQ": sample_fields[index + 2],
                }
            )

        records.append(
            {
                "chrom": fields[0],
                "pos": int(fields[1]),
                "ref": fields[2],
                "alt": fields[3],
                "calls": calls,
            }
        )

    return records, ""


def parse_number(value):
    if value in {"", ".", None}:
        return None

    try:
        return float(value)
    except ValueError:
        return None


def parse_gt(value):
    if not value or "." in value:
        return None

    alleles = re.split(r"[/|]", value)

    try:
        return tuple(
            sorted(int(allele) for allele in alleles)
        )
    except ValueError:
        return None


def accepted_gt(call, scheme):
    gt = parse_gt(call["GT"])

    if gt is None:
        return None

    minimum_dp = scheme["minimum_dp"]
    minimum_gq = scheme["minimum_gq"]

    if minimum_dp is not None:
        dp = parse_number(call["DP"])

        if dp is None or dp < minimum_dp:
            return None

    if minimum_gq is not None:
        gq = parse_number(call["GQ"])

        if gq is None or gq < minimum_gq:
            return None

    return gt


def safe_fraction(numerator, denominator):
    if denominator == 0:
        return None

    return numerator / denominator


def mcc(tp, tn, fp, fn):
    denominator = math.sqrt(
        (tp + fp)
        * (tp + fn)
        * (tn + fp)
        * (tn + fn)
    )

    if denominator == 0:
        return 0.0

    return (
        (tp * tn) - (fp * fn)
    ) / denominator


def hypergeometric_probability(
    x,
    row_1,
    row_2,
    column_1,
    total,
):
    return (
        math.comb(row_1, x)
        * math.comb(row_2, column_1 - x)
        / math.comb(total, column_1)
    )


def fisher_two_sided(a, b, c, d):
    row_1 = a + b
    row_2 = c + d
    column_1 = a + c
    total = row_1 + row_2

    if total == 0:
        return None

    minimum_x = max(0, column_1 - row_2)
    maximum_x = min(row_1, column_1)

    observed = hypergeometric_probability(
        a,
        row_1,
        row_2,
        column_1,
        total,
    )

    probability = 0.0

    for x in range(minimum_x, maximum_x + 1):
        candidate = hypergeometric_probability(
            x,
            row_1,
            row_2,
            column_1,
            total,
        )

        if candidate <= observed + 1e-12:
            probability += candidate

    return min(probability, 1.0)


def confusion_metrics(
    phenotypes,
    genotypes,
    orientation,
):
    tp = tn = fp = fn = 0

    for phenotype, genotype in zip(
        phenotypes,
        genotypes,
    ):
        if orientation == "ALT_carrier_predicts_presence":
            predicted = int(any(
                allele > 0 for allele in genotype
            ))
        else:
            predicted = int(any(
                allele == 0 for allele in genotype
            ))

        if phenotype == 1 and predicted == 1:
            tp += 1
        elif phenotype == 1 and predicted == 0:
            fn += 1
        elif phenotype == 0 and predicted == 1:
            fp += 1
        else:
            tn += 1

    sensitivity = safe_fraction(tp, tp + fn)
    specificity = safe_fraction(tn, tn + fp)

    if sensitivity is None or specificity is None:
        balanced_accuracy = None
    else:
        balanced_accuracy = (
            sensitivity + specificity
        ) / 2

    return {
        "orientation": orientation,
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "predicted_1": tp + fp,
        "predicted_0": tn + fn,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "balanced_accuracy": balanced_accuracy,
        "MCC": mcc(tp, tn, fp, fn),
        "fisher_p": fisher_two_sided(
            tp,
            fn,
            fp,
            tn,
        ),
        "misclassified_n": fp + fn,
    }


def choose_orientation(
    phenotypes,
    genotypes,
):
    candidates = [
        confusion_metrics(
            phenotypes,
            genotypes,
            "ALT_carrier_predicts_presence",
        ),
        confusion_metrics(
            phenotypes,
            genotypes,
            "REF_carrier_predicts_presence",
        ),
    ]

    return max(
        candidates,
        key=lambda result: (
            result["balanced_accuracy"]
            if result["balanced_accuracy"] is not None
            else -1,
            result["MCC"],
            -result["misclassified_n"],
        ),
    )


def bh_adjust(rows, p_column, q_column):
    grouped = defaultdict(list)

    for index, row in enumerate(rows):
        value = row.get(p_column)

        if value is None:
            continue

        key = (
            row["panel"],
            row["trait"],
            row["scheme"],
        )

        grouped[key].append(
            (index, float(value))
        )

    for group in grouped.values():
        ordered = sorted(
            group,
            key=lambda item: item[1],
        )

        total = len(ordered)
        adjusted = [None] * total
        running = 1.0

        for reverse_index in range(
            total - 1,
            -1,
            -1,
        ):
            original_index, probability = ordered[
                reverse_index
            ]

            rank = reverse_index + 1
            candidate = min(
                probability * total / rank,
                1.0,
            )
            running = min(running, candidate)
            adjusted[reverse_index] = running

        for (
            original_index,
            probability,
        ), adjusted_probability in zip(
            ordered,
            adjusted,
        ):
            rows[original_index][q_column] = (
                adjusted_probability
            )


def main():
    CHECKPOINT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    missing_inputs = [
        str(panel["vcf"])
        for panel in PANELS.values()
        if not panel["vcf"].is_file()
    ]

    if not WORKBOOK.is_file():
        missing_inputs.append(str(WORKBOOK))

    if missing_inputs:
        print("ERROR: required inputs are missing:")
        for path in missing_inputs:
            print(path)
        return

    metadata, metadata_error = load_metadata()

    if metadata_error:
        print(f"ERROR: {metadata_error}")
        return

    score_rows = []
    cohort_rows = []

    canonical_samples = None

    for panel_name, panel in PANELS.items():
        vcf = panel["vcf"]
        trait = panel["trait"]

        samples, sample_error = get_samples(vcf)

        if sample_error:
            print(
                f"ERROR reading {panel_name}: "
                f"{sample_error}"
            )
            return

        if canonical_samples is None:
            canonical_samples = samples
        elif samples != canonical_samples:
            print(
                "ERROR: reconstructed Optix VCF "
                f"sample order differs in {panel_name}."
            )
            return

        records, record_error = query_records(vcf)

        if record_error:
            print(
                f"ERROR querying {panel_name}: "
                f"{record_error}"
            )
            return

        comparative_indices = []
        phenotype_by_index = {}

        for index, sample in enumerate(samples):
            if sample in FOCAL_SAMPLES:
                continue

            phenotype = metadata.get(
                sample,
                {},
            ).get(trait)

            if phenotype is None:
                continue

            comparative_indices.append(index)
            phenotype_by_index[index] = phenotype

        phenotype_0_n = sum(
            phenotype_by_index[index] == 0
            for index in comparative_indices
        )
        phenotype_1_n = sum(
            phenotype_by_index[index] == 1
            for index in comparative_indices
        )

        cohort_rows.append(
            {
                "panel": panel_name,
                "trait": trait,
                "VCF_samples": len(samples),
                "comparative_scoring_samples": len(
                    comparative_indices
                ),
                "phenotype_0_n": phenotype_0_n,
                "phenotype_1_n": phenotype_1_n,
                "focal_samples_excluded": sum(
                    sample in FOCAL_SAMPLES
                    for sample in samples
                ),
                "unmatched_or_missing_phenotype": (
                    len(samples)
                    - len(comparative_indices)
                    - sum(
                        sample in FOCAL_SAMPLES
                        for sample in samples
                    )
                ),
            }
        )

        for record in records:
            if (
                len(record["ref"]) != 1
                or len(record["alt"]) != 1
                or "," in record["alt"]
            ):
                continue

            for scheme_name, scheme in SCHEMES.items():
                phenotypes = []
                genotypes = []

                called_0 = called_1 = 0
                missing_0 = missing_1 = 0

                for index in comparative_indices:
                    phenotype = phenotype_by_index[
                        index
                    ]

                    genotype = accepted_gt(
                        record["calls"][index],
                        scheme,
                    )

                    if genotype is None:
                        if phenotype == 0:
                            missing_0 += 1
                        else:
                            missing_1 += 1
                        continue

                    phenotypes.append(phenotype)
                    genotypes.append(genotype)

                    if phenotype == 0:
                        called_0 += 1
                    else:
                        called_1 += 1

                total_0 = called_0 + missing_0
                total_1 = called_1 + missing_1
                called_total = called_0 + called_1
                cohort_total = total_0 + total_1

                call_rate = safe_fraction(
                    called_total,
                    cohort_total,
                )
                call_rate_0 = safe_fraction(
                    called_0,
                    total_0,
                )
                call_rate_1 = safe_fraction(
                    called_1,
                    total_1,
                )

                metrics = choose_orientation(
                    phenotypes,
                    genotypes,
                )

                callability_eligible = (
                    called_0 >= 5
                    and called_1 >= 5
                    and call_rate is not None
                    and call_rate >= 0.70
                    and call_rate_0 is not None
                    and call_rate_0 >= 0.50
                    and call_rate_1 is not None
                    and call_rate_1 >= 0.50
                )

                score_eligible = (
                    callability_eligible
                    and metrics["predicted_0"] >= 3
                    and metrics["predicted_1"] >= 3
                )

                near_perfect = (
                    score_eligible
                    and metrics[
                        "balanced_accuracy"
                    ] is not None
                    and metrics[
                        "balanced_accuracy"
                    ] >= 0.90
                    and metrics["MCC"] >= 0.80
                )

                row = {
                    "panel": panel_name,
                    "trait": trait,
                    "scheme": scheme_name,
                    "chrom": record["chrom"],
                    "pos": record["pos"],
                    "ref": record["ref"],
                    "alt": record["alt"],
                    "comparative_cohort_n": (
                        cohort_total
                    ),
                    "phenotype_0_total": total_0,
                    "phenotype_1_total": total_1,
                    "phenotype_0_called": called_0,
                    "phenotype_1_called": called_1,
                    "call_rate": call_rate,
                    "phenotype_0_call_rate": (
                        call_rate_0
                    ),
                    "phenotype_1_call_rate": (
                        call_rate_1
                    ),
                    "call_rate_state_difference": (
                        abs(call_rate_0 - call_rate_1)
                        if (
                            call_rate_0 is not None
                            and call_rate_1 is not None
                        )
                        else None
                    ),
                    "missingness_fisher_p": (
                        fisher_two_sided(
                            called_1,
                            missing_1,
                            called_0,
                            missing_0,
                        )
                    ),
                    "selected_orientation": (
                        metrics["orientation"]
                    ),
                    "reference_expectation": (
                        panel[
                            "reference_expectation"
                        ]
                    ),
                    "matches_reference_expectation": (
                        metrics["orientation"]
                        == panel[
                            "reference_expectation"
                        ]
                    ),
                    "predicted_0": (
                        metrics["predicted_0"]
                    ),
                    "predicted_1": (
                        metrics["predicted_1"]
                    ),
                    "tp": metrics["tp"],
                    "tn": metrics["tn"],
                    "fp": metrics["fp"],
                    "fn": metrics["fn"],
                    "misclassified_n": (
                        metrics["misclassified_n"]
                    ),
                    "sensitivity": (
                        metrics["sensitivity"]
                    ),
                    "specificity": (
                        metrics["specificity"]
                    ),
                    "balanced_accuracy": (
                        metrics[
                            "balanced_accuracy"
                        ]
                    ),
                    "MCC": metrics["MCC"],
                    "fisher_p": metrics["fisher_p"],
                    "BH_q": None,
                    "missingness_BH_q": None,
                    "callability_eligible": (
                        callability_eligible
                    ),
                    "score_eligible": score_eligible,
                    "near_perfect": near_perfect,
                }

                score_rows.append(row)

    bh_adjust(
        score_rows,
        "fisher_p",
        "BH_q",
    )
    bh_adjust(
        score_rows,
        "missingness_fisher_p",
        "missingness_BH_q",
    )

    score_fieldnames = list(score_rows[0].keys())

    with SCORE_OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=score_fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(score_rows)

    cohort_fieldnames = list(cohort_rows[0].keys())

    with COHORT_OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=cohort_fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(cohort_rows)

    eligibility_rows = []

    grouped = defaultdict(list)

    for row in score_rows:
        grouped[
            (
                row["panel"],
                row["trait"],
                row["scheme"],
            )
        ].append(row)

    for (
        panel,
        trait,
        scheme,
    ), rows in sorted(grouped.items()):
        eligibility_rows.append(
            {
                "panel": panel,
                "trait": trait,
                "scheme": scheme,
                "biallelic_SNPs": len(rows),
                "callability_eligible": sum(
                    row["callability_eligible"]
                    for row in rows
                ),
                "score_eligible": sum(
                    row["score_eligible"]
                    for row in rows
                ),
                "balanced_accuracy_ge_0.80": sum(
                    row["score_eligible"]
                    and row[
                        "balanced_accuracy"
                    ] is not None
                    and row[
                        "balanced_accuracy"
                    ] >= 0.80
                    for row in rows
                ),
                "balanced_accuracy_ge_0.90": sum(
                    row["score_eligible"]
                    and row[
                        "balanced_accuracy"
                    ] is not None
                    and row[
                        "balanced_accuracy"
                    ] >= 0.90
                    for row in rows
                ),
                "near_perfect": sum(
                    row["near_perfect"]
                    for row in rows
                ),
                "BH_q_le_0.05": sum(
                    row["score_eligible"]
                    and row["BH_q"] is not None
                    and row["BH_q"] <= 0.05
                    for row in rows
                ),
                "missingness_q_le_0.05": sum(
                    row["missingness_BH_q"]
                    is not None
                    and row[
                        "missingness_BH_q"
                    ] <= 0.05
                    for row in rows
                ),
            }
        )

    eligibility_fieldnames = list(
        eligibility_rows[0].keys()
    )

    with ELIGIBILITY_OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=eligibility_fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(eligibility_rows)

    print("RECONSTRUCTED OPTIX COMPARATIVE SCORING")
    print("\nSCORING COHORTS")
    print(
        "panel\ttrait\tVCF_samples\tcomparative_n\t"
        "phenotype_0\tphenotype_1\tfocal_excluded\t"
        "unmatched_or_missing"
    )

    for row in cohort_rows:
        print(
            f"{row['panel']}\t"
            f"{row['trait']}\t"
            f"{row['VCF_samples']}\t"
            f"{row['comparative_scoring_samples']}\t"
            f"{row['phenotype_0_n']}\t"
            f"{row['phenotype_1_n']}\t"
            f"{row['focal_samples_excluded']}\t"
            f"{row['unmatched_or_missing_phenotype']}"
        )

    print("\nUSABLE SITE SUMMARY")
    print(
        "panel\ttrait\tscheme\tSNPs\tcallable\t"
        "score_eligible\tBA_ge_0.80\tBA_ge_0.90\t"
        "near_perfect\tq_le_0.05\t"
        "missingness_q_le_0.05"
    )

    for row in eligibility_rows:
        print(
            f"{row['panel']}\t"
            f"{row['trait']}\t"
            f"{row['scheme']}\t"
            f"{row['biallelic_SNPs']}\t"
            f"{row['callability_eligible']}\t"
            f"{row['score_eligible']}\t"
            f"{row['balanced_accuracy_ge_0.80']}\t"
            f"{row['balanced_accuracy_ge_0.90']}\t"
            f"{row['near_perfect']}\t"
            f"{row['BH_q_le_0.05']}\t"
            f"{row['missingness_q_le_0.05']}"
        )

    print("\nTOP SCORE-ELIGIBLE SITES")

    top_rows = sorted(
        (
            row
            for row in score_rows
            if row["score_eligible"]
        ),
        key=lambda row: (
            row["scheme"] != "DP2_GQ20",
            -row["balanced_accuracy"],
            -row["MCC"],
            row["fisher_p"],
        ),
    )[:25]

    print(
        "panel\ttrait\tscheme\tposition\torientation\t"
        "called_0\tcalled_1\tBA\tMCC\t"
        "errors\tp\tq\tmissingness_q"
    )

    for row in top_rows:
        print(
            f"{row['panel']}\t"
            f"{row['trait']}\t"
            f"{row['scheme']}\t"
            f"{row['chrom']}:{row['pos']}\t"
            f"{row['selected_orientation']}\t"
            f"{row['phenotype_0_called']}\t"
            f"{row['phenotype_1_called']}\t"
            f"{row['balanced_accuracy']:.4f}\t"
            f"{row['MCC']:.4f}\t"
            f"{row['misclassified_n']}\t"
            f"{row['fisher_p']:.4g}\t"
            f"{row['BH_q']:.4g}\t"
            f"{row['missingness_BH_q']:.4g}"
        )

    print(f"\nScores: {SCORE_OUTPUT}")
    print(f"Eligibility: {ELIGIBILITY_OUTPUT}")
    print(f"Cohorts: {COHORT_OUTPUT}")
    print(
        "\nFocal leg genotypes and pooled genotypes "
        "were not used in any score."
    )


if __name__ == "__main__":
    main()
