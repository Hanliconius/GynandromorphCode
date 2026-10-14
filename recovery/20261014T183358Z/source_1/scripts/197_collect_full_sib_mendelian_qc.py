#!/usr/bin/env python3

import csv
import glob
import json
import math
import os
import statistics
import time
from collections import defaultdict
from pathlib import Path


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/"
    "GynandroCalling"
)

INPUT_ROOT = (
    PROJECT
    / "analysis/genomewide_mixing_corrected_v1/"
    "full_sib_validation/mendelian_qc/by_chromosome"
)

OUTPUT_ROOT = (
    PROJECT
    / "analysis/genomewide_mixing_corrected_v1/"
    "full_sib_validation/mendelian_qc/genomewide"
)


OFFSPRING_INTEGER_FIELDS = [
    "parent_confident_sites",
    "evaluated_sites",
    "mendelian_consistent",
    "mendelian_errors",
    "offspring_qc_failures",
    "GT_0",
    "GT_1",
    "GT_2",
]

PAIR_INTEGER_FIELDS = [
    "evaluated_sites",
    "same_genotype",
    "dosage_distance_one",
    "both_heterozygous",
    "opposing_homozygotes",
    "pair_qc_failures",
]


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fieldnames, rows):
    path = Path(path)
    temporary = Path(str(path) + ".tmp")

    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    os.replace(temporary, path)


def robust_z(value, values):
    median = statistics.median(values)
    deviations = [
        abs(item - median)
        for item in values
    ]
    mad = statistics.median(deviations)

    if mad == 0:
        return 0.0

    return (value - median) / (1.4826 * mad)


def safe_rate(numerator, denominator):
    if denominator == 0:
        return ""
    return numerator / denominator


def main():
    started = time.time()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    completion_files = glob.glob(
        str(
            INPUT_ROOT
            / "chr*_*"
            / "mendelian_qc.complete.tsv"
        )
    )

    offspring_files = glob.glob(
        str(
            INPUT_ROOT
            / "chr*_*"
            / "offspring_mendelian_qc.tsv"
        )
    )

    pair_files = glob.glob(
        str(
            INPUT_ROOT
            / "chr*_*"
            / "sibling_pair_genotype_states.tsv"
        )
    )

    if len(completion_files) != 20:
        raise RuntimeError(
            "Expected 20 chromosome completion files; found "
            f"{len(completion_files)}"
        )

    if len(offspring_files) != 20:
        raise RuntimeError(
            "Expected 20 offspring-QC files; found "
            f"{len(offspring_files)}"
        )

    if len(pair_files) != 20:
        raise RuntimeError(
            "Expected 20 sibling-pair files; found "
            f"{len(pair_files)}"
        )

    completed_chromosomes = set()

    for path in completion_files:
        rows = read_tsv(path)

        if len(rows) != 1:
            raise RuntimeError(
                f"Malformed completion file: {path}"
            )

        row = rows[0]

        if row["status"] != "PASS":
            raise RuntimeError(
                f"Non-PASS chromosome completion: {path}"
            )

        completed_chromosomes.add(
            int(row["chromosome"])
        )

    expected_chromosomes = set(range(1, 21))

    if completed_chromosomes != expected_chromosomes:
        raise RuntimeError(
            "Completed chromosome set is incorrect: "
            f"{sorted(completed_chromosomes)}"
        )

    offspring_totals = defaultdict(
        lambda: defaultdict(int)
    )
    offspring_chromosomes = defaultdict(set)

    for path in offspring_files:
        for row in read_tsv(path):
            sample = row["offspring"]
            chromosome = int(row["chromosome"])
            offspring_chromosomes[sample].add(chromosome)

            for field in OFFSPRING_INTEGER_FIELDS:
                offspring_totals[sample][field] += int(
                    row[field]
                )

    if len(offspring_totals) != 12:
        raise RuntimeError(
            "Expected 12 distinct offspring; found "
            f"{len(offspring_totals)}"
        )

    for sample, chromosomes in offspring_chromosomes.items():
        if chromosomes != expected_chromosomes:
            raise RuntimeError(
                f"Incomplete chromosome coverage for {sample}: "
                f"{sorted(chromosomes)}"
            )

    offspring_error_rates = []

    for sample in sorted(offspring_totals):
        counts = offspring_totals[sample]
        rate = safe_rate(
            counts["mendelian_errors"],
            counts["evaluated_sites"],
        )

        if rate != "":
            offspring_error_rates.append(rate)

    median_error_rate = statistics.median(
        offspring_error_rates
    )

    error_mad = statistics.median(
        [
            abs(rate - median_error_rate)
            for rate in offspring_error_rates
        ]
    )

    review_threshold = max(
        0.005,
        median_error_rate + 6.0 * 1.4826 * error_mad,
    )

    offspring_rows = []

    for sample in sorted(offspring_totals):
        counts = offspring_totals[sample]
        evaluated = counts["evaluated_sites"]
        error_rate = safe_rate(
            counts["mendelian_errors"],
            evaluated,
        )

        z_score = (
            robust_z(
                error_rate,
                offspring_error_rates,
            )
            if error_rate != ""
            else ""
        )

        review_flag = (
            "REVIEW"
            if (
                error_rate != ""
                and error_rate > review_threshold
            )
            else "PASS"
        )

        offspring_rows.append(
            {
                "offspring": sample,
                "autosomes": 20,
                "parent_confident_sites":
                    counts["parent_confident_sites"],
                "evaluated_sites": evaluated,
                "mendelian_consistent":
                    counts["mendelian_consistent"],
                "mendelian_errors":
                    counts["mendelian_errors"],
                "mendelian_error_rate": error_rate,
                "robust_error_z": z_score,
                "review_threshold":
                    review_threshold,
                "review_flag": review_flag,
                "offspring_qc_failures":
                    counts["offspring_qc_failures"],
                "GT_0": counts["GT_0"],
                "GT_1": counts["GT_1"],
                "GT_2": counts["GT_2"],
                "heterozygous_fraction":
                    safe_rate(counts["GT_1"], evaluated),
            }
        )

    pair_totals = defaultdict(
        lambda: defaultdict(int)
    )
    pair_chromosomes = defaultdict(set)

    for path in pair_files:
        for row in read_tsv(path):
            pair = tuple(
                sorted(
                    (
                        row["offspring_a"],
                        row["offspring_b"],
                    )
                )
            )
            chromosome = int(row["chromosome"])
            pair_chromosomes[pair].add(chromosome)

            for field in PAIR_INTEGER_FIELDS:
                pair_totals[pair][field] += int(
                    row[field]
                )

    if len(pair_totals) != 66:
        raise RuntimeError(
            "Expected 66 sibling pairs; found "
            f"{len(pair_totals)}"
        )

    for pair, chromosomes in pair_chromosomes.items():
        if chromosomes != expected_chromosomes:
            raise RuntimeError(
                "Incomplete chromosome coverage for pair "
                f"{pair}: {sorted(chromosomes)}"
            )

    maximum_pair_sites = max(
        counts["evaluated_sites"]
        for counts in pair_totals.values()
    )

    opposing_rates = [
        counts["opposing_homozygotes"]
        / counts["evaluated_sites"]
        for counts in pair_totals.values()
        if counts["evaluated_sites"] > 0
    ]

    median_opposing_rate = statistics.median(
        opposing_rates
    )

    median_pair_sites = statistics.median(
        [
            counts["evaluated_sites"]
            for counts in pair_totals.values()
        ]
    )

    provisional_pair_rows = []

    for pair in sorted(pair_totals):
        counts = pair_totals[pair]
        evaluated = counts["evaluated_sites"]

        opposing_rate = safe_rate(
            counts["opposing_homozygotes"],
            evaluated,
        )

        same_rate = safe_rate(
            counts["same_genotype"],
            evaluated,
        )

        distance_one_rate = safe_rate(
            counts["dosage_distance_one"],
            evaluated,
        )

        both_heterozygous_rate = safe_rate(
            counts["both_heterozygous"],
            evaluated,
        )

        coverage_fraction = (
            evaluated / maximum_pair_sites
            if maximum_pair_sites
            else 0.0
        )

        contrast_score = (
            opposing_rate
            if opposing_rate != ""
            else 0.0
        )

        combined_information_score = (
            math.sqrt(
                coverage_fraction * contrast_score
            )
            if contrast_score > 0
            else 0.0
        )

        information_tier = (
            "high_information"
            if (
                evaluated >= median_pair_sites
                and opposing_rate != ""
                and opposing_rate >= median_opposing_rate
            )
            else "standard_information"
        )

        provisional_pair_rows.append(
            {
                "offspring_a": pair[0],
                "offspring_b": pair[1],
                "autosomes": 20,
                "evaluated_sites": evaluated,
                "coverage_fraction_of_max":
                    coverage_fraction,
                "same_genotype":
                    counts["same_genotype"],
                "same_genotype_rate": same_rate,
                "dosage_distance_one":
                    counts["dosage_distance_one"],
                "dosage_distance_one_rate":
                    distance_one_rate,
                "both_heterozygous":
                    counts["both_heterozygous"],
                "both_heterozygous_rate":
                    both_heterozygous_rate,
                "opposing_homozygotes":
                    counts["opposing_homozygotes"],
                "opposing_homozygote_rate":
                    opposing_rate,
                "pair_qc_failures":
                    counts["pair_qc_failures"],
                "combined_information_score":
                    combined_information_score,
                "information_tier":
                    information_tier,
            }
        )

    ranked_rows = sorted(
        provisional_pair_rows,
        key=lambda row: (
            -row["combined_information_score"],
            -row["evaluated_sites"],
            row["offspring_a"],
            row["offspring_b"],
        ),
    )

    for rank, row in enumerate(ranked_rows, start=1):
        row["information_rank"] = rank

    reviewed_samples = [
        row["offspring"]
        for row in offspring_rows
        if row["review_flag"] == "REVIEW"
    ]

    high_information_pairs = sum(
        row["information_tier"] == "high_information"
        for row in ranked_rows
    )

    offspring_output = (
        OUTPUT_ROOT
        / "genomewide_offspring_mendelian_qc.tsv"
    )
    pair_output = (
        OUTPUT_ROOT
        / "genomewide_sibling_pair_genotype_states.tsv"
    )
    report_output = (
        OUTPUT_ROOT
        / "pedigree_qc_report.txt"
    )
    audit_output = (
        OUTPUT_ROOT
        / "pedigree_qc.audit.json"
    )
    completion_output = (
        OUTPUT_ROOT
        / "pedigree_qc.complete.tsv"
    )

    write_tsv(
        offspring_output,
        [
            "offspring",
            "autosomes",
            "parent_confident_sites",
            "evaluated_sites",
            "mendelian_consistent",
            "mendelian_errors",
            "mendelian_error_rate",
            "robust_error_z",
            "review_threshold",
            "review_flag",
            "offspring_qc_failures",
            "GT_0",
            "GT_1",
            "GT_2",
            "heterozygous_fraction",
        ],
        offspring_rows,
    )

    write_tsv(
        pair_output,
        [
            "information_rank",
            "offspring_a",
            "offspring_b",
            "autosomes",
            "evaluated_sites",
            "coverage_fraction_of_max",
            "same_genotype",
            "same_genotype_rate",
            "dosage_distance_one",
            "dosage_distance_one_rate",
            "both_heterozygous",
            "both_heterozygous_rate",
            "opposing_homozygotes",
            "opposing_homozygote_rate",
            "pair_qc_failures",
            "combined_information_score",
            "information_tier",
        ],
        ranked_rows,
    )

    report_lines = [
        "FULL-SIB PEDIGREE QC SUMMARY",
        "",
        "Scope:",
        "  Autosomes 1-20",
        "  Diagnostic aggregation only",
        "  No samples or pairs excluded",
        "  No founder model changed",
        "",
        f"Chromosomes completed: {len(completed_chromosomes)}",
        f"Offspring evaluated: {len(offspring_rows)}",
        f"Sibling pairs evaluated: {len(ranked_rows)}",
        (
            "Median offspring Mendelian error rate: "
            f"{median_error_rate:.8f}"
        ),
        (
            "Diagnostic review threshold: "
            f"{review_threshold:.8f}"
        ),
        (
            "Offspring flagged for review: "
            + (
                ", ".join(reviewed_samples)
                if reviewed_samples
                else "none"
            )
        ),
        (
            "High-information sibling pairs: "
            f"{high_information_pairs}"
        ),
        "",
        "Highest-ranked empirical sibling pairs:",
    ]

    for row in ranked_rows[:10]:
        report_lines.append(
            "  "
            f"{row['information_rank']}. "
            f"{row['offspring_a']} / "
            f"{row['offspring_b']} "
            f"sites={row['evaluated_sites']} "
            "opposing_homozygote_rate="
            f"{row['opposing_homozygote_rate']:.6f}"
        )

    report_lines.extend(
        [
            "",
            "Interpretation:",
            (
                "  Rankings describe empirical information "
                "content, not pedigree acceptance or rejection."
            ),
            (
                "  Mendelian-review flags are diagnostics and "
                "do not automatically remove samples."
            ),
            (
                "  Formal IBD-state inference and full-sibling "
                "model calibration occur in later stages."
            ),
            "",
        ]
    )

    temporary_report = Path(str(report_output) + ".tmp")
    temporary_report.write_text(
        "\n".join(report_lines)
    )
    os.replace(temporary_report, report_output)

    audit = {
        "chromosomes": 20,
        "offspring": len(offspring_rows),
        "sibling_pairs": len(ranked_rows),
        "median_mendelian_error_rate":
            median_error_rate,
        "mendelian_review_threshold":
            review_threshold,
        "offspring_flagged_for_review":
            reviewed_samples,
        "high_information_pairs":
            high_information_pairs,
        "automatic_exclusions": False,
        "founder_model_changed": False,
        "elapsed_seconds":
            round(time.time() - started, 3),
        "next_stage": (
            "chromosome-scale sibling IBD-state "
            "characterization"
        ),
    }

    temporary_audit = Path(str(audit_output) + ".tmp")
    temporary_audit.write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    os.replace(temporary_audit, audit_output)

    write_tsv(
        completion_output,
        [
            "status",
            "chromosomes",
            "offspring",
            "sibling_pairs",
            "reviewed_offspring",
            "elapsed_seconds",
        ],
        [
            {
                "status": "PASS",
                "chromosomes": 20,
                "offspring": len(offspring_rows),
                "sibling_pairs": len(ranked_rows),
                "reviewed_offspring":
                    len(reviewed_samples),
                "elapsed_seconds":
                    round(time.time() - started, 3),
            }
        ],
    )

    print("\n".join(report_lines))
    print(
        "ENDPOINT="
        "FULL_SIB_PEDIGREE_QC_AGGREGATED"
    )


if __name__ == "__main__":
    main()
