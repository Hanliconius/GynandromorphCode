#!/usr/bin/env python3

import argparse
import csv
import glob
import itertools
import json
import os
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chromosome", type=int, required=True)
    parser.add_argument("--bcftools", required=True)
    parser.add_argument("--call-root", required=True)
    parser.add_argument("--sample-file", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--minimum-depth", type=int, default=5)
    parser.add_argument("--minimum-pl-gap", type=int, default=20)
    return parser.parse_args()


def parse_genotype(field):
    parts = field.split(":", 3)

    while len(parts) < 4:
        parts.append(".")

    gt_text, dp_text, ad_text, pl_text = parts

    if gt_text in {".", "./.", ".|."}:
        return None

    alleles = gt_text.replace("|", "/").split("/")

    if len(alleles) != 2:
        return None

    try:
        allele_numbers = [int(value) for value in alleles]
    except ValueError:
        return None

    if any(value not in {0, 1} for value in allele_numbers):
        return None

    try:
        depth = int(dp_text)
    except ValueError:
        return None

    try:
        allele_depths = [int(value) for value in ad_text.split(",")]
    except ValueError:
        return None

    try:
        likelihoods = [int(value) for value in pl_text.split(",")]
    except ValueError:
        return None

    if len(allele_depths) != 2 or len(likelihoods) < 3:
        return None

    likelihoods = likelihoods[:3]
    dosage = sum(allele_numbers)
    ordered = sorted((value, index) for index, value in enumerate(likelihoods))

    if len(ordered) < 2:
        return None

    return {
        "dosage": dosage,
        "depth": depth,
        "ad_sum": sum(allele_depths),
        "best_pl_genotype": ordered[0][1],
        "pl_gap": ordered[1][0] - ordered[0][0],
    }


def confident_genotype(parsed, minimum_depth, minimum_pl_gap):
    if parsed is None:
        return False

    return (
        parsed["depth"] >= minimum_depth
        and parsed["ad_sum"] > 0
        and parsed["best_pl_genotype"] == parsed["dosage"]
        and parsed["pl_gap"] >= minimum_pl_gap
    )


def possible_child_dosages(mother_dosage, father_dosage):
    genotype_alleles = {
        0: (0, 0),
        1: (0, 1),
        2: (1, 1),
    }

    mother = genotype_alleles[mother_dosage]
    father = genotype_alleles[father_dosage]

    return {
        maternal_allele + paternal_allele
        for maternal_allele in mother
        for paternal_allele in father
    }


def write_tsv(path, fieldnames, rows):
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


def main():
    args = parse_args()
    started = time.time()

    if args.chromosome < 1 or args.chromosome > 20:
        raise RuntimeError(
            "Stage 196 is autosomal only: chromosome must be 1-20."
        )

    sample_file = Path(args.sample_file)

    if not sample_file.is_file() or sample_file.stat().st_size == 0:
        raise RuntimeError(f"Missing sample file: {sample_file}")

    samples = [
        line.strip()
        for line in sample_file.read_text().splitlines()
        if line.strip()
    ]

    if len(samples) != 14:
        raise RuntimeError(
            f"Expected 14 pedigree samples; found {len(samples)}"
        )

    mother = samples[0]
    father = samples[1]
    offspring = samples[2:]

    if len(offspring) != 12:
        raise RuntimeError(
            f"Expected 12 offspring; found {len(offspring)}"
        )

    pattern = str(
        Path(args.call_root)
        / f"chr{args.chromosome}_*"
        / "*.bcf"
    )

    bcf_candidates = [
        Path(path)
        for path in glob.glob(pattern)
        if not path.endswith(".bcf.new")
    ]

    if len(bcf_candidates) != 1:
        raise RuntimeError(
            "Expected exactly one completed family BCF for "
            f"chromosome {args.chromosome}; found "
            f"{len(bcf_candidates)}: {bcf_candidates}"
        )

    bcf = bcf_candidates[0]
    chromosome_directory = bcf.parent
    output_directory = (
        Path(args.output_root)
        / chromosome_directory.name
    )
    output_directory.mkdir(parents=True, exist_ok=True)

    observed_samples = subprocess.check_output(
        [args.bcftools, "query", "-l", str(bcf)],
        text=True,
    ).splitlines()

    if observed_samples != samples:
        raise RuntimeError(
            "BCF sample order differs from the validated pedigree sample file."
        )

    offspring_counts = {
        sample: Counter()
        for sample in offspring
    }

    pair_counts = {
        pair: Counter()
        for pair in itertools.combinations(offspring, 2)
    }

    site_counts = Counter()

    query_format = (
        "%CHROM\\t%POS\\t%REF\\t%ALT"
        "[\\t%GT:%DP:%AD:%PL]\\n"
    )

    command = [
        args.bcftools,
        "query",
        "-f",
        query_format,
        str(bcf),
    ]

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )

    assert process.stdout is not None

    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        site_counts["queried"] += 1

        if len(fields) != 4 + len(samples):
            site_counts["malformed_sample_fields"] += 1
            continue

        ref = fields[2]
        alt = fields[3]

        if (
            len(ref) != 1
            or len(alt) != 1
            or ref not in "ACGT"
            or alt not in "ACGT"
            or "," in alt
        ):
            site_counts["non_biallelic_snp"] += 1
            continue

        site_counts["biallelic_snp"] += 1

        parsed = [
            parse_genotype(value)
            for value in fields[4:]
        ]

        confident = [
            confident_genotype(
                value,
                args.minimum_depth,
                args.minimum_pl_gap,
            )
            for value in parsed
        ]

        if not confident[0] or not confident[1]:
            site_counts["parent_qc_failure"] += 1
            continue

        mother_dosage = parsed[0]["dosage"]
        father_dosage = parsed[1]["dosage"]
        allowed = possible_child_dosages(
            mother_dosage,
            father_dosage,
        )

        site_counts["parents_confident"] += 1

        if len(allowed) < 3:
            site_counts["parent_informative"] += 1

        offspring_dosages = {}

        for index, sample in enumerate(offspring, start=2):
            counts = offspring_counts[sample]
            counts["parent_confident_sites"] += 1

            if not confident[index]:
                counts["offspring_qc_failure"] += 1
                continue

            dosage = parsed[index]["dosage"]
            offspring_dosages[sample] = dosage
            counts["evaluated"] += 1
            counts[f"offspring_GT_{dosage}"] += 1

            if dosage in allowed:
                counts["mendelian_consistent"] += 1
            else:
                counts["mendelian_error"] += 1

        for first, second in pair_counts:
            counts = pair_counts[(first, second)]

            if (
                first not in offspring_dosages
                or second not in offspring_dosages
            ):
                counts["pair_qc_failure"] += 1
                continue

            first_dosage = offspring_dosages[first]
            second_dosage = offspring_dosages[second]

            counts["evaluated"] += 1

            if first_dosage == second_dosage:
                counts["same_genotype"] += 1

            if abs(first_dosage - second_dosage) == 1:
                counts["dosage_distance_one"] += 1

            if {first_dosage, second_dosage} == {0, 2}:
                counts["opposing_homozygotes"] += 1

            if first_dosage == 1 and second_dosage == 1:
                counts["both_heterozygous"] += 1

    stderr_text = ""

    if process.stderr is not None:
        stderr_text = process.stderr.read()

    return_code = process.wait()

    if return_code != 0:
        if stderr_text:
            print(stderr_text, file=sys.stderr)
        raise RuntimeError(
            f"bcftools query failed with status {return_code}"
        )

    offspring_rows = []

    for sample in offspring:
        counts = offspring_counts[sample]
        evaluated = counts["evaluated"]
        errors = counts["mendelian_error"]

        offspring_rows.append(
            {
                "chromosome": args.chromosome,
                "mother": mother,
                "father": father,
                "offspring": sample,
                "parent_confident_sites":
                    counts["parent_confident_sites"],
                "evaluated_sites": evaluated,
                "mendelian_consistent":
                    counts["mendelian_consistent"],
                "mendelian_errors": errors,
                "mendelian_error_rate":
                    errors / evaluated if evaluated else "",
                "offspring_qc_failures":
                    counts["offspring_qc_failure"],
                "GT_0": counts["offspring_GT_0"],
                "GT_1": counts["offspring_GT_1"],
                "GT_2": counts["offspring_GT_2"],
            }
        )

    pair_rows = []

    for first, second in itertools.combinations(offspring, 2):
        counts = pair_counts[(first, second)]
        evaluated = counts["evaluated"]
        opposing = counts["opposing_homozygotes"]

        pair_rows.append(
            {
                "chromosome": args.chromosome,
                "offspring_a": first,
                "offspring_b": second,
                "evaluated_sites": evaluated,
                "same_genotype": counts["same_genotype"],
                "dosage_distance_one":
                    counts["dosage_distance_one"],
                "both_heterozygous":
                    counts["both_heterozygous"],
                "opposing_homozygotes": opposing,
                "opposing_homozygote_rate":
                    opposing / evaluated if evaluated else "",
                "pair_qc_failures":
                    counts["pair_qc_failure"],
            }
        )

    offspring_path = (
        output_directory
        / "offspring_mendelian_qc.tsv"
    )
    pair_path = (
        output_directory
        / "sibling_pair_genotype_states.tsv"
    )
    audit_path = (
        output_directory
        / "mendelian_qc.audit.json"
    )
    completion_path = (
        output_directory
        / "mendelian_qc.complete.tsv"
    )

    write_tsv(
        offspring_path,
        [
            "chromosome",
            "mother",
            "father",
            "offspring",
            "parent_confident_sites",
            "evaluated_sites",
            "mendelian_consistent",
            "mendelian_errors",
            "mendelian_error_rate",
            "offspring_qc_failures",
            "GT_0",
            "GT_1",
            "GT_2",
        ],
        offspring_rows,
    )

    write_tsv(
        pair_path,
        [
            "chromosome",
            "offspring_a",
            "offspring_b",
            "evaluated_sites",
            "same_genotype",
            "dosage_distance_one",
            "both_heterozygous",
            "opposing_homozygotes",
            "opposing_homozygote_rate",
            "pair_qc_failures",
        ],
        pair_rows,
    )

    audit = {
        "chromosome": args.chromosome,
        "input_bcf": str(bcf),
        "mother": mother,
        "father": father,
        "offspring": offspring,
        "minimum_depth": args.minimum_depth,
        "minimum_PL_gap": args.minimum_pl_gap,
        "site_counts": dict(site_counts),
        "offspring_rows": len(offspring_rows),
        "sibling_pairs": len(pair_rows),
        "elapsed_seconds": round(time.time() - started, 3),
        "scope": (
            "autosomal pedigree and sibling-pair diagnostic; "
            "no founder-model calibration performed"
        ),
    }

    temporary_audit = Path(str(audit_path) + ".tmp")
    temporary_audit.write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    os.replace(temporary_audit, audit_path)

    write_tsv(
        completion_path,
        [
            "chromosome",
            "status",
            "offspring",
            "sibling_pairs",
            "queried_sites",
            "evaluated_parent_sites",
            "elapsed_seconds",
        ],
        [
            {
                "chromosome": args.chromosome,
                "status": "PASS",
                "offspring": len(offspring_rows),
                "sibling_pairs": len(pair_rows),
                "queried_sites": site_counts["queried"],
                "evaluated_parent_sites":
                    site_counts["parents_confident"],
                "elapsed_seconds":
                    round(time.time() - started, 3),
            }
        ],
    )

    print(
        "MENDELIAN_QC_COMPLETE"
        f"\tchromosome={args.chromosome}"
        f"\tqueried={site_counts['queried']}"
        f"\tparents_confident={site_counts['parents_confident']}"
        f"\toffspring={len(offspring_rows)}"
        f"\tsibling_pairs={len(pair_rows)}"
    )


if __name__ == "__main__":
    main()
