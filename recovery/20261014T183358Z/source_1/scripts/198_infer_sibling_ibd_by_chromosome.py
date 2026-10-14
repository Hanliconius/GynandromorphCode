#!/usr/bin/env python3

import argparse
import csv
import glob
import itertools
import json
import math
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chromosome", type=int, required=True)
    parser.add_argument("--bcftools", required=True)
    parser.add_argument("--call-root", required=True)
    parser.add_argument("--sample-file", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--window-size", type=int, default=1_000_000)
    parser.add_argument("--minimum-depth", type=int, default=5)
    parser.add_argument("--minimum-pl-gap", type=int, default=20)
    parser.add_argument("--classification-error", type=float, default=0.01)
    parser.add_argument("--minimum-informative-sites", type=int, default=10)
    parser.add_argument("--minimum-absolute-llr", type=float, default=5.0)
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
        depth = int(dp_text)
        allele_depths = [
            int(value)
            for value in ad_text.split(",")
        ]
        likelihoods = [
            int(value)
            for value in pl_text.split(",")
        ]
    except ValueError:
        return None

    if any(value not in {0, 1} for value in allele_numbers):
        return None

    if len(allele_depths) != 2 or len(likelihoods) < 3:
        return None

    dosage = sum(allele_numbers)
    ordered_pl = sorted(
        (value, index)
        for index, value in enumerate(likelihoods[:3])
    )

    return {
        "dosage": dosage,
        "depth": depth,
        "ad_sum": sum(allele_depths),
        "best_pl_genotype": ordered_pl[0][1],
        "pl_gap": ordered_pl[1][0] - ordered_pl[0][0],
    }


def is_confident(parsed, minimum_depth, minimum_pl_gap):
    if parsed is None:
        return False

    return (
        parsed["depth"] >= minimum_depth
        and parsed["ad_sum"] > 0
        and parsed["best_pl_genotype"] == parsed["dosage"]
        and parsed["pl_gap"] >= minimum_pl_gap
    )


def transmitted_allele(
    heterozygous_parent_dosage,
    homozygous_parent_dosage,
    offspring_dosage,
):
    if heterozygous_parent_dosage != 1:
        return None

    if homozygous_parent_dosage == 0:
        transmitted = offspring_dosage
    elif homozygous_parent_dosage == 2:
        transmitted = offspring_dosage - 1
    else:
        return None

    if transmitted not in {0, 1}:
        return None

    return transmitted


def classify_sharing(
    same_count,
    different_count,
    error_rate,
    minimum_sites,
    minimum_absolute_llr,
):
    total = same_count + different_count

    if total < minimum_sites:
        return {
            "sites": total,
            "same_fraction": (
                same_count / total if total else ""
            ),
            "llr_shared_vs_not": "",
            "state": "insufficient",
        }

    log_shared = (
        same_count * math.log(1.0 - error_rate)
        + different_count * math.log(error_rate)
    )

    log_not_shared = (
        same_count * math.log(error_rate)
        + different_count * math.log(1.0 - error_rate)
    )

    llr = log_shared - log_not_shared

    if llr >= minimum_absolute_llr:
        state = "shared"
    elif llr <= -minimum_absolute_llr:
        state = "not_shared"
    else:
        state = "uncertain"

    return {
        "sites": total,
        "same_fraction": same_count / total,
        "llr_shared_vs_not": llr,
        "state": state,
    }


def provisional_ibd(maternal_state, paternal_state):
    resolved = {"shared": 1, "not_shared": 0}

    if (
        maternal_state not in resolved
        or paternal_state not in resolved
    ):
        return ""

    return (
        resolved[maternal_state]
        + resolved[paternal_state]
    )


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
            "Stage 198 is autosomal only: chromosome must be 1-20."
        )

    if not 0.0 < args.classification_error < 0.5:
        raise RuntimeError(
            "classification-error must be between 0 and 0.5."
        )

    sample_file = Path(args.sample_file)

    if not sample_file.is_file() or sample_file.stat().st_size == 0:
        raise RuntimeError(
            f"Missing pedigree sample file: {sample_file}"
        )

    samples = [
        line.strip()
        for line in sample_file.read_text().splitlines()
        if line.strip()
    ]

    if len(samples) != 14:
        raise RuntimeError(
            f"Expected 14 samples; found {len(samples)}"
        )

    mother = samples[0]
    father = samples[1]
    offspring = samples[2:]
    pairs = list(itertools.combinations(offspring, 2))

    if len(offspring) != 12 or len(pairs) != 66:
        raise RuntimeError(
            "Expected 12 offspring and 66 sibling pairs."
        )

    bcf_pattern = str(
        Path(args.call_root)
        / f"chr{args.chromosome}_*"
        / "*.bcf"
    )

    bcf_candidates = [
        Path(path)
        for path in glob.glob(bcf_pattern)
        if not path.endswith(".bcf.new")
    ]

    if len(bcf_candidates) != 1:
        raise RuntimeError(
            "Expected one family BCF for chromosome "
            f"{args.chromosome}; found {bcf_candidates}"
        )

    bcf = bcf_candidates[0]
    contig = bcf.parent.name.split("_", 1)[1]

    output_directory = (
        Path(args.output_root)
        / bcf.parent.name
    )
    output_directory.mkdir(parents=True, exist_ok=True)

    observed_samples = subprocess.check_output(
        [args.bcftools, "query", "-l", str(bcf)],
        text=True,
    ).splitlines()

    if observed_samples != samples:
        raise RuntimeError(
            "BCF sample order differs from pedigree sample order."
        )

    window_counts = defaultdict(
        lambda: {
            pair: {
                "maternal_same": 0,
                "maternal_different": 0,
                "paternal_same": 0,
                "paternal_different": 0,
            }
            for pair in pairs
        }
    )

    chromosome_counts = {
        pair: Counter()
        for pair in pairs
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
            site_counts["malformed"] += 1
            continue

        position = int(fields[1])
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

        parsed = [
            parse_genotype(value)
            for value in fields[4:]
        ]

        confident = [
            is_confident(
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

        maternal_informative = (
            mother_dosage == 1
            and father_dosage in {0, 2}
        )

        paternal_informative = (
            father_dosage == 1
            and mother_dosage in {0, 2}
        )

        if not maternal_informative and not paternal_informative:
            site_counts["parent_uninformative"] += 1
            continue

        if maternal_informative:
            site_counts["maternal_informative"] += 1

        if paternal_informative:
            site_counts["paternal_informative"] += 1

        transmitted = {}

        for index, sample in enumerate(offspring, start=2):
            if not confident[index]:
                continue

            offspring_dosage = parsed[index]["dosage"]
            maternal_allele = None
            paternal_allele = None

            if maternal_informative:
                maternal_allele = transmitted_allele(
                    mother_dosage,
                    father_dosage,
                    offspring_dosage,
                )

            if paternal_informative:
                paternal_allele = transmitted_allele(
                    father_dosage,
                    mother_dosage,
                    offspring_dosage,
                )

            transmitted[sample] = (
                maternal_allele,
                paternal_allele,
            )

        window_start = (
            ((position - 1) // args.window_size)
            * args.window_size
        )
        window_key = window_start

        for first, second in pairs:
            if first not in transmitted or second not in transmitted:
                continue

            first_maternal, first_paternal = transmitted[first]
            second_maternal, second_paternal = transmitted[second]

            window_pair = window_counts[window_key][
                (first, second)
            ]
            chromosome_pair = chromosome_counts[
                (first, second)
            ]

            if (
                first_maternal is not None
                and second_maternal is not None
            ):
                if first_maternal == second_maternal:
                    window_pair["maternal_same"] += 1
                    chromosome_pair["maternal_same"] += 1
                else:
                    window_pair["maternal_different"] += 1
                    chromosome_pair["maternal_different"] += 1

            if (
                first_paternal is not None
                and second_paternal is not None
            ):
                if first_paternal == second_paternal:
                    window_pair["paternal_same"] += 1
                    chromosome_pair["paternal_same"] += 1
                else:
                    window_pair["paternal_different"] += 1
                    chromosome_pair["paternal_different"] += 1

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

    window_rows = []

    for window_start in sorted(window_counts):
        window_end = window_start + args.window_size

        for first, second in pairs:
            counts = window_counts[window_start][
                (first, second)
            ]

            maternal = classify_sharing(
                counts["maternal_same"],
                counts["maternal_different"],
                args.classification_error,
                args.minimum_informative_sites,
                args.minimum_absolute_llr,
            )

            paternal = classify_sharing(
                counts["paternal_same"],
                counts["paternal_different"],
                args.classification_error,
                args.minimum_informative_sites,
                args.minimum_absolute_llr,
            )

            window_rows.append(
                {
                    "chromosome": args.chromosome,
                    "contig": contig,
                    "window_start": window_start,
                    "window_end": window_end,
                    "offspring_a": first,
                    "offspring_b": second,
                    "maternal_sites": maternal["sites"],
                    "maternal_same":
                        counts["maternal_same"],
                    "maternal_different":
                        counts["maternal_different"],
                    "maternal_same_fraction":
                        maternal["same_fraction"],
                    "maternal_LLR_shared_vs_not":
                        maternal["llr_shared_vs_not"],
                    "maternal_state": maternal["state"],
                    "paternal_sites": paternal["sites"],
                    "paternal_same":
                        counts["paternal_same"],
                    "paternal_different":
                        counts["paternal_different"],
                    "paternal_same_fraction":
                        paternal["same_fraction"],
                    "paternal_LLR_shared_vs_not":
                        paternal["llr_shared_vs_not"],
                    "paternal_state": paternal["state"],
                    "provisional_IBD":
                        provisional_ibd(
                            maternal["state"],
                            paternal["state"],
                        ),
                }
            )

    chromosome_rows = []

    for first, second in pairs:
        counts = chromosome_counts[(first, second)]

        maternal = classify_sharing(
            counts["maternal_same"],
            counts["maternal_different"],
            args.classification_error,
            args.minimum_informative_sites,
            args.minimum_absolute_llr,
        )

        paternal = classify_sharing(
            counts["paternal_same"],
            counts["paternal_different"],
            args.classification_error,
            args.minimum_informative_sites,
            args.minimum_absolute_llr,
        )

        chromosome_rows.append(
            {
                "chromosome": args.chromosome,
                "contig": contig,
                "offspring_a": first,
                "offspring_b": second,
                "maternal_sites": maternal["sites"],
                "maternal_same":
                    counts["maternal_same"],
                "maternal_different":
                    counts["maternal_different"],
                "maternal_same_fraction":
                    maternal["same_fraction"],
                "maternal_LLR_shared_vs_not":
                    maternal["llr_shared_vs_not"],
                "maternal_state": maternal["state"],
                "paternal_sites": paternal["sites"],
                "paternal_same":
                    counts["paternal_same"],
                "paternal_different":
                    counts["paternal_different"],
                "paternal_same_fraction":
                    paternal["same_fraction"],
                "paternal_LLR_shared_vs_not":
                    paternal["llr_shared_vs_not"],
                "paternal_state": paternal["state"],
                "provisional_IBD":
                    provisional_ibd(
                        maternal["state"],
                        paternal["state"],
                    ),
            }
        )

    window_output = (
        output_directory
        / "sibling_transmission_IBD_windows.tsv"
    )
    chromosome_output = (
        output_directory
        / "sibling_transmission_IBD_chromosome.tsv"
    )
    audit_output = (
        output_directory
        / "sibling_transmission_IBD.audit.json"
    )
    completion_output = (
        output_directory
        / "sibling_transmission_IBD.complete.tsv"
    )

    common_fields = [
        "chromosome",
        "contig",
        "offspring_a",
        "offspring_b",
        "maternal_sites",
        "maternal_same",
        "maternal_different",
        "maternal_same_fraction",
        "maternal_LLR_shared_vs_not",
        "maternal_state",
        "paternal_sites",
        "paternal_same",
        "paternal_different",
        "paternal_same_fraction",
        "paternal_LLR_shared_vs_not",
        "paternal_state",
        "provisional_IBD",
    ]

    write_tsv(
        window_output,
        [
            "chromosome",
            "contig",
            "window_start",
            "window_end",
            "offspring_a",
            "offspring_b",
            "maternal_sites",
            "maternal_same",
            "maternal_different",
            "maternal_same_fraction",
            "maternal_LLR_shared_vs_not",
            "maternal_state",
            "paternal_sites",
            "paternal_same",
            "paternal_different",
            "paternal_same_fraction",
            "paternal_LLR_shared_vs_not",
            "paternal_state",
            "provisional_IBD",
        ],
        window_rows,
    )

    write_tsv(
        chromosome_output,
        common_fields,
        chromosome_rows,
    )

    audit = {
        "chromosome": args.chromosome,
        "contig": contig,
        "input_bcf": str(bcf),
        "mother": mother,
        "father": father,
        "offspring": len(offspring),
        "sibling_pairs": len(pairs),
        "window_size": args.window_size,
        "minimum_depth": args.minimum_depth,
        "minimum_PL_gap": args.minimum_pl_gap,
        "classification_error":
            args.classification_error,
        "minimum_informative_sites":
            args.minimum_informative_sites,
        "minimum_absolute_LLR":
            args.minimum_absolute_llr,
        "site_counts": dict(site_counts),
        "window_rows": len(window_rows),
        "chromosome_rows": len(chromosome_rows),
        "scope": (
            "parental-transmission IBD diagnostic; "
            "unsmoothed provisional classifications"
        ),
        "elapsed_seconds":
            round(time.time() - started, 3),
    }

    temporary_audit = Path(str(audit_output) + ".tmp")
    temporary_audit.write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    os.replace(temporary_audit, audit_output)

    write_tsv(
        completion_output,
        [
            "chromosome",
            "status",
            "sibling_pairs",
            "window_rows",
            "maternal_informative_sites",
            "paternal_informative_sites",
            "elapsed_seconds",
        ],
        [
            {
                "chromosome": args.chromosome,
                "status": "PASS",
                "sibling_pairs": len(pairs),
                "window_rows": len(window_rows),
                "maternal_informative_sites":
                    site_counts["maternal_informative"],
                "paternal_informative_sites":
                    site_counts["paternal_informative"],
                "elapsed_seconds":
                    round(time.time() - started, 3),
            }
        ],
    )

    print(
        "SIBLING_IBD_COMPLETE"
        f"\tchromosome={args.chromosome}"
        f"\twindows={len(window_counts)}"
        f"\trows={len(window_rows)}"
        "\tmaternal_informative="
        f"{site_counts['maternal_informative']}"
        "\tpaternal_informative="
        f"{site_counts['paternal_informative']}"
    )


if __name__ == "__main__":
    main()
