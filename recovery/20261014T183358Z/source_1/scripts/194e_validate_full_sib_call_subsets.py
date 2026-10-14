#!/usr/bin/env python3

import argparse
import csv
import gzip
import itertools
import subprocess
from pathlib import Path


ROOT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/"
    "GynandroCalling"
)

BCFTOOLS = Path(
    "/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/"
    ".local/share/mamba/envs/colias_w_popgen/bin/bcftools"
)

SOURCE_ROOT = (
    ROOT
    / "analysis/genomewide_mixing_input/variant_calls/"
    "joint_24/by_chromosome"
)

CALL_ROOT = (
    ROOT
    / "analysis/genomewide_mixing_corrected_v1/"
    "full_sib_validation/variant_calls/by_chromosome"
)

SAMPLE_FILE = (
    ROOT
    / "analysis/genomewide_mixing_corrected_v1/"
    "full_sib_validation/metadata/bam_set/"
    "full_sib_14_samples.samples"
)

AUDIT_DIRECTORY = (
    ROOT
    / "analysis/genomewide_mixing_corrected_v1/"
    "full_sib_validation/audit"
)


def one_matching_directory(root, chromosome):
    matches = sorted(root.glob(f"chr{chromosome}_*"))

    if len(matches) != 1:
        raise RuntimeError(
            f"Chromosome {chromosome}: expected one directory under "
            f"{root}; found {len(matches)}"
        )

    return matches[0]


def target_records(path, counts):
    previous_key = None

    with gzip.open(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 3:
                raise RuntimeError(
                    f"{path}:{line_number}: expected three columns; "
                    f"found {len(fields)}"
                )

            chrom, position_text, allele_token = fields

            if "," not in allele_token:
                raise RuntimeError(
                    f"{path}:{line_number}: column 3 is not REF,ALT: "
                    f"{allele_token!r}"
                )

            position = int(position_text)
            key = (chrom, position)

            if previous_key is not None and key < previous_key:
                raise RuntimeError(
                    f"{path} is not coordinate sorted: "
                    f"{previous_key} followed by {key}"
                )

            previous_key = key
            counts["target_sites"] += 1

            yield chrom, position, allele_token


def target_groups(path, counts):
    records = target_records(path, counts)

    for key, group in itertools.groupby(
        records,
        key=lambda record: (record[0], record[1]),
    ):
        allele_tokens = {record[2] for record in group}
        yield key, allele_tokens


def contig_from_directory(directory, chromosome):
    prefix = f"chr{chromosome}_"

    if not directory.name.startswith(prefix):
        raise RuntimeError(
            f"Unexpected chromosome directory: {directory}"
        )

    return directory.name[len(prefix):]


def validate_chromosome(task_id):
    chromosome = task_id + 1

    source_directory = one_matching_directory(
        SOURCE_ROOT,
        chromosome,
    )
    call_directory = one_matching_directory(
        CALL_ROOT,
        chromosome,
    )

    contig = contig_from_directory(
        call_directory,
        chromosome,
    )

    alleles = (
        source_directory
        / "controls.discovery.alleles.tsv.gz"
    )

    bcf = (
        call_directory
        / "full_sib_14.joint_at_existing_control_sites.bcf"
    )

    index = Path(str(bcf) + ".csi")
    completion = call_directory / "call.complete.tsv"

    for path in (
        alleles,
        bcf,
        index,
        SAMPLE_FILE,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(
                f"Missing or empty required input: {path}"
            )

    expected_samples = SAMPLE_FILE.read_text().splitlines()

    sample_result = subprocess.run(
        [
            str(BCFTOOLS),
            "query",
            "-l",
            str(bcf),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    observed_samples = sample_result.stdout.splitlines()

    if observed_samples != expected_samples:
        raise RuntimeError(
            f"Chromosome {chromosome}: sample order differs "
            "from the BAM-set sample order"
        )

    if len(observed_samples) != 14:
        raise RuntimeError(
            f"Chromosome {chromosome}: expected 14 samples; "
            f"found {len(observed_samples)}"
        )

    index_result = subprocess.run(
        [
            str(BCFTOOLS),
            "index",
            "-n",
            str(bcf),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    indexed_sites = int(index_result.stdout.strip())

    counts = {
        "target_sites": 0,
        "emitted_sites": 0,
        "off_target_sites": 0,
    }

    grouped_targets = target_groups(
        alleles,
        counts,
    )

    current_target = next(grouped_targets, None)
    previous_bcf_key = None
    off_target_examples = []

    query_process = subprocess.Popen(
        [
            str(BCFTOOLS),
            "query",
            "-f",
            r"%CHROM\t%POS\t%REF\t%ALT\n",
            str(bcf),
        ],
        stdout=subprocess.PIPE,
        text=True,
    )

    if query_process.stdout is None:
        raise RuntimeError(
            "Failed to read bcftools query output"
        )

    for line_number, line in enumerate(
        query_process.stdout,
        start=1,
    ):
        fields = line.rstrip("\n").split("\t")

        if len(fields) != 4:
            query_process.kill()

            raise RuntimeError(
                f"BCF query line {line_number}: expected four "
                f"columns; found {len(fields)}"
            )

        chrom, position_text, ref, alt = fields
        position = int(position_text)
        bcf_key = (chrom, position)
        allele_token = f"{ref},{alt}"

        if (
            previous_bcf_key is not None
            and bcf_key < previous_bcf_key
        ):
            query_process.kill()

            raise RuntimeError(
                f"BCF is not coordinate sorted: "
                f"{previous_bcf_key} followed by {bcf_key}"
            )

        previous_bcf_key = bcf_key
        counts["emitted_sites"] += 1

        while (
            current_target is not None
            and current_target[0] < bcf_key
        ):
            current_target = next(
                grouped_targets,
                None,
            )

        valid = (
            current_target is not None
            and current_target[0] == bcf_key
            and allele_token in current_target[1]
        )

        if not valid:
            counts["off_target_sites"] += 1

            if len(off_target_examples) < 10:
                off_target_examples.append(
                    f"{chrom}:{position}:{allele_token}"
                )

    query_return_code = query_process.wait()

    if query_return_code != 0:
        raise RuntimeError(
            f"bcftools query failed with status "
            f"{query_return_code}"
        )

    # Finish reading the target so target_sites is the complete count.
    for _ in grouped_targets:
        pass

    if counts["emitted_sites"] != indexed_sites:
        raise RuntimeError(
            f"Chromosome {chromosome}: query returned "
            f"{counts['emitted_sites']} sites but the BCF index "
            f"reports {indexed_sites}"
        )

    if counts["emitted_sites"] > counts["target_sites"]:
        raise RuntimeError(
            f"Chromosome {chromosome}: emitted-site count "
            "exceeds target-site count"
        )

    if counts["off_target_sites"] != 0:
        examples = ", ".join(off_target_examples)

        raise RuntimeError(
            f"Chromosome {chromosome}: found "
            f"{counts['off_target_sites']} emitted sites absent "
            f"from the canonical CHROM/POS/REF,ALT target. "
            f"Examples: {examples}"
        )

    retained_fraction = (
        counts["emitted_sites"]
        / counts["target_sites"]
    )

    temporary = Path(str(completion) + ".new")

    with temporary.open("w", newline="") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            [
                "chromosome",
                "contig",
                "samples",
                "target_sites",
                "emitted_sites",
                "retained_fraction",
                "off_target_sites",
                "subset_status",
                "bcf",
                "bcf_bytes",
            ]
        )

        writer.writerow(
            [
                chromosome,
                contig,
                len(observed_samples),
                counts["target_sites"],
                counts["emitted_sites"],
                f"{retained_fraction:.12f}",
                counts["off_target_sites"],
                "PASS",
                str(bcf),
                bcf.stat().st_size,
            ]
        )

    temporary.replace(completion)

    print(
        "PASS"
        f"\tchromosome={chromosome}"
        f"\tcontig={contig}"
        f"\ttarget_sites={counts['target_sites']}"
        f"\temitted_sites={counts['emitted_sites']}"
        f"\tretained_fraction={retained_fraction:.6f}"
        "\toff_target_sites=0"
    )


def aggregate_validation():
    rows = []

    for chromosome in range(1, 22):
        call_directory = one_matching_directory(
            CALL_ROOT,
            chromosome,
        )

        completion = (
            call_directory
            / "call.complete.tsv"
        )

        bcf = (
            call_directory
            / "full_sib_14.joint_at_existing_control_sites.bcf"
        )

        index = Path(str(bcf) + ".csi")

        for path in (
            completion,
            bcf,
            index,
        ):
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError(
                    f"Missing or empty required file: {path}"
                )

        with completion.open() as handle:
            completion_rows = list(
                csv.DictReader(
                    handle,
                    delimiter="\t",
                )
            )

        if len(completion_rows) != 1:
            raise RuntimeError(
                f"{completion}: expected exactly one data row"
            )

        row = completion_rows[0]

        if row.get("subset_status") != "PASS":
            raise RuntimeError(
                f"Chromosome {chromosome}: subset validation "
                "did not pass"
            )

        if int(row["samples"]) != 14:
            raise RuntimeError(
                f"Chromosome {chromosome}: expected 14 samples"
            )

        if int(row["off_target_sites"]) != 0:
            raise RuntimeError(
                f"Chromosome {chromosome}: off-target sites found"
            )

        if (
            int(row["emitted_sites"])
            > int(row["target_sites"])
        ):
            raise RuntimeError(
                f"Chromosome {chromosome}: emitted count "
                "exceeds target count"
            )

        rows.append(row)

    AUDIT_DIRECTORY.mkdir(
        parents=True,
        exist_ok=True,
    )

    report = (
        AUDIT_DIRECTORY
        / "full_sib_joint_call_validation.tsv"
    )

    report_temporary = Path(str(report) + ".new")

    with report_temporary.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=list(rows[0].keys()),
        )

        writer.writeheader()
        writer.writerows(rows)

    report_temporary.replace(report)

    total_target = sum(
        int(row["target_sites"])
        for row in rows
    )

    total_emitted = sum(
        int(row["emitted_sites"])
        for row in rows
    )

    retained_fraction = (
        total_emitted / total_target
    )

    marker = (
        AUDIT_DIRECTORY
        / "full_sib_joint_calls.complete.tsv"
    )

    marker_temporary = Path(str(marker) + ".new")

    marker_temporary.write_text(
        "item\tvalue\n"
        "chromosomes\t21\n"
        "samples\t14\n"
        f"target_sites\t{total_target}\n"
        f"emitted_sites\t{total_emitted}\n"
        f"retained_fraction\t{retained_fraction:.12f}\n"
        "off_target_sites\t0\n"
        "status\tPASS\n"
    )

    marker_temporary.replace(marker)

    print("FULL-SIB JOINT-CALL VALIDATION PASSED")
    print("chromosomes\t21")
    print("samples\t14")
    print(f"target_sites\t{total_target}")
    print(f"emitted_sites\t{total_emitted}")
    print(
        f"retained_fraction\t"
        f"{retained_fraction:.6f}"
    )
    print("off_target_sites\t0")
    print(f"report\t{report}")


def main():
    parser = argparse.ArgumentParser()

    mode = parser.add_mutually_exclusive_group(
        required=True,
    )

    mode.add_argument(
        "--task-id",
        type=int,
    )

    mode.add_argument(
        "--aggregate",
        action="store_true",
    )

    arguments = parser.parse_args()

    if arguments.aggregate:
        aggregate_validation()
    else:
        if not 0 <= arguments.task_id <= 20:
            raise RuntimeError(
                "--task-id must be between 0 and 20"
            )

        validate_chromosome(
            arguments.task_id
        )


if __name__ == "__main__":
    main()
