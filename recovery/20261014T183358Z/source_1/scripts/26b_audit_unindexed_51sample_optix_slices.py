#!/usr/bin/env python3

import csv
import re
import statistics
import subprocess
from pathlib import Path

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

VCF = (
    PROJECT
    / "data"
    / "raw"
    / "Optix_Hmr18-938195-1239421.vcf"
)

OUTPUT = (
    PROJECT
    / "analysis"
    / "checkpoints"
    / "candidate_51sample_optix_slice_summary.tsv"
)

INTERVALS = {
    "band": ("Hmr18", 1112651, 1114314),
    "dennis1": ("Hmr18", 1215210, 1216124),
    "dennis2": ("Hmr18", 1218557, 1221494),
    "ray1": ("Hmr18", 1156746, 1160220),
    "ray2": ("Hmr18", 1174547, 1187833),
}

FOCAL_LEGS = [
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
]


def canonical_gt(value):
    if not value or "." in value:
        return None

    alleles = re.split(r"[/|]", value)

    try:
        alleles = sorted(int(allele) for allele in alleles)
    except ValueError:
        alleles = sorted(alleles)

    return "/".join(str(allele) for allele in alleles)


def get_samples():
    result = subprocess.run(
        ["bcftools", "query", "-l", str(VCF)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        return [], result.stderr.strip()

    return [
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip()
    ], ""


def find_interval(chrom, position):
    for name, (target_chrom, start, end) in INTERVALS.items():
        if (
            chrom == target_chrom
            and start <= position <= end
        ):
            return name

    return None


def main():
    if not VCF.is_file():
        print(f"ERROR: missing VCF: {VCF}")
        return

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)

    samples, error = get_samples()

    if error:
        print(f"ERROR: could not read samples: {error}")
        return

    sample_index = {
        sample: index
        for index, sample in enumerate(samples)
    }

    focal_present = [
        sample for sample in FOCAL_LEGS
        if sample in sample_index
    ]

    records = {
        name: []
        for name in INTERVALS
    }

    process = subprocess.Popen(
        [
            "bcftools",
            "query",
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
            str(VCF),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if process.stdout is not None:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 4:
                continue

            chrom = fields[0]

            try:
                position = int(fields[1])
            except ValueError:
                continue

            interval_name = find_interval(chrom, position)

            if interval_name is None:
                continue

            records[interval_name].append(
                {
                    "chrom": chrom,
                    "pos": position,
                    "ref": fields[2],
                    "alt": fields[3],
                    "genotypes": [
                        canonical_gt(gt)
                        for gt in fields[4:]
                    ],
                }
            )

    stderr = ""

    if process.stderr is not None:
        stderr = process.stderr.read().strip()

    return_code = process.wait()

    if return_code != 0:
        print("ERROR: bcftools query failed.")
        if stderr:
            print(stderr)
        return

    rows = []

    for interval_name, interval in INTERVALS.items():
        chrom, start, end = interval
        interval_records = records[interval_name]

        call_rates = []
        biallelic_snps = 0
        polymorphic_hard_gt_sites = 0
        focal_variable_sites = 0

        left_both_called = 0
        left_legs_agree = 0
        right_both_called = 0
        right_legs_agree = 0
        both_sides_internally_agree = 0
        consistent_left_right_differences = 0

        for record in interval_records:
            genotypes = record["genotypes"]

            if (
                len(record["ref"]) == 1
                and len(record["alt"]) == 1
                and "," not in record["alt"]
            ):
                biallelic_snps += 1

            called_genotypes = [
                gt for gt in genotypes
                if gt is not None
            ]

            if samples:
                call_rates.append(
                    len(called_genotypes) / len(samples)
                )

            if len(set(called_genotypes)) > 1:
                polymorphic_hard_gt_sites += 1

            if all(
                sample in sample_index
                for sample in FOCAL_LEGS
            ):
                left_1 = genotypes[
                    sample_index["TN2L_S1"]
                ]
                right_1 = genotypes[
                    sample_index["TN2R_S2"]
                ]
                left_2 = genotypes[
                    sample_index["TN3L_S3"]
                ]
                right_2 = genotypes[
                    sample_index["TN3R_S4"]
                ]

                focal_called = [
                    gt
                    for gt in [
                        left_1,
                        right_1,
                        left_2,
                        right_2,
                    ]
                    if gt is not None
                ]

                if len(set(focal_called)) > 1:
                    focal_variable_sites += 1

                left_called = (
                    left_1 is not None
                    and left_2 is not None
                )
                right_called = (
                    right_1 is not None
                    and right_2 is not None
                )

                if left_called:
                    left_both_called += 1

                    if left_1 == left_2:
                        left_legs_agree += 1

                if right_called:
                    right_both_called += 1

                    if right_1 == right_2:
                        right_legs_agree += 1

                if (
                    left_called
                    and right_called
                    and left_1 == left_2
                    and right_1 == right_2
                ):
                    both_sides_internally_agree += 1

                    if left_1 != right_1:
                        consistent_left_right_differences += 1

        rows.append(
            {
                "interval": interval_name,
                "region": f"{chrom}:{start}-{end}",
                "samples": len(samples),
                "records": len(interval_records),
                "biallelic_snps": biallelic_snps,
                "polymorphic_hard_GT_sites": (
                    polymorphic_hard_gt_sites
                ),
                "focal_variable_sites": focal_variable_sites,
                "call_rate_min": (
                    min(call_rates) if call_rates else ""
                ),
                "call_rate_median": (
                    statistics.median(call_rates)
                    if call_rates else ""
                ),
                "call_rate_max": (
                    max(call_rates) if call_rates else ""
                ),
                "left_both_called": left_both_called,
                "left_legs_agree": left_legs_agree,
                "right_both_called": right_both_called,
                "right_legs_agree": right_legs_agree,
                "both_sides_internally_agree": (
                    both_sides_internally_agree
                ),
                "consistent_left_right_differences": (
                    consistent_left_right_differences
                ),
            }
        )

    fieldnames = list(rows[0].keys())

    with OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print("UNINDEXED 51-SAMPLE OPTIX SLICE AUDIT")
    print(f"Input samples: {len(samples)}")
    print(
        "Focal samples present: "
        + ",".join(focal_present)
    )

    print(
        "\ninterval\trecords\tbiallelic_SNPs\t"
        "polymorphic_sites\tfocal_variable_sites\t"
        "call_rate_median\tleft_agree/called\t"
        "right_agree/called\tboth_sides_agree\t"
        "consistent_LR_differences"
    )

    for row in rows:
        print(
            f"{row['interval']}\t"
            f"{row['records']}\t"
            f"{row['biallelic_snps']}\t"
            f"{row['polymorphic_hard_GT_sites']}\t"
            f"{row['focal_variable_sites']}\t"
            f"{row['call_rate_median']}\t"
            f"{row['left_legs_agree']}/"
            f"{row['left_both_called']}\t"
            f"{row['right_legs_agree']}/"
            f"{row['right_both_called']}\t"
            f"{row['both_sides_internally_agree']}\t"
            f"{row['consistent_left_right_differences']}"
        )

    print(f"\nSummary: {OUTPUT}")


if __name__ == "__main__":
    main()
