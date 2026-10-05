#!/usr/bin/env python3

import csv
import re
import statistics
import subprocess
from pathlib import Path

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
CHECKPOINT_DIR = PROJECT / "analysis" / "checkpoints"

PANELS = {
    "broad_51_sample_optix": (
        PROJECT
        / "data"
        / "raw"
        / "Optix_Hmr18-938195-1239421.vcf"
    ),
    "filtered_502_sample_optix": (
        PROJECT
        / "data"
        / "raw"
        / "region_panel"
        / "optix_filtered.vcf.gz"
    ),
}

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

POOLED_SIDES = [
    "TN_Left",
    "TN_Right",
]

SUMMARY_OUTPUT = (
    CHECKPOINT_DIR
    / "candidate_plot_specific_optix_slice_summary.tsv"
)

CODE_OUTPUT = (
    CHECKPOINT_DIR
    / "plot_specific_optix_plotting_commands.tsv"
)

RELEVANT_LINES = (
    CHECKPOINT_DIR
    / "original_plotting_relevant_lines.tsv"
)


def run_capture(command):
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    return result.returncode, result.stdout, result.stderr


def canonical_gt(value):
    if not value or "." in value:
        return None

    alleles = re.split(r"[/|]", value)

    try:
        alleles = sorted(int(allele) for allele in alleles)
    except ValueError:
        alleles = sorted(alleles)

    return "/".join(str(allele) for allele in alleles)


def get_samples(vcf):
    code, stdout, stderr = run_capture(
        ["bcftools", "query", "-l", str(vcf)]
    )

    if code != 0:
        return [], stderr.strip()

    samples = [
        line.strip()
        for line in stdout.splitlines()
        if line.strip()
    ]

    return samples, ""


def query_interval(vcf, region):
    command = [
        "bcftools",
        "query",
        "-r",
        region,
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
        str(vcf),
    ]

    code, stdout, stderr = run_capture(command)

    if code != 0:
        return [], stderr.strip()

    records = []

    for line in stdout.splitlines():
        fields = line.split("\t")

        if len(fields) < 4:
            continue

        records.append(
            {
                "chrom": fields[0],
                "pos": int(fields[1]),
                "ref": fields[2],
                "alt": fields[3],
                "genotypes": fields[4:],
            }
        )

    return records, ""


def summarize_interval(panel, vcf, interval_name, interval):
    chrom, start, end = interval
    region = f"{chrom}:{start}-{end}"

    samples, sample_error = get_samples(vcf)

    if sample_error:
        return {
            "panel": panel,
            "input_path": vcf.relative_to(PROJECT).as_posix(),
            "interval_name": interval_name,
            "region": region,
            "status": "sample_error",
            "error": sample_error[:300],
        }

    records, record_error = query_interval(vcf, region)

    if record_error:
        return {
            "panel": panel,
            "input_path": vcf.relative_to(PROJECT).as_posix(),
            "interval_name": interval_name,
            "region": region,
            "status": "query_error",
            "error": record_error[:300],
        }

    sample_index = {
        sample: index
        for index, sample in enumerate(samples)
    }

    focal_present = [
        sample for sample in FOCAL_LEGS
        if sample in sample_index
    ]

    pooled_present = [
        sample for sample in POOLED_SIDES
        if sample in sample_index
    ]

    call_rates = []
    biallelic_snps = 0

    left_both_called = 0
    left_legs_agree = 0
    right_both_called = 0
    right_legs_agree = 0
    both_sides_internally_agree = 0
    consistent_left_right_differences = 0

    pooled_both_called = 0
    pooled_left_right_differences = 0

    for record in records:
        genotypes = [
            canonical_gt(gt)
            for gt in record["genotypes"]
        ]

        if (
            len(record["ref"]) == 1
            and len(record["alt"]) == 1
            and "," not in record["alt"]
        ):
            biallelic_snps += 1

        called = sum(gt is not None for gt in genotypes)

        if samples:
            call_rates.append(called / len(samples))

        if all(sample in sample_index for sample in FOCAL_LEGS):
            left_1 = genotypes[sample_index["TN2L_S1"]]
            right_1 = genotypes[sample_index["TN2R_S2"]]
            left_2 = genotypes[sample_index["TN3L_S3"]]
            right_2 = genotypes[sample_index["TN3R_S4"]]

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

        if all(sample in sample_index for sample in POOLED_SIDES):
            pooled_left = genotypes[sample_index["TN_Left"]]
            pooled_right = genotypes[sample_index["TN_Right"]]

            if (
                pooled_left is not None
                and pooled_right is not None
            ):
                pooled_both_called += 1

                if pooled_left != pooled_right:
                    pooled_left_right_differences += 1

    return {
        "panel": panel,
        "input_path": vcf.relative_to(PROJECT).as_posix(),
        "interval_name": interval_name,
        "region": region,
        "status": "ok",
        "samples": len(samples),
        "records": len(records),
        "biallelic_snps": biallelic_snps,
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
        "focal_legs_present": ",".join(focal_present),
        "pooled_sides_present": ",".join(pooled_present),
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
        "pooled_both_called": pooled_both_called,
        "pooled_left_right_differences": (
            pooled_left_right_differences
        ),
        "error": "",
    }


def extract_plotting_commands():
    if not RELEVANT_LINES.is_file():
        return []

    selected = []

    pattern = re.compile(
        r"""(?ix)
        genotype_plot
        |our_popmap
        |read\.vcfR
        |metadata\s*<-
        |metadata\$
        |filter\(
        |mutate\(
        |arrange\(
        |pop\s*=
        |ind\s*=
        |facet
        |sort
        |order
        |colour
        |color
        """
    )

    with RELEVANT_LINES.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            if (
                row.get("source")
                != "gynandromorph_genotype_plots.html"
            ):
                continue

            text = row.get("text", "")

            if pattern.search(text):
                selected.append(
                    {
                        "source": row.get("source", ""),
                        "location": row.get("location", ""),
                        "text": text,
                    }
                )

    with CODE_OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source", "location", "text"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(selected)

    return selected


def main():
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    missing = [
        str(vcf)
        for vcf in PANELS.values()
        if not vcf.is_file()
    ]

    if missing:
        print("ERROR: missing Optix inputs:")
        for path in missing:
            print(path)
        return

    rows = []

    for panel, vcf in PANELS.items():
        for interval_name, interval in INTERVALS.items():
            rows.append(
                summarize_interval(
                    panel,
                    vcf,
                    interval_name,
                    interval,
                )
            )

    fieldnames = [
        "panel",
        "input_path",
        "interval_name",
        "region",
        "status",
        "samples",
        "records",
        "biallelic_snps",
        "call_rate_min",
        "call_rate_median",
        "call_rate_max",
        "focal_legs_present",
        "pooled_sides_present",
        "left_both_called",
        "left_legs_agree",
        "right_both_called",
        "right_legs_agree",
        "both_sides_internally_agree",
        "consistent_left_right_differences",
        "pooled_both_called",
        "pooled_left_right_differences",
        "error",
    ]

    with SUMMARY_OUTPUT.open(
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

    commands = extract_plotting_commands()

    print("CANDIDATE PLOT-SPECIFIC OPTIX SLICE AUDIT")
    print(
        "panel\tinterval\trecords\tbiallelic_SNPs\t"
        "call_rate_median\tlegs\tpooled\t"
        "internally_consistent_LR_differences\t"
        "pooled_LR_differences"
    )

    for row in rows:
        print(
            f"{row['panel']}\t"
            f"{row['interval_name']}\t"
            f"{row.get('records', '')}\t"
            f"{row.get('biallelic_snps', '')}\t"
            f"{row.get('call_rate_median', '')}\t"
            f"{row.get('focal_legs_present', '')}\t"
            f"{row.get('pooled_sides_present', '')}\t"
            f"{row.get('consistent_left_right_differences', '')}\t"
            f"{row.get('pooled_left_right_differences', '')}"
        )

    print("\nPLOT-CONSTRUCTION COMMANDS")
    print(f"Relevant lines: {len(commands)}")

    for row in commands:
        print(
            f"{row['location']}\t{row['text']}"
        )

    print(f"\nSlice summary: {SUMMARY_OUTPUT}")
    print(f"Plot command table: {CODE_OUTPUT}")


if __name__ == "__main__":
    main()
