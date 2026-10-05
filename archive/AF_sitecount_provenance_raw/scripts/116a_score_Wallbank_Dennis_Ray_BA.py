#!/usr/bin/env python3

from pathlib import Path
import csv
import math
import re
from collections import Counter


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

VCF = (
    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
    / "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

WALLBANK_MANIFEST = (
    PROJECT
    / "analysis"
    / "original_502"
    / "fixed_SNP_replication"
    / "Wallbank2016_fixed_SNP_replication_samples.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis"
    / "original_502"
    / "Wallbank_Dennis_Ray_BA"
)

SCORE_FILE = OUTDIR / "Wallbank_Dennis_Ray_dominant_BA.tsv"
SAMPLE_FILE = OUTDIR / "Wallbank_Dennis_Ray_BA_samples.tsv"
FOCAL_FILE = OUTDIR / "Wallbank_Dennis_Ray_focal_calls.tsv"
SUMMARY_FILE = OUTDIR / "Wallbank_Dennis_Ray_BA.summary.tsv"
KNOWN_FILE = OUTDIR / "Wallbank_Dennis_Ray_known_sites.tsv"

CHROM = "Hmr18"
START = 1060000
END = 1240000

FOCAL_SAMPLES = ("TN_Left", "TN_Right")

RAY_ABSENT = [
    "ERR3653537",
    "ERR3653539",
    "ERR1143606",
    "ERR1143605",
    "ERR3653533",
]

RAY_PRESENT = [
    "ERR3653298",
    "ERR3653300",
    "ERR3653301",
    "ERR3653302",
    "ERR3653303",
    "ERR3653305",
    "ERR3653306",
    "ERR3653307",
    "ERR3653309",
    "ERR3653543",
]

KNOWN_RAY_POSITIONS = {
    1162302,
    1167014,
    1170389,
    1210276,
}


def canonical_gt(value):
    if value is None:
        return None

    value = value.replace("|", "/")

    if value in {".", "./.", ".|."}:
        return None

    alleles = value.split("/")

    if len(alleles) != 2:
        return None

    if any(allele not in {"0", "1"} for allele in alleles):
        return None

    return "/".join(sorted(alleles))


def numeric(value):
    if value in {None, "", "."}:
        return None

    try:
        return float(value)
    except ValueError:
        return None


def read_header():
    with VCF.open() as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                return fields[9:]

    raise RuntimeError("VCF #CHROM header was not found.")


def extract_dennis_groups(vcf_samples):
    vcf_sample_set = set(vcf_samples)

    with WALLBANK_MANIFEST.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    absent = []
    present = []
    diagnostic_rows = []

    for row in rows:
        values = [
            str(value).strip()
            for value in row.values()
            if value is not None
        ]

        blob = " | ".join(values)
        blob_lower = blob.lower()

        if "dennis_exact_specimen_overlap" not in blob_lower:
            continue

        matching_samples = sorted(
            {
                value
                for value in values
                if value in vcf_sample_set
            }
        )

        if len(matching_samples) != 1:
            continue

        sample = matching_samples[0]

        is_absent = bool(
            re.search(
                r"(^|[^a-z0-9])no[_ -]?d([^a-z0-9]|$)",
                blob_lower,
            )
        )

        is_present = bool(
            re.search(
                r"(^|[^a-z0-9])d[_ -]?present([^a-z0-9]|$)",
                blob_lower,
            )
        )

        if is_absent and not is_present:
            absent.append(sample)
            group = "Dennis_absent"
        elif is_present and not is_absent:
            present.append(sample)
            group = "Dennis_present"
        else:
            continue

        diagnostic_rows.append(
            {
                "module": "Dennis",
                "phenotype": 0 if group == "Dennis_absent" else 1,
                "group": group,
                "sample": sample,
            }
        )

    absent = sorted(set(absent))
    present = sorted(set(present))

    if len(absent) != 4 or len(present) != 2:
        raise RuntimeError(
            "Expected Dennis groups of 4 absent and 2 present; "
            f"recovered {len(absent)} absent and {len(present)} present. "
            "Inspect the Wallbank sample manifest."
        )

    return absent, present, diagnostic_rows


def parse_sample_call(sample_field, format_index):
    values = sample_field.split(":")

    def get_field(name):
        index = format_index.get(name)

        if index is None or index >= len(values):
            return None

        return values[index]

    gt = canonical_gt(get_field("GT"))
    dp = numeric(get_field("DP"))
    gq = numeric(get_field("GQ"))
    ad = get_field("AD")

    return {
        "GT": gt,
        "DP": dp,
        "GQ": gq,
        "AD": ad if ad not in {None, "", "."} else "",
    }


def paper_qc(call):
    return (
        call["GT"] is not None
        and call["DP"] is not None
        and call["GQ"] is not None
        and 4 <= call["DP"] <= 300
        and call["GQ"] >= 30
    )


def focal_display_qc(call):
    return (
        call["GT"] is not None
        and call["DP"] is not None
        and call["GQ"] is not None
        and call["DP"] >= 2
        and call["GQ"] >= 20
    )


def dominant_prediction(gt, associated_allele):
    if gt is None:
        return None

    alleles = gt.split("/")
    return int(associated_allele in alleles)


def evaluate_orientation(
    absent_calls,
    present_calls,
    associated_allele,
):
    absent_predictions = [
        dominant_prediction(call["GT"], associated_allele)
        for call in absent_calls
    ]

    present_predictions = [
        dominant_prediction(call["GT"], associated_allele)
        for call in present_calls
    ]

    absent_correct = sum(
        prediction == 0
        for prediction in absent_predictions
    )

    present_correct = sum(
        prediction == 1
        for prediction in present_predictions
    )

    specificity = (
        absent_correct / len(absent_predictions)
        if absent_predictions
        else None
    )

    sensitivity = (
        present_correct / len(present_predictions)
        if present_predictions
        else None
    )

    balanced_accuracy = (
        (specificity + sensitivity) / 2
        if specificity is not None and sensitivity is not None
        else None
    )

    return {
        "associated_allele": associated_allele,
        "absent_correct": absent_correct,
        "present_correct": present_correct,
        "specificity": specificity,
        "sensitivity": sensitivity,
        "balanced_accuracy": balanced_accuracy,
        "misclassified_n": (
            len(absent_predictions)
            + len(present_predictions)
            - absent_correct
            - present_correct
        ),
    }


def score_module(
    module,
    absent_samples,
    present_samples,
    calls,
):
    absent_qc = [
        calls[sample]
        for sample in absent_samples
        if paper_qc(calls[sample])
    ]

    present_qc = [
        calls[sample]
        for sample in present_samples
        if paper_qc(calls[sample])
    ]

    orientation_alt = evaluate_orientation(
        absent_qc,
        present_qc,
        "1",
    )

    orientation_ref = evaluate_orientation(
        absent_qc,
        present_qc,
        "0",
    )

    alternatives = [
        orientation_alt,
        orientation_ref,
    ]

    alternatives.sort(
        key=lambda row: (
            row["balanced_accuracy"]
            if row["balanced_accuracy"] is not None
            else -1,
            -row["misclassified_n"],
            row["associated_allele"] == "1",
        ),
        reverse=True,
    )

    selected = alternatives[0]

    absent_called = len(absent_qc)
    present_called = len(present_qc)

    absent_n = len(absent_samples)
    present_n = len(present_samples)

    absent_call_rate = absent_called / absent_n
    present_call_rate = present_called / present_n
    minimum_group_call_rate = min(
        absent_call_rate,
        present_call_rate,
    )

    eligible_min2 = (
        absent_called >= 2
        and present_called >= 2
    )

    callable_80pct = (
        absent_called >= math.ceil(0.8 * absent_n)
        and present_called >= math.ceil(0.8 * present_n)
    )

    complete_callability = (
        absent_called == absent_n
        and present_called == present_n
    )

    orientation = (
        "ALT_carrier_predicts_presence"
        if selected["associated_allele"] == "1"
        else "REF_carrier_predicts_presence"
    )

    return {
        "module": module,
        "absent_n": absent_n,
        "present_n": present_n,
        "absent_called": absent_called,
        "present_called": present_called,
        "absent_call_rate": absent_call_rate,
        "present_call_rate": present_call_rate,
        "minimum_group_call_rate": minimum_group_call_rate,
        "selected_orientation": orientation,
        "associated_allele": selected["associated_allele"],
        "specificity": selected["specificity"],
        "sensitivity": selected["sensitivity"],
        "balanced_accuracy": selected["balanced_accuracy"],
        "absent_correct": selected["absent_correct"],
        "present_correct": selected["present_correct"],
        "misclassified_n": selected["misclassified_n"],
        "eligible_min2": eligible_min2,
        "callable_80pct": callable_80pct,
        "complete_callability": complete_callability,
    }


def format_value(value):
    if value is None:
        return ""

    if isinstance(value, bool):
        return str(value)

    return value


def main():
    if not VCF.is_file():
        raise RuntimeError(f"Missing broad Optix VCF: {VCF}")

    if not WALLBANK_MANIFEST.is_file():
        raise RuntimeError(
            f"Missing Wallbank sample manifest: {WALLBANK_MANIFEST}"
        )

    OUTDIR.mkdir(parents=True, exist_ok=True)

    vcf_samples = read_header()
    vcf_sample_set = set(vcf_samples)
    sample_to_index = {
        sample: index
        for index, sample in enumerate(vcf_samples)
    }

    dennis_absent, dennis_present, dennis_manifest = (
        extract_dennis_groups(vcf_samples)
    )

    required_samples = set(
        dennis_absent
        + dennis_present
        + RAY_ABSENT
        + RAY_PRESENT
        + list(FOCAL_SAMPLES)
    )

    missing_samples = sorted(
        required_samples - vcf_sample_set
    )

    if missing_samples:
        raise RuntimeError(
            "Samples missing from VCF: "
            + ", ".join(missing_samples)
        )

    sample_manifest = list(dennis_manifest)

    for sample in RAY_ABSENT:
        sample_manifest.append(
            {
                "module": "Ray",
                "phenotype": 0,
                "group": "Ray_absent",
                "sample": sample,
            }
        )

    for sample in RAY_PRESENT:
        sample_manifest.append(
            {
                "module": "Ray",
                "phenotype": 1,
                "group": "Ray_present_contigua",
                "sample": sample,
            }
        )

    with SAMPLE_FILE.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "module",
                "phenotype",
                "group",
                "sample",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(sample_manifest)

    module_groups = {
        "Dennis": (
            dennis_absent,
            dennis_present,
        ),
        "Ray": (
            RAY_ABSENT,
            RAY_PRESENT,
        ),
    }

    required_order = sorted(required_samples)
    required_indices = {
        sample: sample_to_index[sample]
        for sample in required_order
    }

    score_rows = []
    focal_rows = []
    processed_records = 0

    with VCF.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                continue

            chrom = fields[0]

            if chrom != CHROM:
                continue

            try:
                pos = int(fields[1])
            except ValueError:
                continue

            if pos < START or pos > END:
                continue

            ref = fields[3]
            alt = fields[4]

            if (
                len(ref) != 1
                or len(alt) != 1
                or "," in alt
                or alt in {".", "*"}
            ):
                continue

            format_names = fields[8].split(":")
            format_index = {
                name: index
                for index, name in enumerate(format_names)
            }

            if "GT" not in format_index:
                continue

            calls = {}

            for sample, sample_index in required_indices.items():
                calls[sample] = parse_sample_call(
                    fields[9 + sample_index],
                    format_index,
                )

            left_call = calls["TN_Left"]
            right_call = calls["TN_Right"]

            focal_rows.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "TN_Left_GT": left_call["GT"],
                    "TN_Left_AD": left_call["AD"],
                    "TN_Left_DP": left_call["DP"],
                    "TN_Left_GQ": left_call["GQ"],
                    "TN_Left_DP2_GQ20": focal_display_qc(left_call),
                    "TN_Right_GT": right_call["GT"],
                    "TN_Right_AD": right_call["AD"],
                    "TN_Right_DP": right_call["DP"],
                    "TN_Right_GQ": right_call["GQ"],
                    "TN_Right_DP2_GQ20": focal_display_qc(right_call),
                }
            )

            for module, (absent, present) in module_groups.items():
                result = score_module(
                    module,
                    absent,
                    present,
                    calls,
                )

                result.update(
                    {
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                    }
                )

                score_rows.append(result)

            processed_records += 1

            if processed_records % 10000 == 0:
                print(
                    f"Processed {processed_records:,} "
                    "biallelic SNPs in plotting interval",
                    flush=True,
                )

    score_fields = [
        "module",
        "chrom",
        "pos",
        "ref",
        "alt",
        "absent_n",
        "present_n",
        "absent_called",
        "present_called",
        "absent_call_rate",
        "present_call_rate",
        "minimum_group_call_rate",
        "selected_orientation",
        "associated_allele",
        "specificity",
        "sensitivity",
        "balanced_accuracy",
        "absent_correct",
        "present_correct",
        "misclassified_n",
        "eligible_min2",
        "callable_80pct",
        "complete_callability",
    ]

    with SCORE_FILE.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=score_fields,
            delimiter="\t",
        )
        writer.writeheader()

        for row in score_rows:
            writer.writerow(
                {
                    key: format_value(row.get(key))
                    for key in score_fields
                }
            )

    focal_fields = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "TN_Left_GT",
        "TN_Left_AD",
        "TN_Left_DP",
        "TN_Left_GQ",
        "TN_Left_DP2_GQ20",
        "TN_Right_GT",
        "TN_Right_AD",
        "TN_Right_DP",
        "TN_Right_GQ",
        "TN_Right_DP2_GQ20",
    ]

    with FOCAL_FILE.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=focal_fields,
            delimiter="\t",
        )
        writer.writeheader()

        for row in focal_rows:
            writer.writerow(
                {
                    key: format_value(row.get(key))
                    for key in focal_fields
                }
            )

    summary_rows = []

    for module in ("Dennis", "Ray"):
        rows = [
            row
            for row in score_rows
            if row["module"] == module
        ]

        eligible = [
            row
            for row in rows
            if row["eligible_min2"]
            and row["balanced_accuracy"] is not None
        ]

        summary_rows.append(
            {
                "module": module,
                "absent_n": rows[0]["absent_n"],
                "present_n": rows[0]["present_n"],
                "biallelic_SNPs": len(rows),
                "eligible_min2": len(eligible),
                "BA_ge_0.80": sum(
                    row["balanced_accuracy"] >= 0.80
                    for row in eligible
                ),
                "BA_ge_0.90": sum(
                    row["balanced_accuracy"] >= 0.90
                    for row in eligible
                ),
                "BA_equal_1_among_called": sum(
                    abs(row["balanced_accuracy"] - 1) < 1e-12
                    for row in eligible
                ),
                "BA_equal_1_80pct_callable": sum(
                    abs(row["balanced_accuracy"] - 1) < 1e-12
                    and row["callable_80pct"]
                    for row in eligible
                ),
                "BA_equal_1_complete": sum(
                    abs(row["balanced_accuracy"] - 1) < 1e-12
                    and row["complete_callability"]
                    for row in eligible
                ),
            }
        )

    with SUMMARY_FILE.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(summary_rows[0]),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    known_positions = (
        set(range(1215000, 1220000))
        | KNOWN_RAY_POSITIONS
    )

    known_rows = [
        row
        for row in score_rows
        if row["pos"] in known_positions
        and row["eligible_min2"]
        and row["balanced_accuracy"] is not None
        and abs(row["balanced_accuracy"] - 1) < 1e-12
    ]

    with KNOWN_FILE.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=score_fields,
            delimiter="\t",
        )
        writer.writeheader()

        for row in known_rows:
            writer.writerow(
                {
                    key: format_value(row.get(key))
                    for key in score_fields
                }
            )

    print()
    print("WALLBANK-COHORT DOMINANT D/R BALANCED ACCURACY")
    print()
    print("SCORING GROUPS")
    print("module\tabsent\tpresent")

    for row in summary_rows:
        print(
            f'{row["module"]}\t'
            f'{row["absent_n"]}\t'
            f'{row["present_n"]}'
        )

    print()
    print("SITE SUMMARY")
    print(
        "module\tSNPs\teligible_min2\tBA80\tBA90\t"
        "BA1_among_called\tBA1_80pct\tBA1_complete"
    )

    for row in summary_rows:
        print(
            f'{row["module"]}\t'
            f'{row["biallelic_SNPs"]}\t'
            f'{row["eligible_min2"]}\t'
            f'{row["BA_ge_0.80"]}\t'
            f'{row["BA_ge_0.90"]}\t'
            f'{row["BA_equal_1_among_called"]}\t'
            f'{row["BA_equal_1_80pct_callable"]}\t'
            f'{row["BA_equal_1_complete"]}'
        )

    print()
    print("BA model: sitewise best dominant allele orientation.")
    print("Absent: homozygous non-pattern state.")
    print("Present: heterozygous or pattern-associated homozygote.")
    print("Comparative QC: DP 4-300 and GQ >=30.")
    print("TN_Left and TN_Right were not used in BA calculation.")
    print()
    print(f"Scores: {SCORE_FILE}")
    print(f"Focal calls: {FOCAL_FILE}")
    print(f"Samples: {SAMPLE_FILE}")
    print(f"Summary: {SUMMARY_FILE}")
    print(f"Known-site audit: {KNOWN_FILE}")


if __name__ == "__main__":
    main()
