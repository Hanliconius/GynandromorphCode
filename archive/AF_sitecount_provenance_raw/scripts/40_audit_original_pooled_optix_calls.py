#!/usr/bin/env python3

from pathlib import Path
import csv
import re


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

SOURCE_DIR = (
    PROJECT
    / "external/krista_genotype_plots_original/"
    "Genotype_Plots"
)

PRIOR_CANDIDATES = (
    PROJECT
    / "analysis/checkpoints/"
    "optix_51sample_stable_candidate_sites.tsv"
)

SUMMARY_OUTPUT = (
    PROJECT
    / "analysis/checkpoints/"
    "original_502_pooled_optix_summary.tsv"
)

PRIOR_SITE_OUTPUT = (
    PROJECT
    / "analysis/checkpoints/"
    "original_502_calls_at_provisional_51sample_candidates.tsv"
)

VCFS = {
    "band": SOURCE_DIR / "optix_band-1112651-1114314.vcf",
    "dennis1": SOURCE_DIR / "dennis1-1215210-1216124.vcf",
    "dennis2": SOURCE_DIR / "dennis2-1218557-1221494.vcf",
    "ray1": SOURCE_DIR / "ray1-1156746-1160220.vcf",
    "ray2": SOURCE_DIR / "ray2-1174547-1187833.vcf",
}

POOLED = ["TN_Left", "TN_Right"]


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

    ad = parse_integer_vector(fields.get("AD"))

    alt_fraction = None

    if ad is not None and len(ad) == 2:
        total = ad[0] + ad[1]

        if total > 0:
            alt_fraction = ad[1] / total

    return {
        "GT": canonical_gt(fields.get("GT")),
        "AD": ad,
        "DP": parse_number(fields.get("DP")),
        "GQ": parse_number(fields.get("GQ")),
        "alt_fraction": alt_fraction,
    }


def qc_pass(call):
    return (
        call["GT"] is not None
        and call["DP"] is not None
        and call["DP"] >= 2
        and call["GQ"] is not None
        and call["GQ"] >= 20
    )


def vector_text(value):
    if value is None:
        return ""
    return ",".join(str(part) for part in value)


def number_text(value, digits=5):
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def load_prior_candidates():
    rows = []

    with PRIOR_CANDIDATES.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            if str(row["tier_1"]).lower() == "true":
                tier = "tier_1"
            elif str(row["tier_2"]).lower() == "true":
                tier = "tier_2"
            else:
                continue

            rows.append(
                {
                    "panel": row["panel"],
                    "trait": row["trait"],
                    "chrom": row["chrom"],
                    "pos": int(row["pos"]),
                    "ref": row["ref"],
                    "alt": row["alt"],
                    "tier": tier,
                    "minimum_BA_51sample": row["minimum_BA"],
                    "minimum_MCC_51sample": row["minimum_MCC"],
                }
            )

    return rows


def main():
    required = [*VCFS.values(), PRIOR_CANDIDATES]
    missing = [path for path in required if not path.is_file()]

    if missing:
        print("ERROR: required files are missing:")
        for path in missing:
            print(path)
        return

    prior_candidates = load_prior_candidates()

    summary_rows = []
    prior_rows = []

    reference_sample_order = None

    for panel, path in VCFS.items():
        samples = []
        pooled_indices = {}
        records = 0
        biallelic_snps = 0
        both_called = 0
        both_qc = 0
        different_gt = 0
        different_gt_qc = 0
        left_11_right_01 = 0
        left_11_right_01_qc = 0

        records_by_exact_site = {}
        records_by_position = {}

        with path.open(errors="replace") as handle:
            for line in handle:
                if line.startswith("##"):
                    continue

                if line.startswith("#CHROM"):
                    fields = line.rstrip("\n").split("\t")
                    samples = fields[9:]

                    pooled_indices = {
                        sample: samples.index(sample)
                        for sample in POOLED
                        if sample in samples
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

                records += 1

                if not (
                    "," not in alt
                    and len(ref) == 1
                    and len(alt) == 1
                ):
                    continue

                biallelic_snps += 1

                format_ids = fields[8].split(":")
                sample_values = fields[9:]

                calls = {}

                for sample in POOLED:
                    if sample not in pooled_indices:
                        calls[sample] = {
                            "GT": None,
                            "AD": None,
                            "DP": None,
                            "GQ": None,
                            "alt_fraction": None,
                        }
                    else:
                        calls[sample] = parse_call(
                            format_ids,
                            sample_values[
                                pooled_indices[sample]
                            ],
                        )

                left = calls["TN_Left"]
                right = calls["TN_Right"]

                if (
                    left["GT"] is not None
                    and right["GT"] is not None
                ):
                    both_called += 1

                    if left["GT"] != right["GT"]:
                        different_gt += 1

                    if (
                        left["GT"] == "1/1"
                        and right["GT"] == "0/1"
                    ):
                        left_11_right_01 += 1

                if qc_pass(left) and qc_pass(right):
                    both_qc += 1

                    if left["GT"] != right["GT"]:
                        different_gt_qc += 1

                    if (
                        left["GT"] == "1/1"
                        and right["GT"] == "0/1"
                    ):
                        left_11_right_01_qc += 1

                site_data = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "left": left,
                    "right": right,
                }

                records_by_exact_site[
                    (chrom, pos, ref, alt)
                ] = site_data

                records_by_position.setdefault(
                    (chrom, pos),
                    [],
                ).append(site_data)

        if reference_sample_order is None:
            reference_sample_order = samples

        sample_order_identical = (
            samples == reference_sample_order
        )

        summary_rows.append(
            {
                "panel": panel,
                "samples": len(samples),
                "sample_order_identical": (
                    sample_order_identical
                ),
                "records": records,
                "biallelic_SNPs": biallelic_snps,
                "both_pooled_called": both_called,
                "both_pooled_DP2_GQ20": both_qc,
                "different_pooled_GT_raw": different_gt,
                "different_pooled_GT_DP2_GQ20": (
                    different_gt_qc
                ),
                "left_11_right_01_raw": (
                    left_11_right_01
                ),
                "left_11_right_01_DP2_GQ20": (
                    left_11_right_01_qc
                ),
            }
        )

        panel_prior = [
            row for row in prior_candidates
            if row["panel"] == panel
        ]

        for candidate in panel_prior:
            exact_key = (
                candidate["chrom"],
                candidate["pos"],
                candidate["ref"],
                candidate["alt"],
            )

            position_key = (
                candidate["chrom"],
                candidate["pos"],
            )

            site = records_by_exact_site.get(exact_key)
            exact_match = site is not None
            position_match = (
                position_key in records_by_position
            )

            if site is None and position_match:
                site = records_by_position[
                    position_key
                ][0]

            left = (
                site["left"]
                if site is not None
                else {}
            )

            right = (
                site["right"]
                if site is not None
                else {}
            )

            prior_rows.append(
                {
                    **candidate,
                    "exact_REF_ALT_match": exact_match,
                    "position_match": position_match,
                    "original_ref": (
                        site["ref"]
                        if site is not None
                        else ""
                    ),
                    "original_alt": (
                        site["alt"]
                        if site is not None
                        else ""
                    ),
                    "TN_Left_GT": left.get("GT") or "",
                    "TN_Left_AD": vector_text(
                        left.get("AD")
                    ),
                    "TN_Left_DP": number_text(
                        left.get("DP")
                    ),
                    "TN_Left_GQ": number_text(
                        left.get("GQ")
                    ),
                    "TN_Left_ALT_fraction": number_text(
                        left.get("alt_fraction")
                    ),
                    "TN_Right_GT": right.get("GT") or "",
                    "TN_Right_AD": vector_text(
                        right.get("AD")
                    ),
                    "TN_Right_DP": number_text(
                        right.get("DP")
                    ),
                    "TN_Right_GQ": number_text(
                        right.get("GQ")
                    ),
                    "TN_Right_ALT_fraction": number_text(
                        right.get("alt_fraction")
                    ),
                    "pooled_GT_different": (
                        site is not None
                        and left.get("GT") is not None
                        and right.get("GT") is not None
                        and left.get("GT") != right.get("GT")
                    ),
                    "left_11_right_01": (
                        site is not None
                        and left.get("GT") == "1/1"
                        and right.get("GT") == "0/1"
                    ),
                    "both_DP2_GQ20": (
                        site is not None
                        and qc_pass(left)
                        and qc_pass(right)
                    ),
                }
            )

    SUMMARY_OUTPUT.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with SUMMARY_OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(summary_rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    with PRIOR_SITE_OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(prior_rows[0].keys()),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(prior_rows)

    print("ORIGINAL 502-SAMPLE POOLED OPTIX AUDIT")
    print()
    print(
        "panel\tSNPs\tboth_called\tboth_DP2GQ20\t"
        "different_GT_QC\tL11_R01_QC\t"
        "prior_candidates\tprior_exact_overlap"
    )

    for summary in summary_rows:
        panel = summary["panel"]

        panel_prior = [
            row for row in prior_rows
            if row["panel"] == panel
        ]

        exact_overlap = sum(
            row["exact_REF_ALT_match"]
            for row in panel_prior
        )

        print(
            f'{panel}\t'
            f'{summary["biallelic_SNPs"]}\t'
            f'{summary["both_pooled_called"]}\t'
            f'{summary["both_pooled_DP2_GQ20"]}\t'
            f'{summary["different_pooled_GT_DP2_GQ20"]}\t'
            f'{summary["left_11_right_01_DP2_GQ20"]}\t'
            f'{len(panel_prior)}\t'
            f'{exact_overlap}'
        )

    print()
    print("OVERLAPPING PROVISIONAL CANDIDATES")
    print(
        "panel\tposition\ttier\tleft_GT\tleft_AD\t"
        "right_GT\tright_AD\tL11_R01\tQC"
    )

    for row in prior_rows:
        if not row["position_match"]:
            continue

        print(
            f'{row["panel"]}\t'
            f'{row["chrom"]}:{row["pos"]}\t'
            f'{row["tier"]}\t'
            f'{row["TN_Left_GT"]}\t'
            f'{row["TN_Left_AD"]}\t'
            f'{row["TN_Right_GT"]}\t'
            f'{row["TN_Right_AD"]}\t'
            f'{row["left_11_right_01"]}\t'
            f'{row["both_DP2_GQ20"]}'
        )

    print()
    print(f"Summary: {SUMMARY_OUTPUT}")
    print(f"Prior-site audit: {PRIOR_SITE_OUTPUT}")
    print()
    print(
        "The provisional 51-sample and pseudo-pool files "
        "were not modified."
    )


if __name__ == "__main__":
    main()
