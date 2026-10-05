#!/usr/bin/env python3

from pathlib import Path
import csv
import math


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

SCORE_FILE = (
    PROJECT
    / "analysis/original_502/checkpoints/"
    "optix_original_502_comparative_scores.tsv"
)

OUTPUT_ALL = (
    PROJECT
    / "analysis/original_502/checkpoints/"
    "optix_original_502_paired_site_evidence.tsv"
)

OUTPUT_JOINT = (
    PROJECT
    / "analysis/original_502/checkpoints/"
    "optix_original_502_joint_focal_candidates.tsv"
)


def as_bool(value):
    return str(value).strip().lower() in {
        "true", "t", "1", "yes"
    }


def as_float(value):
    if value in {None, "", ".", "NA", "None"}:
        return None

    try:
        number = float(value)
    except ValueError:
        return None

    if not math.isfinite(number):
        return None

    return number


def as_int(value):
    number = as_float(value)

    if number is None:
        return None

    return int(number)


def missingness_ok(value):
    return value is None or value > 0.05


def minimum_defined(first, second):
    if first is None or second is None:
        return None

    return min(first, second)


def load_scores():
    rows = {}

    with SCORE_FILE.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            key = (
                row["panel"],
                row["trait"],
                row["chrom"],
                int(row["pos"]),
                row["ref"],
                row["alt"],
            )

            scheme = row["scheme"]

            parsed = dict(row)
            parsed["pos"] = int(row["pos"])

            for field in [
                "balanced_accuracy",
                "MCC",
                "BH_q",
                "missingness_BH_q",
                "TN_Left_DP",
                "TN_Left_GQ",
                "TN_Left_ALT_fraction",
                "TN_Right_DP",
                "TN_Right_GQ",
                "TN_Right_ALT_fraction",
            ]:
                parsed[field] = as_float(row[field])

            for field in [
                "misclassified_n",
                "phenotype_0_called",
                "phenotype_1_called",
                "predicted_0",
                "predicted_1",
            ]:
                parsed[field] = as_int(row[field])

            for field in [
                "callability_eligible",
                "score_eligible",
                "matches_reference_expectation",
                "pooled_both_DP2_GQ20",
                "pooled_GT_different",
                "left_11_right_01",
                "dominant_focal_model_match",
            ]:
                parsed[field] = as_bool(row[field])

            rows.setdefault(key, {})[scheme] = parsed

    return rows


def main():
    if not SCORE_FILE.is_file():
        print("ERROR: score file is missing:")
        print(SCORE_FILE)
        return

    scores = load_scores()
    paired_rows = []

    for key, schemes in scores.items():
        raw = schemes.get("called_GT")
        qc = schemes.get("DP2_GQ20")

        if raw is None or qc is None:
            continue

        panel, trait, chrom, pos, ref, alt = key

        stable_orientation = (
            raw["selected_orientation"]
            == qc["selected_orientation"]
        )

        both_eligible = (
            raw["score_eligible"]
            and qc["score_eligible"]
        )

        both_q_significant = (
            raw["BH_q"] is not None
            and qc["BH_q"] is not None
            and raw["BH_q"] <= 0.05
            and qc["BH_q"] <= 0.05
        )

        both_missingness_ok = (
            missingness_ok(
                raw["missingness_BH_q"]
            )
            and missingness_ok(
                qc["missingness_BH_q"]
            )
        )

        minimum_ba = minimum_defined(
            raw["balanced_accuracy"],
            qc["balanced_accuracy"],
        )

        minimum_mcc = minimum_defined(
            raw["MCC"],
            qc["MCC"],
        )

        comparative_stable = (
            both_eligible
            and stable_orientation
            and both_q_significant
        )

        comparative_high = (
            comparative_stable
            and minimum_ba is not None
            and minimum_mcc is not None
            and minimum_ba >= 0.80
            and minimum_mcc >= 0.60
        )

        comparative_very_high = (
            comparative_stable
            and minimum_ba is not None
            and minimum_mcc is not None
            and minimum_ba >= 0.90
            and minimum_mcc >= 0.80
        )

        focal_exact = (
            qc["pooled_both_DP2_GQ20"]
            and qc["left_11_right_01"]
        )

        focal_dominant_match = (
            qc["pooled_both_DP2_GQ20"]
            and qc["dominant_focal_model_match"]
        )

        joint_exact = (
            comparative_stable
            and focal_exact
        )

        joint_exact_missingness_ok = (
            joint_exact
            and both_missingness_ok
        )

        joint_dominant = (
            comparative_stable
            and focal_dominant_match
        )

        if focal_exact:
            focal_class = "left_11_right_01"
        elif focal_dominant_match:
            focal_class = "dominant_model_only"
        elif qc["pooled_both_DP2_GQ20"]:
            focal_class = "pooled_called_no_match"
        else:
            focal_class = "pooled_not_high_confidence"

        paired_rows.append(
            {
                "panel": panel,
                "trait": trait,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "selected_orientation": (
                    qc["selected_orientation"]
                ),
                "stable_orientation": (
                    stable_orientation
                ),
                "matches_reference_expectation": (
                    raw[
                        "matches_reference_expectation"
                    ]
                    and qc[
                        "matches_reference_expectation"
                    ]
                ),
                "both_score_eligible": both_eligible,
                "both_q_significant": (
                    both_q_significant
                ),
                "both_missingness_ok": (
                    both_missingness_ok
                ),
                "comparative_stable": (
                    comparative_stable
                ),
                "comparative_high_BA80": (
                    comparative_high
                ),
                "comparative_very_high_BA90": (
                    comparative_very_high
                ),
                "raw_BA": raw["balanced_accuracy"],
                "QC_BA": qc["balanced_accuracy"],
                "minimum_BA": minimum_ba,
                "raw_MCC": raw["MCC"],
                "QC_MCC": qc["MCC"],
                "minimum_MCC": minimum_mcc,
                "raw_errors": raw["misclassified_n"],
                "QC_errors": qc["misclassified_n"],
                "raw_q": raw["BH_q"],
                "QC_q": qc["BH_q"],
                "raw_missingness_q": (
                    raw["missingness_BH_q"]
                ),
                "QC_missingness_q": (
                    qc["missingness_BH_q"]
                ),
                "TN_Left_GT": qc["TN_Left_GT"],
                "TN_Left_AD": qc["TN_Left_AD"],
                "TN_Left_DP": qc["TN_Left_DP"],
                "TN_Left_GQ": qc["TN_Left_GQ"],
                "TN_Left_ALT_fraction": (
                    qc["TN_Left_ALT_fraction"]
                ),
                "TN_Right_GT": qc["TN_Right_GT"],
                "TN_Right_AD": qc["TN_Right_AD"],
                "TN_Right_DP": qc["TN_Right_DP"],
                "TN_Right_GQ": qc["TN_Right_GQ"],
                "TN_Right_ALT_fraction": (
                    qc["TN_Right_ALT_fraction"]
                ),
                "pooled_both_DP2_GQ20": (
                    qc["pooled_both_DP2_GQ20"]
                ),
                "pooled_GT_different": (
                    qc["pooled_GT_different"]
                ),
                "left_11_right_01": (
                    qc["left_11_right_01"]
                ),
                "dominant_focal_model_match": (
                    qc["dominant_focal_model_match"]
                ),
                "focal_class": focal_class,
                "joint_exact": joint_exact,
                "joint_exact_missingness_ok": (
                    joint_exact_missingness_ok
                ),
                "joint_dominant": joint_dominant,
            }
        )

    paired_rows.sort(
        key=lambda row: (
            row["panel"],
            -int(row["joint_exact"]),
            -int(row["both_missingness_ok"]),
            -(
                row["minimum_BA"]
                if row["minimum_BA"] is not None
                else -1
            ),
            -(
                row["minimum_MCC"]
                if row["minimum_MCC"] is not None
                else -2
            ),
            row["pos"],
        )
    )

    joint_rows = [
        row for row in paired_rows
        if row["joint_exact"]
    ]

    OUTPUT_ALL.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fields = list(paired_rows[0].keys())

    with OUTPUT_ALL.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(paired_rows)

    with OUTPUT_JOINT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(joint_rows)

    panel_order = [
        "band",
        "dennis1",
        "dennis2",
        "ray1",
        "ray2",
    ]

    print("ORIGINAL 502-SAMPLE JOINT-EVIDENCE AUDIT")
    print()
    print(
        "panel\tpaired_SNPs\tcomparative_stable\t"
        "BA80\tBA90\tfocal_exact\tjoint_exact\t"
        "joint_exact_missingness_ok"
    )

    for panel in panel_order:
        rows = [
            row for row in paired_rows
            if row["panel"] == panel
        ]

        print(
            f"{panel}\t"
            f"{len(rows)}\t"
            f"{sum(row['comparative_stable'] for row in rows)}\t"
            f"{sum(row['comparative_high_BA80'] for row in rows)}\t"
            f"{sum(row['comparative_very_high_BA90'] for row in rows)}\t"
            f"{sum(row['left_11_right_01'] and row['pooled_both_DP2_GQ20'] for row in rows)}\t"
            f"{sum(row['joint_exact'] for row in rows)}\t"
            f"{sum(row['joint_exact_missingness_ok'] for row in rows)}"
        )

    print()
    print("JOINT EXACT SITES")
    print(
        "panel\tposition\torientation\tmin_BA\tmin_MCC\t"
        "raw_errors\tQC_errors\tmissingness_ok\t"
        "left_GT\tleft_AD\tright_GT\tright_AD"
    )

    for row in joint_rows:
        print(
            f'{row["panel"]}\t'
            f'{row["chrom"]}:{row["pos"]}\t'
            f'{row["selected_orientation"]}\t'
            f'{row["minimum_BA"]:.4f}\t'
            f'{row["minimum_MCC"]:.4f}\t'
            f'{row["raw_errors"]}\t'
            f'{row["QC_errors"]}\t'
            f'{row["both_missingness_ok"]}\t'
            f'{row["TN_Left_GT"]}\t'
            f'{row["TN_Left_AD"]}\t'
            f'{row["TN_Right_GT"]}\t'
            f'{row["TN_Right_AD"]}'
        )

    print()
    print(f"All paired sites: {OUTPUT_ALL}")
    print(f"Joint focal candidates: {OUTPUT_JOINT}")


if __name__ == "__main__":
    main()
