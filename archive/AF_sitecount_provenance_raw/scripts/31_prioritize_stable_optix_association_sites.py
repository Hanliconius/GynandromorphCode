#!/usr/bin/env python3

import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

from openpyxl import load_workbook

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

SCORES = (
    PROJECT
    / "analysis"
    / "checkpoints"
    / "optix_51sample_comparative_scores.tsv"
)

WORKBOOK = (
    PROJECT
    / "metadata"
    / "Wing_Pattern_MetaData.xlsx"
)

BAND_VCF = (
    PROJECT
    / "data"
    / "derived"
    / "reconstructed_plot_inputs"
    / "optix_band-1112651-1114314.vcf"
)

OUTPUT = (
    PROJECT
    / "analysis"
    / "checkpoints"
    / "optix_51sample_stable_candidate_sites.tsv"
)

EXCLUDED_OUTPUT = (
    PROJECT
    / "analysis"
    / "checkpoints"
    / "optix_51sample_excluded_comparative_samples.tsv"
)

FOCAL = {
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
}


def as_float(value):
    if value in {"", "None", None}:
        return None

    return float(value)


def as_bool(value):
    return str(value).lower() == "true"


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

    index = {
        header: position
        for position, header in enumerate(headers)
    }

    metadata = {}

    for values in rows:
        accession_value = values[index["Accession"]]

        if accession_value is None:
            continue

        accession = str(accession_value).strip()

        if not accession:
            continue

        metadata[accession] = {
            "Taxon Name": (
                str(values[index["Taxon Name"]]).strip()
                if values[index["Taxon Name"]] is not None
                else ""
            ),
            "B": parse_binary(values[index["B"]]),
            "R": parse_binary(values[index["R"]]),
            "D": parse_binary(values[index["D"]]),
        }

    return metadata


def get_samples():
    result = subprocess.run(
        ["bcftools", "query", "-l", str(BAND_VCF)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        return []

    return [
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip()
    ]


def missingness_ok(row):
    q_value = row.get("missingness_BH_q")

    return q_value is None or q_value > 0.05


def main():
    missing = [
        str(path)
        for path in [SCORES, WORKBOOK, BAND_VCF]
        if not path.is_file()
    ]

    if missing:
        print("ERROR: required inputs are missing:")
        for path in missing:
            print(path)
        return

    with SCORES.open(
        "r",
        encoding="utf-8",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    numeric_columns = [
        "call_rate",
        "phenotype_0_call_rate",
        "phenotype_1_call_rate",
        "call_rate_state_difference",
        "missingness_fisher_p",
        "sensitivity",
        "specificity",
        "balanced_accuracy",
        "MCC",
        "fisher_p",
        "BH_q",
        "missingness_BH_q",
    ]

    integer_columns = [
        "phenotype_0_called",
        "phenotype_1_called",
        "predicted_0",
        "predicted_1",
        "misclassified_n",
    ]

    boolean_columns = [
        "callability_eligible",
        "score_eligible",
        "near_perfect",
        "matches_reference_expectation",
    ]

    for row in rows:
        for column in numeric_columns:
            row[column] = as_float(row.get(column))

        for column in integer_columns:
            row[column] = int(row[column])

        for column in boolean_columns:
            row[column] = as_bool(row[column])

    paired = defaultdict(dict)

    for row in rows:
        key = (
            row["panel"],
            row["trait"],
            row["chrom"],
            int(row["pos"]),
            row["ref"],
            row["alt"],
        )

        paired[key][row["scheme"]] = row

    candidate_rows = []

    for key, schemes in paired.items():
        raw = schemes.get("called_GT")
        qc = schemes.get("DP2_GQ20")

        if raw is None or qc is None:
            continue

        stable_orientation = (
            raw["selected_orientation"]
            == qc["selected_orientation"]
        )

        both_score_eligible = (
            raw["score_eligible"]
            and qc["score_eligible"]
        )

        both_missingness_ok = (
            missingness_ok(raw)
            and missingness_ok(qc)
        )

        both_q_significant = (
            raw["BH_q"] is not None
            and qc["BH_q"] is not None
            and raw["BH_q"] <= 0.05
            and qc["BH_q"] <= 0.05
        )

        ba_values = [
            raw["balanced_accuracy"],
            qc["balanced_accuracy"],
        ]
        mcc_values = [
            raw["MCC"],
            qc["MCC"],
        ]

        minimum_ba = (
            min(ba_values)
            if all(value is not None for value in ba_values)
            else None
        )

        minimum_mcc = (
            min(mcc_values)
            if all(value is not None for value in mcc_values)
            else None
        )

        tier_1 = (
            both_score_eligible
            and stable_orientation
            and both_missingness_ok
            and both_q_significant
            and minimum_ba is not None
            and minimum_mcc is not None
            and minimum_ba >= 0.90
            and minimum_mcc >= 0.80
        )

        tier_2 = (
            not tier_1
            and both_score_eligible
            and stable_orientation
            and both_missingness_ok
            and both_q_significant
            and minimum_ba is not None
            and minimum_mcc is not None
            and minimum_ba >= 0.80
            and minimum_mcc >= 0.60
        )

        panel, trait, chrom, pos, ref, alt = key

        candidate_rows.append(
            {
                "panel": panel,
                "trait": trait,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "orientation": (
                    qc["selected_orientation"]
                ),
                "stable_orientation": stable_orientation,
                "raw_score_eligible": (
                    raw["score_eligible"]
                ),
                "QC_score_eligible": (
                    qc["score_eligible"]
                ),
                "raw_called_n": (
                    raw["phenotype_0_called"]
                    + raw["phenotype_1_called"]
                ),
                "QC_called_n": (
                    qc["phenotype_0_called"]
                    + qc["phenotype_1_called"]
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
                "both_missingness_ok": (
                    both_missingness_ok
                ),
                "both_q_significant": (
                    both_q_significant
                ),
                "tier_1": tier_1,
                "tier_2": tier_2,
            }
        )

    candidate_rows.sort(
        key=lambda row: (
            row["panel"],
            not row["tier_1"],
            not row["tier_2"],
            -(row["minimum_BA"] if row["minimum_BA"] is not None else -1.0),
            -(row["minimum_MCC"] if row["minimum_MCC"] is not None else -2.0),
            row["pos"],
        )
    )

    fieldnames = list(candidate_rows[0].keys())

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
        writer.writerows(candidate_rows)

    metadata = load_metadata()
    samples = get_samples()

    excluded_rows = []

    for sample in samples:
        if sample in FOCAL:
            continue

        sample_metadata = metadata.get(sample)

        if sample_metadata is None:
            reason = "no_metadata_match"
            values = {
                "Taxon Name": "",
                "B": None,
                "R": None,
                "D": None,
            }
        else:
            values = sample_metadata
            missing_traits = [
                trait
                for trait in ["B", "R", "D"]
                if values[trait] is None
            ]

            if missing_traits:
                reason = (
                    "missing_phenotype:"
                    + ",".join(missing_traits)
                )
            else:
                continue

        excluded_rows.append(
            {
                "sample": sample,
                "reason": reason,
                "Taxon Name": values["Taxon Name"],
                "B": values["B"],
                "R": values["R"],
                "D": values["D"],
            }
        )

    with EXCLUDED_OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        excluded_fields = [
            "sample",
            "reason",
            "Taxon Name",
            "B",
            "R",
            "D",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=excluded_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(excluded_rows)

    print("STABLE OPTIX CANDIDATE PRIORITIZATION")

    print("\nCANDIDATE COUNTS")
    print(
        "panel\tpaired_SNPs\tboth_eligible\t"
        "tier_1\ttier_2\tmissingness_flagged"
    )

    panels = sorted({
        row["panel"]
        for row in candidate_rows
    })

    for panel in panels:
        panel_rows = [
            row
            for row in candidate_rows
            if row["panel"] == panel
        ]

        print(
            f"{panel}\t"
            f"{len(panel_rows)}\t"
            f"{sum(row['raw_score_eligible'] and row['QC_score_eligible'] for row in panel_rows)}\t"
            f"{sum(row['tier_1'] for row in panel_rows)}\t"
            f"{sum(row['tier_2'] for row in panel_rows)}\t"
            f"{sum(not row['both_missingness_ok'] for row in panel_rows)}"
        )

    print("\nTIER-1 AND TIER-2 SITES")
    print(
        "panel\ttrait\tposition\torientation\t"
        "raw_BA\tQC_BA\traw_MCC\tQC_MCC\t"
        "raw_errors\tQC_errors\ttier"
    )

    prioritized = [
        row
        for row in candidate_rows
        if row["tier_1"] or row["tier_2"]
    ]

    for row in prioritized:
        tier = "tier_1" if row["tier_1"] else "tier_2"

        print(
            f"{row['panel']}\t"
            f"{row['trait']}\t"
            f"{row['chrom']}:{row['pos']}\t"
            f"{row['orientation']}\t"
            f"{row['raw_BA']:.4f}\t"
            f"{row['QC_BA']:.4f}\t"
            f"{row['raw_MCC']:.4f}\t"
            f"{row['QC_MCC']:.4f}\t"
            f"{row['raw_errors']}\t"
            f"{row['QC_errors']}\t"
            f"{tier}"
        )

    print("\nEXCLUDED NONFOCAL SAMPLES")
    print("sample\treason\tTaxon Name\tB\tR\tD")

    for row in excluded_rows:
        print(
            f"{row['sample']}\t"
            f"{row['reason']}\t"
            f"{row['Taxon Name']}\t"
            f"{row['B']}\t"
            f"{row['R']}\t"
            f"{row['D']}"
        )

    print(f"\nCandidate table: {OUTPUT}")
    print(f"Excluded samples: {EXCLUDED_OUTPUT}")


if __name__ == "__main__":
    main()
