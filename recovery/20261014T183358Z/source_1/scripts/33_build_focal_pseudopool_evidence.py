#!/usr/bin/env python3

from pathlib import Path
import csv
import re


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

CANDIDATES = (
    PROJECT
    / "analysis/checkpoints/optix_51sample_stable_candidate_sites.tsv"
)

PATTERN_GROUPS = (
    PROJECT
    / "analysis/checkpoints/optix_51sample_candidate_pattern_groups.tsv"
)

OUTPUT_LONG = (
    PROJECT
    / "analysis/checkpoints/"
    "optix_focal_pseudopool_evidence.long.tsv"
)

OUTPUT_WIDE = (
    PROJECT
    / "analysis/checkpoints/"
    "optix_focal_pseudopool_evidence.wide.tsv"
)

VCF_BY_PANEL = {
    "band": (
        PROJECT
        / "data/derived/reconstructed_plot_inputs/"
        "optix_band-1112651-1114314.vcf"
    ),
    "dennis1": (
        PROJECT
        / "data/derived/reconstructed_plot_inputs/"
        "dennis1-1215210-1216124.vcf"
    ),
    "dennis2": (
        PROJECT
        / "data/derived/reconstructed_plot_inputs/"
        "dennis2-1218557-1221494.vcf"
    ),
    "ray1": (
        PROJECT
        / "data/derived/reconstructed_plot_inputs/"
        "ray1-1156746-1160220.vcf"
    ),
    "ray2": (
        PROJECT
        / "data/derived/reconstructed_plot_inputs/"
        "ray2-1174547-1187833.vcf"
    ),
}

SIDE_SAMPLES = {
    "left": ["TN2L_S1", "TN3L_S3"],
    "right": ["TN2R_S2", "TN3R_S4"],
}

GT_FROM_PL_INDEX = {
    0: "0/0",
    1: "0/1",
    2: "1/1",
}

MIN_COMBINED_AD_DEPTH = 10
MIN_PSEUDOPOOL_GQ = 20
BALANCED_LOWER = 0.30
BALANCED_UPPER = 0.70
PRIMARY_ALT_ENRICHED_THRESHOLD = 0.80


def as_bool(value):
    return str(value).strip().lower() in {
        "true", "t", "1", "yes"
    }


def as_float(value):
    if value in {None, "", ".", "NA", "None"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def as_int(value):
    number = as_float(value)
    if number is None:
        return None
    return int(number)


def canonical_gt(gt):
    if gt in {None, "", ".", "./.", ".|."}:
        return None

    alleles = re.split(r"[/|]", gt)

    if any(allele in {"", "."} for allele in alleles):
        return None

    try:
        numeric = sorted(int(allele) for allele in alleles)
    except ValueError:
        return None

    return "/".join(str(value) for value in numeric)


def parse_integer_vector(value):
    if value in {None, "", "."}:
        return None

    try:
        return [int(part) for part in value.split(",")]
    except ValueError:
        return None


def parse_call(format_ids, sample_field):
    values = sample_field.split(":")
    fields = dict(zip(format_ids, values))

    return {
        "GT": canonical_gt(fields.get("GT")),
        "AD": parse_integer_vector(fields.get("AD")),
        "DP": as_int(fields.get("DP")),
        "GQ": as_float(fields.get("GQ")),
        "PL": parse_integer_vector(fields.get("PL")),
    }


def load_candidate_rows():
    rows = []

    with CANDIDATES.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            if not (
                as_bool(row["tier_1"])
                or as_bool(row["tier_2"])
            ):
                continue

            row["pos"] = int(row["pos"])
            row["tier"] = (
                "tier_1"
                if as_bool(row["tier_1"])
                else "tier_2"
            )

            for field in [
                "minimum_BA",
                "minimum_MCC",
                "raw_BA",
                "QC_BA",
                "raw_MCC",
                "QC_MCC",
            ]:
                row[field] = as_float(row[field])

            for field in [
                "raw_errors",
                "QC_errors",
                "raw_called_n",
                "QC_called_n",
            ]:
                row[field] = as_int(row[field])

            rows.append(row)

    return rows


def load_group_assignments():
    assignments = {}

    with PATTERN_GROUPS.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            key = (
                row["panel"],
                row["chrom"],
                int(row["pos"]),
                row["ref"],
                row["alt"],
            )

            assignments[key] = {
                "group_id": row["group_id"],
                "is_representative": as_bool(
                    row["is_representative"]
                ),
            }

    return assignments


def load_vcf_records(vcf_path, wanted_keys):
    samples = []
    records = {}

    with vcf_path.open() as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                continue

            key = (
                fields[0],
                int(fields[1]),
                fields[3],
                fields[4],
            )

            if key not in wanted_keys:
                continue

            format_ids = fields[8].split(":")
            calls = {}

            for sample, sample_field in zip(
                samples, fields[9:]
            ):
                calls[sample] = parse_call(
                    format_ids,
                    sample_field,
                )

            records[key] = calls

    return samples, records


def combine_side(calls, samples):
    leg_calls = [calls[sample] for sample in samples]

    valid_ad = [
        call["AD"]
        for call in leg_calls
        if call["AD"] is not None
        and len(call["AD"]) == 2
    ]

    if valid_ad:
        ref_depth = sum(ad[0] for ad in valid_ad)
        alt_depth = sum(ad[1] for ad in valid_ad)
        ad_depth = ref_depth + alt_depth
        alt_fraction = (
            alt_depth / ad_depth
            if ad_depth > 0
            else None
        )
    else:
        ref_depth = None
        alt_depth = None
        ad_depth = None
        alt_fraction = None

    valid_dp = [
        call["DP"]
        for call in leg_calls
        if call["DP"] is not None
    ]

    combined_dp = (
        sum(valid_dp)
        if valid_dp
        else None
    )

    valid_pl = [
        call["PL"]
        for call in leg_calls
        if call["PL"] is not None
        and len(call["PL"]) == 3
    ]

    combined_pl = None
    normalized_pl = None
    pseudo_gt = None
    pseudo_gq = None

    if len(valid_pl) == 2:
        combined_pl = [
            valid_pl[0][index] + valid_pl[1][index]
            for index in range(3)
        ]

        minimum = min(combined_pl)

        normalized_pl = [
            value - minimum
            for value in combined_pl
        ]

        ordered = sorted(
            enumerate(normalized_pl),
            key=lambda item: (item[1], item[0]),
        )

        best_index, best_value = ordered[0]
        second_value = ordered[1][1]

        pseudo_gt = GT_FROM_PL_INDEX[best_index]
        pseudo_gq = second_value - best_value

    high_confidence = (
        len(valid_pl) == 2
        and len(valid_ad) == 2
        and ad_depth is not None
        and ad_depth >= MIN_COMBINED_AD_DEPTH
        and pseudo_gq is not None
        and pseudo_gq >= MIN_PSEUDOPOOL_GQ
    )

    return {
        "leg1_GT": leg_calls[0]["GT"],
        "leg2_GT": leg_calls[1]["GT"],
        "leg1_AD": leg_calls[0]["AD"],
        "leg2_AD": leg_calls[1]["AD"],
        "leg1_DP": leg_calls[0]["DP"],
        "leg2_DP": leg_calls[1]["DP"],
        "leg1_GQ": leg_calls[0]["GQ"],
        "leg2_GQ": leg_calls[1]["GQ"],
        "n_legs_with_AD": len(valid_ad),
        "n_legs_with_PL": len(valid_pl),
        "combined_ref_depth": ref_depth,
        "combined_alt_depth": alt_depth,
        "combined_AD_depth": ad_depth,
        "combined_DP": combined_dp,
        "combined_alt_fraction": alt_fraction,
        "combined_PL": combined_pl,
        "normalized_PL": normalized_pl,
        "pseudo_GT": pseudo_gt,
        "pseudo_GQ": pseudo_gq,
        "high_confidence_pseudopool": high_confidence,
    }


def dominant_prediction(gt, orientation):
    if gt is None:
        return None

    alleles = [int(value) for value in gt.split("/")]

    if orientation == "ALT_carrier_predicts_presence":
        return 1 if any(value > 0 for value in alleles) else 0

    if orientation == "REF_carrier_predicts_presence":
        return 1 if any(value == 0 for value in alleles) else 0

    return None


def expected_dominant_state(trait, side):
    if trait == "B":
        return 0 if side == "left" else 1

    if trait in {"D", "R"}:
        return 1

    return None


def allele_balance_class(alt_fraction):
    if alt_fraction is None:
        return "missing"

    if alt_fraction <= 0.20:
        return "REF_enriched"

    if (
        BALANCED_LOWER
        <= alt_fraction
        <= BALANCED_UPPER
    ):
        return "balanced"

    if alt_fraction >= PRIMARY_ALT_ENRICHED_THRESHOLD:
        return "ALT_enriched"

    return "intermediate"


def vector_text(values):
    if values is None:
        return ""
    return ",".join(str(value) for value in values)


def number_text(value, digits=5):
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def main():
    required = [
        CANDIDATES,
        PATTERN_GROUPS,
        *VCF_BY_PANEL.values(),
    ]

    missing = [path for path in required if not path.is_file()]

    if missing:
        print("ERROR: required inputs are missing:")
        for path in missing:
            print(path)
        return

    candidates = load_candidate_rows()
    groups = load_group_assignments()

    long_rows = []
    wide_rows = []

    for panel, vcf_path in VCF_BY_PANEL.items():
        panel_candidates = [
            row for row in candidates
            if row["panel"] == panel
        ]

        if not panel_candidates:
            continue

        wanted_keys = {
            (
                row["chrom"],
                row["pos"],
                row["ref"],
                row["alt"],
            )
            for row in panel_candidates
        }

        samples, records = load_vcf_records(
            vcf_path,
            wanted_keys,
        )

        required_samples = {
            sample
            for side_samples in SIDE_SAMPLES.values()
            for sample in side_samples
        }

        absent_samples = sorted(required_samples - set(samples))

        if absent_samples:
            print(
                f"ERROR: {panel} lacks focal samples: "
                + ",".join(absent_samples)
            )
            continue

        for candidate in panel_candidates:
            site_key = (
                candidate["chrom"],
                candidate["pos"],
                candidate["ref"],
                candidate["alt"],
            )

            calls = records.get(site_key)

            if calls is None:
                print(
                    "WARNING: candidate absent from VCF: "
                    f"{panel} {site_key}"
                )
                continue

            group_key = (panel, *site_key)
            group = groups.get(
                group_key,
                {
                    "group_id": "",
                    "is_representative": False,
                },
            )

            side_results = {
                side: combine_side(calls, side_samples)
                for side, side_samples
                in SIDE_SAMPLES.items()
            }

            for side, result in side_results.items():
                predicted = dominant_prediction(
                    result["pseudo_GT"],
                    candidate["orientation"],
                )

                expected = expected_dominant_state(
                    candidate["trait"],
                    side,
                )

                long_rows.append(
                    {
                        "panel": panel,
                        "trait": candidate["trait"],
                        "chrom": candidate["chrom"],
                        "pos": candidate["pos"],
                        "ref": candidate["ref"],
                        "alt": candidate["alt"],
                        "orientation": (
                            candidate["orientation"]
                        ),
                        "tier": candidate["tier"],
                        "group_id": group["group_id"],
                        "is_representative": (
                            group["is_representative"]
                        ),
                        "minimum_BA": number_text(
                            candidate["minimum_BA"]
                        ),
                        "minimum_MCC": number_text(
                            candidate["minimum_MCC"]
                        ),
                        "side": side,
                        "leg1_sample": SIDE_SAMPLES[side][0],
                        "leg2_sample": SIDE_SAMPLES[side][1],
                        "leg1_GT": result["leg1_GT"] or "",
                        "leg2_GT": result["leg2_GT"] or "",
                        "leg1_AD": vector_text(
                            result["leg1_AD"]
                        ),
                        "leg2_AD": vector_text(
                            result["leg2_AD"]
                        ),
                        "leg1_DP": number_text(
                            result["leg1_DP"]
                        ),
                        "leg2_DP": number_text(
                            result["leg2_DP"]
                        ),
                        "leg1_GQ": number_text(
                            result["leg1_GQ"]
                        ),
                        "leg2_GQ": number_text(
                            result["leg2_GQ"]
                        ),
                        "n_legs_with_AD": (
                            result["n_legs_with_AD"]
                        ),
                        "n_legs_with_PL": (
                            result["n_legs_with_PL"]
                        ),
                        "combined_ref_depth": number_text(
                            result["combined_ref_depth"]
                        ),
                        "combined_alt_depth": number_text(
                            result["combined_alt_depth"]
                        ),
                        "combined_AD_depth": number_text(
                            result["combined_AD_depth"]
                        ),
                        "combined_DP": number_text(
                            result["combined_DP"]
                        ),
                        "combined_alt_fraction": number_text(
                            result["combined_alt_fraction"]
                        ),
                        "combined_PL": vector_text(
                            result["combined_PL"]
                        ),
                        "normalized_PL": vector_text(
                            result["normalized_PL"]
                        ),
                        "pseudo_GT": (
                            result["pseudo_GT"] or ""
                        ),
                        "pseudo_GQ": number_text(
                            result["pseudo_GQ"]
                        ),
                        "high_confidence_pseudopool": (
                            result[
                                "high_confidence_pseudopool"
                            ]
                        ),
                        "dominant_predicted_state": (
                            ""
                            if predicted is None
                            else predicted
                        ),
                        "dominant_expected_state": expected,
                        "dominant_state_matches": (
                            predicted == expected
                            if predicted is not None
                            else False
                        ),
                        "allele_balance_class": (
                            allele_balance_class(
                                result[
                                    "combined_alt_fraction"
                                ]
                            )
                        ),
                    }
                )

            left = side_results["left"]
            right = side_results["right"]

            left_af = left["combined_alt_fraction"]
            right_af = right["combined_alt_fraction"]

            left_minus_right = (
                left_af - right_af
                if left_af is not None
                and right_af is not None
                else None
            )

            exact_haplotype_match = (
                left["high_confidence_pseudopool"]
                and right["high_confidence_pseudopool"]
                and left["pseudo_GT"] == "1/1"
                and right["pseudo_GT"] == "0/1"
            )

            dominant_left = dominant_prediction(
                left["pseudo_GT"],
                candidate["orientation"],
            )

            dominant_right = dominant_prediction(
                right["pseudo_GT"],
                candidate["orientation"],
            )

            dominant_model_match = (
                dominant_left
                == expected_dominant_state(
                    candidate["trait"],
                    "left",
                )
                and dominant_right
                == expected_dominant_state(
                    candidate["trait"],
                    "right",
                )
            )

            def contamination_aware_match(threshold):
                return (
                    left["n_legs_with_AD"] == 2
                    and right["n_legs_with_AD"] == 2
                    and left["combined_AD_depth"] is not None
                    and right["combined_AD_depth"] is not None
                    and left["combined_AD_depth"]
                    >= MIN_COMBINED_AD_DEPTH
                    and right["combined_AD_depth"]
                    >= MIN_COMBINED_AD_DEPTH
                    and left_af is not None
                    and right_af is not None
                    and left_af >= threshold
                    and BALANCED_LOWER
                    <= right_af
                    <= BALANCED_UPPER
                )

            wide_rows.append(
                {
                    "panel": panel,
                    "trait": candidate["trait"],
                    "chrom": candidate["chrom"],
                    "pos": candidate["pos"],
                    "ref": candidate["ref"],
                    "alt": candidate["alt"],
                    "orientation": candidate["orientation"],
                    "tier": candidate["tier"],
                    "group_id": group["group_id"],
                    "is_representative": (
                        group["is_representative"]
                    ),
                    "minimum_BA": number_text(
                        candidate["minimum_BA"]
                    ),
                    "minimum_MCC": number_text(
                        candidate["minimum_MCC"]
                    ),
                    "raw_errors": candidate["raw_errors"],
                    "QC_errors": candidate["QC_errors"],
                    "left_pseudo_GT": (
                        left["pseudo_GT"] or ""
                    ),
                    "left_pseudo_GQ": number_text(
                        left["pseudo_GQ"]
                    ),
                    "left_ref_depth": number_text(
                        left["combined_ref_depth"]
                    ),
                    "left_alt_depth": number_text(
                        left["combined_alt_depth"]
                    ),
                    "left_alt_fraction": number_text(left_af),
                    "left_balance_class": (
                        allele_balance_class(left_af)
                    ),
                    "right_pseudo_GT": (
                        right["pseudo_GT"] or ""
                    ),
                    "right_pseudo_GQ": number_text(
                        right["pseudo_GQ"]
                    ),
                    "right_ref_depth": number_text(
                        right["combined_ref_depth"]
                    ),
                    "right_alt_depth": number_text(
                        right["combined_alt_depth"]
                    ),
                    "right_alt_fraction": number_text(right_af),
                    "right_balance_class": (
                        allele_balance_class(right_af)
                    ),
                    "left_minus_right_alt_fraction": (
                        number_text(left_minus_right)
                    ),
                    "dominant_model_match": (
                        dominant_model_match
                    ),
                    "exact_left_11_right_01_match": (
                        exact_haplotype_match
                    ),
                    "AB_match_threshold_0.75": (
                        contamination_aware_match(0.75)
                    ),
                    "AB_match_threshold_0.80": (
                        contamination_aware_match(0.80)
                    ),
                    "AB_match_threshold_0.85": (
                        contamination_aware_match(0.85)
                    ),
                    "evidence_scope": (
                        "temporary_VCF_PL_AD_pseudopool;"
                        "not_raw_read_pool"
                    ),
                }
            )

    OUTPUT_LONG.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    long_fields = list(long_rows[0].keys())
    wide_fields = list(wide_rows[0].keys())

    with OUTPUT_LONG.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=long_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(long_rows)

    with OUTPUT_WIDE.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=wide_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(wide_rows)

    print("FOCAL LEG PSEUDO-POOL AUDIT")
    print()
    print(
        "Formal pseudo-genotypes require PL from both legs. "
        "Allele-balance evidence sums AD across both legs."
    )
    print(
        "Contamination-aware pattern: left ALT-enriched; "
        "right ALT fraction 0.30-0.70."
    )
    print()

    print(
        "panel\tstable_sites\trepresentatives\t"
        "formal_11_to_01\tAB_match_0.75\t"
        "AB_match_0.80\tAB_match_0.85"
    )

    for panel in VCF_BY_PANEL:
        rows = [
            row for row in wide_rows
            if row["panel"] == panel
        ]

        if not rows:
            continue

        print(
            f"{panel}\t"
            f"{len(rows)}\t"
            f"{sum(row['is_representative'] for row in rows)}\t"
            f"{sum(row['exact_left_11_right_01_match'] for row in rows)}\t"
            f"{sum(row['AB_match_threshold_0.75'] for row in rows)}\t"
            f"{sum(row['AB_match_threshold_0.80'] for row in rows)}\t"
            f"{sum(row['AB_match_threshold_0.85'] for row in rows)}"
        )

    print()
    print("REPRESENTATIVE SITES")
    print(
        "group\tposition\ttier\tBA\tleft_GT\tleft_AF\t"
        "right_GT\tright_AF\tformal_match\tAB80_match"
    )

    representatives = [
        row for row in wide_rows
        if row["is_representative"]
    ]

    representatives.sort(
        key=lambda row: (
            list(VCF_BY_PANEL).index(row["panel"]),
            row["pos"],
        )
    )

    for row in representatives:
        print(
            f'{row["group_id"]}\t'
            f'{row["chrom"]}:{row["pos"]}\t'
            f'{row["tier"]}\t'
            f'{row["minimum_BA"]}\t'
            f'{row["left_pseudo_GT"]}\t'
            f'{row["left_alt_fraction"]}\t'
            f'{row["right_pseudo_GT"]}\t'
            f'{row["right_alt_fraction"]}\t'
            f'{row["exact_left_11_right_01_match"]}\t'
            f'{row["AB_match_threshold_0.80"]}'
        )

    print()
    print(f"Long table: {OUTPUT_LONG}")
    print(f"Wide table: {OUTPUT_WIDE}")


if __name__ == "__main__":
    main()
