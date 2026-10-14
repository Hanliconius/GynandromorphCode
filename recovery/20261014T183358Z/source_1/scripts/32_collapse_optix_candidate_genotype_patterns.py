#!/usr/bin/env python3

from pathlib import Path
import csv
import re


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

CANDIDATE_FILE = (
    PROJECT
    / "analysis/checkpoints/optix_51sample_stable_candidate_sites.tsv"
)

OUTPUT_GROUPS = (
    PROJECT
    / "analysis/checkpoints/optix_51sample_candidate_pattern_groups.tsv"
)

OUTPUT_REPRESENTATIVES = (
    PROJECT
    / "analysis/checkpoints/optix_51sample_representative_sites.tsv"
)

OUTPUT_FOCAL = (
    PROJECT
    / "analysis/checkpoints/optix_51sample_representative_focal_legs.tsv"
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

FOCAL_SAMPLES = [
    "TN2L_S1",
    "TN3L_S3",
    "TN2R_S2",
    "TN3R_S4",
]

LEFT_SAMPLES = ["TN2L_S1", "TN3L_S3"]
RIGHT_SAMPLES = ["TN2R_S2", "TN3R_S4"]

EXPECTED_FOCAL_STATES = {
    "B": {"left": "0", "right": "1"},
    "D": {"left": "1", "right": "1"},
    "R": {"left": "1", "right": "1"},
}

MINIMUM_JOINTLY_CALLED_FOR_GROUPING = 20


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

    if not alleles or any(
        allele in {"", "."} for allele in alleles
    ):
        return None

    try:
        numeric = [int(allele) for allele in alleles]
    except ValueError:
        return None

    return "/".join(str(x) for x in sorted(numeric))


def phenotype_prediction(gt, orientation):
    gt = canonical_gt(gt)

    if gt is None:
        return None

    alleles = [int(x) for x in gt.split("/")]

    if orientation == "ALT_carrier_predicts_presence":
        return "1" if any(x > 0 for x in alleles) else "0"

    if orientation == "REF_carrier_predicts_presence":
        return "1" if any(x == 0 for x in alleles) else "0"

    return None


def parse_number(value):
    if value in {None, "", "."}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_sample_call(format_ids, sample_field):
    values = sample_field.split(":")
    parsed = dict(zip(format_ids, values))

    return {
        "GT": canonical_gt(parsed.get("GT")),
        "DP": parse_number(parsed.get("DP")),
        "GQ": parse_number(parsed.get("GQ")),
    }


def qc_prediction(call, orientation):
    if call["GT"] is None:
        return None

    if call["DP"] is None or call["DP"] < 2:
        return None

    if call["GQ"] is None or call["GQ"] < 20:
        return None

    return phenotype_prediction(call["GT"], orientation)


def load_candidate_rows():
    rows = []

    with CANDIDATE_FILE.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            if not (
                as_bool(row["tier_1"])
                or as_bool(row["tier_2"])
            ):
                continue

            row["pos"] = int(row["pos"])
            row["tier_1"] = as_bool(row["tier_1"])
            row["tier_2"] = as_bool(row["tier_2"])

            for field in [
                "raw_called_n",
                "QC_called_n",
                "raw_errors",
                "QC_errors",
            ]:
                row[field] = as_int(row[field])

            for field in [
                "raw_BA",
                "QC_BA",
                "minimum_BA",
                "raw_MCC",
                "QC_MCC",
                "minimum_MCC",
                "raw_q",
                "QC_q",
                "raw_missingness_q",
                "QC_missingness_q",
            ]:
                row[field] = as_float(row[field])

            rows.append(row)

    return rows


def load_vcf_sites(vcf_path, wanted_keys):
    samples = []
    records = {}

    with vcf_path.open() as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            key = (chrom, pos, ref, alt)

            if key not in wanted_keys:
                continue

            format_ids = fields[8].split(":")
            calls = {}

            for sample, sample_field in zip(
                samples, fields[9:]
            ):
                calls[sample] = parse_sample_call(
                    format_ids,
                    sample_field,
                )

            records[key] = calls

    return samples, records


def site_rank(row):
    tier_rank = 0 if row["tier_1"] else 1

    minimum_ba = (
        row["minimum_BA"]
        if row["minimum_BA"] is not None
        else -1.0
    )

    minimum_mcc = (
        row["minimum_MCC"]
        if row["minimum_MCC"] is not None
        else -2.0
    )

    qc_called = (
        row["QC_called_n"]
        if row["QC_called_n"] is not None
        else -1
    )

    qc_errors = (
        row["QC_errors"]
        if row["QC_errors"] is not None
        else 999999
    )

    return (
        tier_rank,
        -minimum_ba,
        -minimum_mcc,
        -qc_called,
        qc_errors,
        row["pos"],
    )


def pattern_comparison(pattern_a, pattern_b):
    jointly_called = 0
    discordant = 0

    for value_a, value_b in zip(pattern_a, pattern_b):
        if value_a is None or value_b is None:
            continue

        jointly_called += 1

        if value_a != value_b:
            discordant += 1

    return jointly_called, discordant


def consensus(predictions):
    called = [
        value for value in predictions
        if value is not None
    ]

    if len(called) == 0:
        return "missing"

    if len(called) < len(predictions):
        if len(set(called)) == 1:
            return f"{called[0]}_partial"
        return "discordant_partial"

    if len(set(called)) == 1:
        return called[0]

    return "discordant"


def fmt(value, digits=4):
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def main():
    required = [CANDIDATE_FILE, *VCF_BY_PANEL.values()]
    missing = [path for path in required if not path.is_file()]

    if missing:
        print("ERROR: required files are missing:")
        for path in missing:
            print(path)
        return

    candidates = load_candidate_rows()

    if not candidates:
        print("ERROR: no Tier-1 or Tier-2 candidates found.")
        return

    all_group_rows = []
    representative_rows = []
    focal_rows = []

    for panel in VCF_BY_PANEL:
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

        samples, records = load_vcf_sites(
            VCF_BY_PANEL[panel],
            wanted_keys,
        )

        background_samples = [
            sample for sample in samples
            if sample not in FOCAL_SAMPLES
        ]

        enriched = []

        for row in panel_candidates:
            key = (
                row["chrom"],
                row["pos"],
                row["ref"],
                row["alt"],
            )

            calls = records.get(key)

            if calls is None:
                print(
                    "WARNING: candidate absent from reconstructed "
                    f"VCF: {panel} {key}"
                )
                continue

            orientation = row["orientation"]

            raw_pattern = tuple(
                phenotype_prediction(
                    calls[sample]["GT"],
                    orientation,
                )
                for sample in background_samples
            )

            qc_pattern = tuple(
                qc_prediction(
                    calls[sample],
                    orientation,
                )
                for sample in background_samples
            )

            enriched.append(
                {
                    **row,
                    "calls": calls,
                    "raw_pattern": raw_pattern,
                    "qc_pattern": qc_pattern,
                }
            )

        enriched.sort(key=site_rank)

        groups = []

        for row in enriched:
            assigned_group = None
            assigned_joint = None
            assigned_discordant = None

            for group in groups:
                representative = group["representative"]

                jointly_called, discordant = pattern_comparison(
                    row["qc_pattern"],
                    representative["qc_pattern"],
                )

                if (
                    jointly_called
                    >= MINIMUM_JOINTLY_CALLED_FOR_GROUPING
                    and discordant == 0
                ):
                    assigned_group = group
                    assigned_joint = jointly_called
                    assigned_discordant = discordant
                    break

            if assigned_group is None:
                groups.append(
                    {
                        "representative": row,
                        "members": [
                            {
                                "row": row,
                                "jointly_called": sum(
                                    value is not None
                                    for value in row["qc_pattern"]
                                ),
                                "discordant": 0,
                            }
                        ],
                    }
                )
            else:
                assigned_group["members"].append(
                    {
                        "row": row,
                        "jointly_called": assigned_joint,
                        "discordant": assigned_discordant,
                    }
                )

        for group_number, group in enumerate(groups, start=1):
            representative = group["representative"]
            trait = representative["trait"]

            group_id = (
                f"{panel}_{trait}_G{group_number:02d}"
            )

            member_positions = sorted(
                member["row"]["pos"]
                for member in group["members"]
            )

            calls = representative["calls"]
            orientation = representative["orientation"]

            focal_predictions = {}

            for sample in FOCAL_SAMPLES:
                call = calls.get(
                    sample,
                    {"GT": None, "DP": None, "GQ": None},
                )

                raw_prediction = phenotype_prediction(
                    call["GT"],
                    orientation,
                )

                filtered_prediction = qc_prediction(
                    call,
                    orientation,
                )

                focal_predictions[sample] = (
                    filtered_prediction
                )

                focal_rows.append(
                    {
                        "group_id": group_id,
                        "panel": panel,
                        "trait": trait,
                        "chrom": representative["chrom"],
                        "pos": representative["pos"],
                        "ref": representative["ref"],
                        "alt": representative["alt"],
                        "orientation": orientation,
                        "sample": sample,
                        "side": (
                            "left"
                            if sample in LEFT_SAMPLES
                            else "right"
                        ),
                        "leg": (
                            "T2"
                            if sample.startswith("TN2")
                            else "T3"
                        ),
                        "GT": call["GT"] or "",
                        "DP": fmt(call["DP"]),
                        "GQ": fmt(call["GQ"]),
                        "raw_predicted_phenotype": (
                            raw_prediction or ""
                        ),
                        "DP2_GQ20_predicted_phenotype": (
                            filtered_prediction or ""
                        ),
                    }
                )

            left_consensus = consensus(
                [
                    focal_predictions[sample]
                    for sample in LEFT_SAMPLES
                ]
            )

            right_consensus = consensus(
                [
                    focal_predictions[sample]
                    for sample in RIGHT_SAMPLES
                ]
            )

            expected = EXPECTED_FOCAL_STATES[trait]

            focal_matches_expectation = (
                left_consensus == expected["left"]
                and right_consensus == expected["right"]
            )

            representative_rows.append(
                {
                    "group_id": group_id,
                    "panel": panel,
                    "trait": trait,
                    "group_sites": len(group["members"]),
                    "group_first_pos": min(member_positions),
                    "group_last_pos": max(member_positions),
                    "representative_chrom": (
                        representative["chrom"]
                    ),
                    "representative_pos": (
                        representative["pos"]
                    ),
                    "ref": representative["ref"],
                    "alt": representative["alt"],
                    "orientation": orientation,
                    "tier": (
                        "tier_1"
                        if representative["tier_1"]
                        else "tier_2"
                    ),
                    "minimum_BA": fmt(
                        representative["minimum_BA"]
                    ),
                    "minimum_MCC": fmt(
                        representative["minimum_MCC"]
                    ),
                    "raw_errors": (
                        representative["raw_errors"]
                    ),
                    "QC_errors": (
                        representative["QC_errors"]
                    ),
                    "raw_called_n": (
                        representative["raw_called_n"]
                    ),
                    "QC_called_n": (
                        representative["QC_called_n"]
                    ),
                    "left_leg_consensus_DP2_GQ20": (
                        left_consensus
                    ),
                    "right_leg_consensus_DP2_GQ20": (
                        right_consensus
                    ),
                    "expected_left_state": expected["left"],
                    "expected_right_state": expected["right"],
                    "focal_legs_match_expected_pattern": (
                        focal_matches_expectation
                    ),
                    "interpretation_scope": (
                        "comparative_association; "
                        "four_legs_QC_only; no_pooled_samples"
                    ),
                }
            )

            for member in group["members"]:
                row = member["row"]

                all_group_rows.append(
                    {
                        "group_id": group_id,
                        "panel": panel,
                        "trait": trait,
                        "chrom": row["chrom"],
                        "pos": row["pos"],
                        "ref": row["ref"],
                        "alt": row["alt"],
                        "orientation": row["orientation"],
                        "tier": (
                            "tier_1"
                            if row["tier_1"]
                            else "tier_2"
                        ),
                        "is_representative": (
                            row is representative
                        ),
                        "minimum_BA": fmt(
                            row["minimum_BA"]
                        ),
                        "minimum_MCC": fmt(
                            row["minimum_MCC"]
                        ),
                        "raw_errors": row["raw_errors"],
                        "QC_errors": row["QC_errors"],
                        "QC_called_n": row["QC_called_n"],
                        "jointly_called_vs_representative": (
                            member["jointly_called"]
                        ),
                        "discordant_calls_vs_representative": (
                            member["discordant"]
                        ),
                    }
                )

    OUTPUT_GROUPS.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    group_fields = [
        "group_id",
        "panel",
        "trait",
        "chrom",
        "pos",
        "ref",
        "alt",
        "orientation",
        "tier",
        "is_representative",
        "minimum_BA",
        "minimum_MCC",
        "raw_errors",
        "QC_errors",
        "QC_called_n",
        "jointly_called_vs_representative",
        "discordant_calls_vs_representative",
    ]

    representative_fields = [
        "group_id",
        "panel",
        "trait",
        "group_sites",
        "group_first_pos",
        "group_last_pos",
        "representative_chrom",
        "representative_pos",
        "ref",
        "alt",
        "orientation",
        "tier",
        "minimum_BA",
        "minimum_MCC",
        "raw_errors",
        "QC_errors",
        "raw_called_n",
        "QC_called_n",
        "left_leg_consensus_DP2_GQ20",
        "right_leg_consensus_DP2_GQ20",
        "expected_left_state",
        "expected_right_state",
        "focal_legs_match_expected_pattern",
        "interpretation_scope",
    ]

    focal_fields = [
        "group_id",
        "panel",
        "trait",
        "chrom",
        "pos",
        "ref",
        "alt",
        "orientation",
        "sample",
        "side",
        "leg",
        "GT",
        "DP",
        "GQ",
        "raw_predicted_phenotype",
        "DP2_GQ20_predicted_phenotype",
    ]

    with OUTPUT_GROUPS.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=group_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(all_group_rows)

    with OUTPUT_REPRESENTATIVES.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=representative_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(representative_rows)

    with OUTPUT_FOCAL.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=focal_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(focal_rows)

    print("OPTIX CANDIDATE GENOTYPE-PATTERN COLLAPSE")
    print()
    print(
        "Sites are grouped when their DP2/GQ20 dominant "
        "genotype predictions have zero disagreements across "
        f"at least {MINIMUM_JOINTLY_CALLED_FOR_GROUPING} "
        "jointly called non-focal samples."
    )
    print(
        "Missing calls may differ within a group. Four leg "
        "samples are reported but do not affect grouping or "
        "representative selection."
    )
    print()

    print(
        "panel\tcandidate_sites\tpattern_groups\t"
        "representative_tier1\t"
        "representatives_matching_focal_expectation"
    )

    panels = sorted({
        row["panel"] for row in representative_rows
    })

    for panel in panels:
        panel_groups = [
            row for row in representative_rows
            if row["panel"] == panel
        ]

        panel_sites = [
            row for row in all_group_rows
            if row["panel"] == panel
        ]

        tier1_n = sum(
            row["tier"] == "tier_1"
            for row in panel_groups
        )

        focal_match_n = sum(
            row["focal_legs_match_expected_pattern"]
            for row in panel_groups
        )

        print(
            f"{panel}\t{len(panel_sites)}\t"
            f"{len(panel_groups)}\t{tier1_n}\t"
            f"{focal_match_n}"
        )

    print()
    print("REPRESENTATIVE SITES")
    print(
        "group\tposition\ttier\tsites\tmin_BA\tmin_MCC\t"
        "raw_errors\tQC_errors\tleft_legs\tright_legs\t"
        "matches_expected_focal_pattern"
    )

    for row in representative_rows:
        print(
            f'{row["group_id"]}\t'
            f'{row["representative_chrom"]}:'
            f'{row["representative_pos"]}\t'
            f'{row["tier"]}\t'
            f'{row["group_sites"]}\t'
            f'{row["minimum_BA"]}\t'
            f'{row["minimum_MCC"]}\t'
            f'{row["raw_errors"]}\t'
            f'{row["QC_errors"]}\t'
            f'{row["left_leg_consensus_DP2_GQ20"]}\t'
            f'{row["right_leg_consensus_DP2_GQ20"]}\t'
            f'{row["focal_legs_match_expected_pattern"]}'
        )

    print()
    print(f"Pattern groups: {OUTPUT_GROUPS}")
    print(f"Representatives: {OUTPUT_REPRESENTATIVES}")
    print(f"Focal leg calls: {OUTPUT_FOCAL}")


if __name__ == "__main__":
    main()
