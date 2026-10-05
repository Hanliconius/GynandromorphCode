#!/usr/bin/env python3

from pathlib import Path
import csv
from collections import OrderedDict


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

VCF = (
    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
    / "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

METADATA = (
    PROJECT
    / "data"
    / "derived"
    / "region_panel_metadata.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis"
    / "original_502"
    / "optix_deltaAF_haplotype_analysis"
)

ALL_SCORE_FILE = OUTDIR / "optix_B_D_R_deltaAF_all_SNPs.tsv"
SELECTED_FILE = OUTDIR / "optix_B_D_R_deltaAF_informative_SNPs.tsv"
SAMPLE_FILE = OUTDIR / "optix_B_D_R_deltaAF_comparative_samples.tsv"
SUMMARY_FILE = OUTDIR / "optix_B_D_R_deltaAF_summary.tsv"

CHROM = "Hmr18"
START = 1060000
END = 1240000

MIN_DP = 2
MIN_GQ = 20
MIN_GROUP_CALL_RATE = 0.80
MIN_ABS_DELTA_AF = 0.80

FOCAL_SAMPLES = (
    "TN_Left",
    "TN_Right",
)


def clean_text(value):
    if value is None:
        return ""

    return " ".join(
        str(value).replace("\r", " ").replace("\n", " ").split()
    )


def parse_binary(value):
    value = clean_text(value)

    if value == "":
        return None

    try:
        numeric = float(value)
    except ValueError:
        return None

    if abs(numeric) < 1e-9:
        return 0

    if abs(numeric - 1) < 1e-9:
        return 1

    return None


def canonical_gt(value):
    if value is None:
        return None

    value = value.replace("|", "/")

    if value in {
        "",
        ".",
        "./.",
        ".|.",
    }:
        return None

    alleles = value.split("/")

    if len(alleles) != 2:
        return None

    if any(
        allele not in {"0", "1"}
        for allele in alleles
    ):
        return None

    return "/".join(
        sorted(alleles)
    )


def numeric(value):
    if value in {
        None,
        "",
        ".",
    }:
        return None

    try:
        return float(value)
    except ValueError:
        return None


def parse_ad(value):
    if value in {
        None,
        "",
        ".",
    }:
        return None

    fields = value.split(",")

    if len(fields) < 2:
        return None

    try:
        ref_depth = int(fields[0])
        alt_depth = int(fields[1])
    except ValueError:
        return None

    if ref_depth < 0 or alt_depth < 0:
        return None

    return ref_depth, alt_depth


def genotype_alt_dosage(gt):
    if gt == "0/0":
        return 0

    if gt == "0/1":
        return 1

    if gt == "1/1":
        return 2

    return None


def read_vcf_samples():
    with VCF.open() as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]

    raise RuntimeError(
        "VCF #CHROM header was not found."
    )


def load_comparative_groups(vcf_samples):
    vcf_sample_set = set(vcf_samples)

    with METADATA.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )
        metadata_rows = list(reader)

    logical_records = OrderedDict()

    for row in metadata_rows:
        accession = clean_text(
            row.get("Accession")
        )

        if accession == "":
            continue

        if accession not in logical_records:
            logical_records[accession] = row

    b_group = []
    d_group = []
    r_group = []
    sample_records = {}

    for accession, row in logical_records.items():
        if accession not in vcf_sample_set:
            continue

        taxon = clean_text(
            row.get("Taxon Name")
        ).lower()

        b_value = parse_binary(
            row.get("B")
        )

        d_value = parse_binary(
            row.get("D")
        )

        r_value = parse_binary(
            row.get("R")
        )

        if (
            (
                "heurippa" in taxon
                or "besckei" in taxon
            )
            and b_value == 1
            and d_value == 0
            and r_value == 0
        ):
            b_group.append(accession)
            group = "B_only"

        elif (
            "meriana" in taxon
            and b_value == 0
            and d_value == 1
            and r_value == 0
        ):
            d_group.append(accession)
            group = "D_only_meriana"

        elif (
            "contigua" in taxon
            and b_value == 0
            and d_value == 0
            and r_value == 1
        ):
            r_group.append(accession)
            group = "R_only_contigua"

        else:
            continue

        sample_records[accession] = {
            "sample": accession,
            "phenotype_group": group,
            "taxon": clean_text(
                row.get("Taxon Name")
            ),
            "data_type": clean_text(
                row.get("Data Type")
            ),
            "B": b_value,
            "D": d_value,
            "R": r_value,
        }

    b_group = sorted(set(b_group))
    d_group = sorted(set(d_group))
    r_group = sorted(set(r_group))

    if len(b_group) < 10:
        raise RuntimeError(
            f"Only {len(b_group)} B-only heurippa/besckei "
            "samples were recovered."
        )

    if len(d_group) < 5:
        raise RuntimeError(
            f"Only {len(d_group)} D-only meriana "
            "samples were recovered."
        )

    if len(r_group) < 8:
        raise RuntimeError(
            f"Only {len(r_group)} R-only contigua "
            "samples were recovered."
        )

    return (
        b_group,
        d_group,
        r_group,
        sample_records,
    )


def parse_sample_call(
    sample_field,
    format_index,
):
    values = sample_field.split(":")

    def get_value(name):
        index = format_index.get(name)

        if index is None:
            return None

        if index >= len(values):
            return None

        return values[index]

    gt = canonical_gt(
        get_value("GT")
    )

    dp = numeric(
        get_value("DP")
    )

    gq = numeric(
        get_value("GQ")
    )

    ad_text = get_value("AD")
    ad = parse_ad(ad_text)

    return {
        "GT": gt,
        "DP": dp,
        "GQ": gq,
        "AD": ad,
        "AD_text": (
            ad_text
            if ad_text not in {
                None,
                ".",
            }
            else ""
        ),
    }


def qc_pass(call):
    return (
        call["GT"] is not None
        and call["DP"] is not None
        and call["GQ"] is not None
        and call["DP"] >= MIN_DP
        and call["GQ"] >= MIN_GQ
    )


def group_allele_frequency(
    samples,
    calls,
):
    called = []

    for sample in samples:
        call = calls[sample]

        if not qc_pass(call):
            continue

        dosage = genotype_alt_dosage(
            call["GT"]
        )

        if dosage is None:
            continue

        called.append(dosage)

    called_n = len(called)
    total_n = len(samples)

    if called_n == 0:
        alt_af = None
    else:
        alt_af = (
            sum(called)
            / (2 * called_n)
        )

    return {
        "samples_n": total_n,
        "called_n": called_n,
        "call_rate": called_n / total_n,
        "alt_AF": alt_af,
    }


def dr_focal_values(
    call,
    dr_allele,
):
    gt_dosage = None
    dr_fraction = None

    alt_dosage = genotype_alt_dosage(
        call["GT"]
    )

    if alt_dosage is not None:
        if dr_allele == "ALT":
            gt_dosage = (
                alt_dosage / 2
            )
        else:
            gt_dosage = (
                2 - alt_dosage
            ) / 2

    if call["AD"] is not None:
        ref_depth, alt_depth = call["AD"]
        total_depth = ref_depth + alt_depth

        if total_depth > 0:
            if dr_allele == "ALT":
                dr_fraction = (
                    alt_depth / total_depth
                )
            else:
                dr_fraction = (
                    ref_depth / total_depth
                )

    return {
        "GT": call["GT"],
        "AD": call["AD_text"],
        "DP": call["DP"],
        "GQ": call["GQ"],
        "QC": qc_pass(call),
        "DR_GT_dosage": gt_dosage,
        "DR_AD_fraction": dr_fraction,
    }


def write_rows(
    filename,
    rows,
    fieldnames,
):
    with filename.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(
                {
                    key: (
                        ""
                        if row.get(key) is None
                        else row.get(key)
                    )
                    for key in fieldnames
                }
            )


def main():
    if not VCF.is_file():
        raise RuntimeError(
            f"Missing VCF: {VCF}"
        )

    if not METADATA.is_file():
        raise RuntimeError(
            f"Missing metadata: {METADATA}"
        )

    OUTDIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    vcf_samples = read_vcf_samples()

    sample_to_index = {
        sample: index
        for index, sample in enumerate(
            vcf_samples
        )
    }

    (
        b_group,
        d_group,
        r_group,
        sample_records,
    ) = load_comparative_groups(
        vcf_samples
    )

    comparisons = {
        "Band": {
            "positive": b_group,
            "negative": sorted(
                set(
                    d_group
                    + r_group
                )
            ),
            "dr_enriched_side": "negative",
            "positive_label": "B-only heurippa/besckei",
            "negative_label": "D-only meriana + R-only contigua",
        },
        "Dennis": {
            "positive": d_group,
            "negative": sorted(
                set(
                    b_group
                    + r_group
                )
            ),
            "dr_enriched_side": "positive",
            "positive_label": "D-only meriana",
            "negative_label": "B-only + R-only",
        },
        "Ray": {
            "positive": r_group,
            "negative": sorted(
                set(
                    b_group
                    + d_group
                )
            ),
            "dr_enriched_side": "positive",
            "positive_label": "R-only contigua",
            "negative_label": "B-only + D-only",
        },
    }

    required_samples = sorted(
        set(
            b_group
            + d_group
            + r_group
            + list(FOCAL_SAMPLES)
        )
    )

    missing_samples = sorted(
        set(required_samples)
        - set(vcf_samples)
    )

    if missing_samples:
        raise RuntimeError(
            "Required samples missing from VCF: "
            + ", ".join(missing_samples)
        )

    sample_manifest_rows = []

    for sample in sorted(sample_records):
        sample_manifest_rows.append(
            sample_records[sample]
        )

    write_rows(
        SAMPLE_FILE,
        sample_manifest_rows,
        [
            "sample",
            "phenotype_group",
            "taxon",
            "data_type",
            "B",
            "D",
            "R",
        ],
    )

    required_indices = {
        sample: sample_to_index[sample]
        for sample in required_samples
    }

    all_rows = []
    selected_rows = []
    physical_snps = 0

    with VCF.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip(
                "\n"
            ).split("\t")

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
                or alt in {
                    ".",
                    "*",
                }
            ):
                continue

            format_names = fields[8].split(
                ":"
            )

            format_index = {
                name: index
                for index, name in enumerate(
                    format_names
                )
            }

            if "GT" not in format_index:
                continue

            calls = {}

            for sample, index in required_indices.items():
                calls[sample] = parse_sample_call(
                    fields[9 + index],
                    format_index,
                )

            physical_snps += 1

            for module, specification in comparisons.items():
                positive = group_allele_frequency(
                    specification["positive"],
                    calls,
                )

                negative = group_allele_frequency(
                    specification["negative"],
                    calls,
                )

                if (
                    positive["alt_AF"] is None
                    or negative["alt_AF"] is None
                ):
                    delta_alt_af = None
                    abs_delta_af = None
                    dr_allele = None
                    dr_delta_af = None
                else:
                    delta_alt_af = (
                        positive["alt_AF"]
                        - negative["alt_AF"]
                    )

                    abs_delta_af = abs(
                        delta_alt_af
                    )

                    if (
                        specification[
                            "dr_enriched_side"
                        ]
                        == "positive"
                    ):
                        target_af = positive[
                            "alt_AF"
                        ]
                        other_af = negative[
                            "alt_AF"
                        ]
                    else:
                        target_af = negative[
                            "alt_AF"
                        ]
                        other_af = positive[
                            "alt_AF"
                        ]

                    if target_af >= other_af:
                        dr_allele = "ALT"
                        dr_delta_af = (
                            target_af - other_af
                        )
                    else:
                        dr_allele = "REF"
                        dr_delta_af = (
                            (
                                1 - target_af
                            )
                            - (
                                1 - other_af
                            )
                        )

                    dr_delta_af = abs(
                        dr_delta_af
                    )

                informative = (
                    abs_delta_af is not None
                    and abs_delta_af
                        >= MIN_ABS_DELTA_AF
                    and positive["call_rate"]
                        >= MIN_GROUP_CALL_RATE
                    and negative["call_rate"]
                        >= MIN_GROUP_CALL_RATE
                )

                row = {
                    "module": module,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "positive_label": specification[
                        "positive_label"
                    ],
                    "negative_label": specification[
                        "negative_label"
                    ],
                    "positive_n": positive[
                        "samples_n"
                    ],
                    "negative_n": negative[
                        "samples_n"
                    ],
                    "positive_called": positive[
                        "called_n"
                    ],
                    "negative_called": negative[
                        "called_n"
                    ],
                    "positive_call_rate": positive[
                        "call_rate"
                    ],
                    "negative_call_rate": negative[
                        "call_rate"
                    ],
                    "positive_ALT_AF": positive[
                        "alt_AF"
                    ],
                    "negative_ALT_AF": negative[
                        "alt_AF"
                    ],
                    "delta_ALT_AF": delta_alt_af,
                    "abs_delta_AF": abs_delta_af,
                    "DR_allele": dr_allele,
                    "DR_delta_AF": dr_delta_af,
                    "informative": informative,
                }

                all_rows.append(row)

                if not informative:
                    continue

                left = dr_focal_values(
                    calls["TN_Left"],
                    dr_allele,
                )

                right = dr_focal_values(
                    calls["TN_Right"],
                    dr_allele,
                )

                selected_row = dict(row)

                selected_row.update(
                    {
                        "TN_Left_GT": left["GT"],
                        "TN_Left_AD": left["AD"],
                        "TN_Left_DP": left["DP"],
                        "TN_Left_GQ": left["GQ"],
                        "TN_Left_QC": left["QC"],
                        "TN_Left_DR_GT_dosage": left[
                            "DR_GT_dosage"
                        ],
                        "TN_Left_DR_AD_fraction": left[
                            "DR_AD_fraction"
                        ],
                        "TN_Right_GT": right["GT"],
                        "TN_Right_AD": right["AD"],
                        "TN_Right_DP": right["DP"],
                        "TN_Right_GQ": right["GQ"],
                        "TN_Right_QC": right["QC"],
                        "TN_Right_DR_GT_dosage": right[
                            "DR_GT_dosage"
                        ],
                        "TN_Right_DR_AD_fraction": right[
                            "DR_AD_fraction"
                        ],
                        "exact_left1_right05": (
                            left["QC"]
                            and right["QC"]
                            and left["DR_GT_dosage"]
                                == 1
                            and right["DR_GT_dosage"]
                                == 0.5
                        ),
                    }
                )

                selected_rows.append(
                    selected_row
                )

            if physical_snps % 10000 == 0:
                print(
                    f"Processed {physical_snps:,} "
                    "biallelic SNPs",
                    flush=True,
                )

    all_fields = [
        "module",
        "chrom",
        "pos",
        "ref",
        "alt",
        "positive_label",
        "negative_label",
        "positive_n",
        "negative_n",
        "positive_called",
        "negative_called",
        "positive_call_rate",
        "negative_call_rate",
        "positive_ALT_AF",
        "negative_ALT_AF",
        "delta_ALT_AF",
        "abs_delta_AF",
        "DR_allele",
        "DR_delta_AF",
        "informative",
    ]

    selected_fields = all_fields + [
        "TN_Left_GT",
        "TN_Left_AD",
        "TN_Left_DP",
        "TN_Left_GQ",
        "TN_Left_QC",
        "TN_Left_DR_GT_dosage",
        "TN_Left_DR_AD_fraction",
        "TN_Right_GT",
        "TN_Right_AD",
        "TN_Right_DP",
        "TN_Right_GQ",
        "TN_Right_QC",
        "TN_Right_DR_GT_dosage",
        "TN_Right_DR_AD_fraction",
        "exact_left1_right05",
    ]

    write_rows(
        ALL_SCORE_FILE,
        all_rows,
        all_fields,
    )

    write_rows(
        SELECTED_FILE,
        selected_rows,
        selected_fields,
    )

    summary_rows = []

    for module in (
        "Band",
        "Dennis",
        "Ray",
    ):
        module_rows = [
            row
            for row in selected_rows
            if row["module"] == module
        ]

        both_qc = [
            row
            for row in module_rows
            if row["TN_Left_QC"]
            and row["TN_Right_QC"]
        ]

        exact = [
            row
            for row in both_qc
            if row["exact_left1_right05"]
        ]

        left_af = [
            row["TN_Left_DR_AD_fraction"]
            for row in both_qc
            if row[
                "TN_Left_DR_AD_fraction"
            ] is not None
        ]

        right_af = [
            row["TN_Right_DR_AD_fraction"]
            for row in both_qc
            if row[
                "TN_Right_DR_AD_fraction"
            ] is not None
        ]

        def median(values):
            if not values:
                return None

            ordered = sorted(values)
            midpoint = len(ordered) // 2

            if len(ordered) % 2:
                return ordered[midpoint]

            return (
                ordered[midpoint - 1]
                + ordered[midpoint]
            ) / 2

        summary_rows.append(
            {
                "module": module,
                "informative_SNPs": len(
                    module_rows
                ),
                "both_focal_QC": len(
                    both_qc
                ),
                "exact_left1_right05": len(
                    exact
                ),
                "exact_fraction": (
                    len(exact) / len(both_qc)
                    if both_qc
                    else None
                ),
                "median_left_DR_AD_fraction": median(
                    left_af
                ),
                "median_right_DR_AD_fraction": median(
                    right_af
                ),
            }
        )

    write_rows(
        SUMMARY_FILE,
        summary_rows,
        [
            "module",
            "informative_SNPs",
            "both_focal_QC",
            "exact_left1_right05",
            "exact_fraction",
            "median_left_DR_AD_fraction",
            "median_right_DR_AD_fraction",
        ],
    )

    print()
    print("OPTIX COMPARATIVE DELTA-AF HAPLOTYPE ANALYSIS")
    print()
    print("COMPARATIVE GROUPS")
    print(
        f"B-only heurippa/besckei: {len(b_group)}"
    )
    print(
        f"D-only meriana: {len(d_group)}"
    )
    print(
        f"R-only contigua: {len(r_group)}"
    )

    print()
    print("SELECTION RULE")
    print(
        f"|delta AF| >= {MIN_ABS_DELTA_AF:.2f}"
    )
    print(
        f"Both comparative groups >= "
        f"{MIN_GROUP_CALL_RATE:.2f} callable"
    )
    print(
        f"Calls require DP >= {MIN_DP} "
        f"and GQ >= {MIN_GQ}"
    )
    print(
        "TN_Left and TN_Right were not used "
        "to select or polarize SNPs."
    )

    print()
    print("FOCAL SUMMARY")
    print(
        "module\tinformative\tboth_QC\texact_L1_R0.5\t"
        "exact_fraction\tmedian_left_AF\tmedian_right_AF"
    )

    for row in summary_rows:
        print(
            f'{row["module"]}\t'
            f'{row["informative_SNPs"]}\t'
            f'{row["both_focal_QC"]}\t'
            f'{row["exact_left1_right05"]}\t'
            f'{row["exact_fraction"]}\t'
            f'{row["median_left_DR_AD_fraction"]}\t'
            f'{row["median_right_DR_AD_fraction"]}'
        )

    print()
    print(f"All SNPs: {ALL_SCORE_FILE}")
    print(f"Informative SNPs: {SELECTED_FILE}")
    print(f"Comparative samples: {SAMPLE_FILE}")
    print(f"Summary: {SUMMARY_FILE}")


if __name__ == "__main__":
    main()
