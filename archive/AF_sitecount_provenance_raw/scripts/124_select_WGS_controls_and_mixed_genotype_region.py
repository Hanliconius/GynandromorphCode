#!/usr/bin/env python3

from pathlib import Path
import re

import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

BA90_TABLE = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "optix_R_all_BA90_comparative_selected_genotypes.long.tsv"
)

VCFS = {
    "Band": (
        PROJECT
        / "external/krista_genotype_plots_original/Genotype_Plots"
        / "optix_band-1112651-1114314.vcf"
    ),
    "Dennis1": (
        PROJECT
        / "external/krista_genotype_plots_original/Genotype_Plots"
        / "dennis1-1215210-1216124.vcf"
    ),
    "Dennis2": (
        PROJECT
        / "external/krista_genotype_plots_original/Genotype_Plots"
        / "dennis2-1218557-1221494.vcf"
    ),
    "Ray1": (
        PROJECT
        / "external/krista_genotype_plots_original/Genotype_Plots"
        / "ray1-1156746-1160220.vcf"
    ),
    "Ray2": (
        PROJECT
        / "external/krista_genotype_plots_original/Genotype_Plots"
        / "ray2-1174547-1187833.vcf"
    ),
}

OUTDIR = (
    PROJECT
    / "analysis/original_502/checkpoints"
)

CONTROL_FILE = (
    OUTDIR
    / "optix_highest_depth_WGS_controls.tsv"
)

REGION_SAMPLE_AUDIT = (
    OUTDIR
    / "optix_WGS_control_region_genotype_audit.tsv"
)

REGION_SUMMARY = (
    OUTDIR
    / "optix_WGS_control_region_selection_summary.tsv"
)

SELECTED_LONG = (
    OUTDIR
    / "optix_selected_region_focal_and_WGS_controls.long.tsv"
)

SELECTION_FILE = (
    OUTDIR
    / "optix_selected_mixed_genotype_region.tsv"
)

MINIMUM_PER_GENOTYPE_CLASS = 5

FOCAL_SAMPLES = [
    "TN_Left",
    "TN_Right",
]


def clean(value):
    return re.sub(r"\s+", " ", str(value)).strip()


def numeric(value):
    try:
        if value in {"", ".", None}:
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_gt(value):
    if value in {"", ".", "./.", ".|."}:
        return None

    alleles = re.split(r"[/|]", value)

    if len(alleles) != 2:
        return None

    if any(allele not in {"0", "1"} for allele in alleles):
        return None

    alleles = sorted(int(allele) for allele in alleles)

    return f"{alleles[0]}/{alleles[1]}"


def parse_ad(value):
    if value in {"", ".", ".,.", None}:
        return None, None

    pieces = str(value).split(",")

    if len(pieces) < 2:
        return None, None

    return numeric(pieces[0]), numeric(pieces[1])


for required in [BA90_TABLE, *VCFS.values()]:
    if not required.is_file():
        raise RuntimeError(f"Missing required input: {required}")

OUTDIR.mkdir(parents=True, exist_ok=True)

existing = pd.read_csv(
    BA90_TABLE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

existing["DP_numeric"] = pd.to_numeric(
    existing["DP"],
    errors="coerce",
)

wgs = existing.loc[
    (existing["sample_type"] == "comparative")
    & (
        existing["data_type"]
        .map(clean)
        .str.upper()
        == "WGS"
    )
].copy()

if wgs.empty:
    raise RuntimeError(
        "No comparative WGS samples were found in the BA90 table."
    )

control_rows = []

for sample, group in wgs.groupby("sample"):
    depths = group["DP_numeric"].dropna()

    if depths.empty:
        continue

    control_rows.append(
        {
            "sample": sample,
            "taxon_name": clean(group["taxon_name"].iloc[0]),
            "phenotype_R": clean(group["phenotype"].iloc[0]),
            "data_type": clean(group["data_type"].iloc[0]),
            "BA90_sites": len(group),
            "sites_with_DP": len(depths),
            "median_DP": float(depths.median()),
            "mean_DP": float(depths.mean()),
            "maximum_DP": float(depths.max()),
        }
    )

control_statistics = pd.DataFrame(control_rows)

control_statistics = control_statistics.sort_values(
    [
        "median_DP",
        "sites_with_DP",
        "mean_DP",
        "sample",
    ],
    ascending=[
        False,
        False,
        False,
        True,
    ],
)

selected_controls = control_statistics.head(3).copy()

if len(selected_controls) != 3:
    raise RuntimeError(
        "Fewer than three WGS control samples were available."
    )

selected_controls["coverage_rank"] = range(
    1,
    len(selected_controls) + 1,
)

selected_controls.to_csv(
    CONTROL_FILE,
    sep="\t",
    index=False,
)

control_samples = selected_controls["sample"].tolist()
wanted_samples = FOCAL_SAMPLES + control_samples

all_region_rows = []
region_record_counts = {}

for region_name, vcf_path in VCFS.items():
    region_rows = []
    biallelic_snp_count = 0

    with vcf_path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        sample_index = None

        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                vcf_samples = fields[9:]

                sample_index = {
                    sample: index
                    for index, sample in enumerate(vcf_samples)
                }

                missing = set(wanted_samples).difference(
                    sample_index
                )

                if missing:
                    raise RuntimeError(
                        f"{region_name}: samples missing from VCF: "
                        + ",".join(sorted(missing))
                    )

                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            # A common 0/0, 0/1, 1/1 allele-fraction scale is
            # unambiguous only for biallelic SNPs.
            if (
                len(ref) != 1
                or len(alt) != 1
                or "," in alt
            ):
                continue

            biallelic_snp_count += 1

            format_names = fields[8].split(":")
            format_index = {
                name: index
                for index, name in enumerate(format_names)
            }

            for sample in wanted_samples:
                values = fields[
                    9 + sample_index[sample]
                ].split(":")

                def get_field(name):
                    index = format_index.get(name)

                    if index is None or index >= len(values):
                        return ""

                    return values[index]

                raw_gt = get_field("GT")
                gt_class = parse_gt(raw_gt)
                ad = get_field("AD")
                dp = numeric(get_field("DP"))
                gq = numeric(get_field("GQ"))
                ref_depth, alt_depth = parse_ad(ad)

                if (
                    ref_depth is not None
                    and alt_depth is not None
                    and ref_depth + alt_depth > 0
                ):
                    alt_fraction = (
                        alt_depth
                        / (ref_depth + alt_depth)
                    )
                    ad_total = ref_depth + alt_depth
                else:
                    alt_fraction = None
                    ad_total = None

                qc_pass = (
                    gt_class is not None
                    and dp is not None
                    and gq is not None
                    and dp >= 2
                    and gq >= 20
                    and alt_fraction is not None
                )

                if sample == "TN_Left":
                    sample_type = "focal"
                    display_name = "TN_Left | pooled left"
                    coverage_rank = 0
                    taxon_name = "gynandromorph"
                    phenotype = "Left"

                elif sample == "TN_Right":
                    sample_type = "focal"
                    display_name = "TN_Right | pooled right"
                    coverage_rank = 0
                    taxon_name = "gynandromorph"
                    phenotype = "Right"

                else:
                    sample_type = "WGS_control"

                    control = selected_controls.loc[
                        selected_controls["sample"] == sample
                    ].iloc[0]

                    coverage_rank = int(
                        control["coverage_rank"]
                    )

                    taxon_name = clean(
                        control["taxon_name"]
                    )

                    phenotype = clean(
                        control["phenotype_R"]
                    )

                    display_name = (
                        f"{sample} | WGS coverage rank "
                        f"{coverage_rank} | "
                        f"{taxon_name} | R={phenotype}"
                    )

                region_rows.append(
                    {
                        "region": region_name,
                        "source_vcf": str(vcf_path),
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "sample": sample,
                        "sample_type": sample_type,
                        "display_name": display_name,
                        "coverage_rank": coverage_rank,
                        "taxon_name": taxon_name,
                        "phenotype": phenotype,
                        "raw_GT": raw_gt,
                        "GT_class": gt_class,
                        "AD": ad,
                        "ref_depth": ref_depth,
                        "alt_depth": alt_depth,
                        "AD_total": ad_total,
                        "DP": dp,
                        "GQ": gq,
                        "alt_fraction": alt_fraction,
                        "QC_pass": qc_pass,
                    }
                )

    region_record_counts[region_name] = biallelic_snp_count
    all_region_rows.extend(region_rows)

all_data = pd.DataFrame(all_region_rows)

audit_rows = []

for region_name in VCFS:
    region_data = all_data.loc[
        all_data["region"] == region_name
    ]

    for sample in wanted_samples:
        sample_data = region_data.loc[
            region_data["sample"] == sample
        ]

        qc_data = sample_data.loc[
            sample_data["QC_pass"]
        ]

        genotype_counts = (
            qc_data["GT_class"]
            .value_counts()
            .to_dict()
        )

        audit_rows.append(
            {
                "region": region_name,
                "sample": sample,
                "sample_type": sample_data[
                    "sample_type"
                ].iloc[0],
                "biallelic_SNPs": region_record_counts[
                    region_name
                ],
                "QC_sites": len(qc_data),
                "call_rate": (
                    len(qc_data)
                    / region_record_counts[region_name]
                ),
                "GT_0_0": genotype_counts.get("0/0", 0),
                "GT_0_1": genotype_counts.get("0/1", 0),
                "GT_1_1": genotype_counts.get("1/1", 0),
            }
        )

audit = pd.DataFrame(audit_rows)

audit.to_csv(
    REGION_SAMPLE_AUDIT,
    sep="\t",
    index=False,
)

region_summary_rows = []

for region_name in VCFS:
    control_audit = audit.loc[
        (audit["region"] == region_name)
        & (audit["sample_type"] == "WGS_control")
    ]

    focal_audit = audit.loc[
        (audit["region"] == region_name)
        & (audit["sample_type"] == "focal")
    ]

    control_class_counts = control_audit[
        [
            "GT_0_0",
            "GT_0_1",
            "GT_1_1",
        ]
    ].to_numpy()

    focal_class_counts = focal_audit[
        [
            "GT_0_0",
            "GT_0_1",
            "GT_1_1",
        ]
    ].to_numpy()

    minimum_control_class_count = int(
        control_class_counts.min()
    )

    minimum_focal_class_count = int(
        focal_class_counts.min()
    )

    region_summary_rows.append(
        {
            "region": region_name,
            "biallelic_SNPs": region_record_counts[
                region_name
            ],
            "minimum_control_genotype_class_count": (
                minimum_control_class_count
            ),
            "all_controls_have_at_least_5_each": (
                minimum_control_class_count
                >= MINIMUM_PER_GENOTYPE_CLASS
            ),
            "minimum_focal_genotype_class_count": (
                minimum_focal_class_count
            ),
            "both_focal_samples_have_all_classes": (
                minimum_focal_class_count >= 1
            ),
            "mean_control_call_rate": (
                control_audit["call_rate"].mean()
            ),
            "mean_focal_call_rate": (
                focal_audit["call_rate"].mean()
            ),
        }
    )

region_summary = pd.DataFrame(region_summary_rows)

eligible = region_summary.loc[
    region_summary["all_controls_have_at_least_5_each"]
].copy()

selection_basis = (
    "shortest_region_with_controls_ge5_each_genotype"
)

if not eligible.empty:
    eligible = eligible.sort_values(
        [
            "biallelic_SNPs",
            "minimum_control_genotype_class_count",
            "mean_control_call_rate",
        ],
        ascending=[
            True,
            False,
            False,
        ],
    )

    selected_region = eligible.iloc[0]["region"]

else:
    selection_basis = (
        "fallback_maximum_minimum_control_genotype_count"
    )

    fallback = region_summary.sort_values(
        [
            "minimum_control_genotype_class_count",
            "mean_control_call_rate",
            "biallelic_SNPs",
        ],
        ascending=[
            False,
            False,
            True,
        ],
    )

    selected_region = fallback.iloc[0]["region"]

region_summary["selected"] = (
    region_summary["region"] == selected_region
)

region_summary["selection_basis"] = selection_basis

region_summary.to_csv(
    REGION_SUMMARY,
    sep="\t",
    index=False,
)

selected_data = all_data.loc[
    all_data["region"] == selected_region
].copy()

sample_order = {
    "TN_Left": 1,
    "TN_Right": 2,
}

for rank, sample in enumerate(
    control_samples,
    start=3,
):
    sample_order[sample] = rank

selected_data["sample_order"] = selected_data[
    "sample"
].map(sample_order)

selected_data = selected_data.sort_values(
    [
        "sample_order",
        "pos",
    ]
)

selected_data.to_csv(
    SELECTED_LONG,
    sep="\t",
    index=False,
)

selection = region_summary.loc[
    region_summary["region"] == selected_region
].copy()

selection.to_csv(
    SELECTION_FILE,
    sep="\t",
    index=False,
)

print("WGS CONTROL AND MIXED-GENOTYPE REGION SELECTION")
print()
print("HIGH-DEPTH WGS CONTROLS")
print(
    selected_controls[
        [
            "coverage_rank",
            "sample",
            "taxon_name",
            "phenotype_R",
            "median_DP",
            "mean_DP",
        ]
    ].to_string(index=False)
)

print()
print("REGION SUMMARY")
print(region_summary.to_string(index=False))

print()
print(f"Selected region: {selected_region}")
print(f"Selection basis: {selection_basis}")
print()
print("Selection used WGS controls only.")
print("TN_Left and TN_Right did not influence region selection.")
print()
print(f"Plot table: {SELECTED_LONG}")
print(f"Control table: {CONTROL_FILE}")
print(f"Region audit: {REGION_SAMPLE_AUDIT}")
print(f"Region summary: {REGION_SUMMARY}")
