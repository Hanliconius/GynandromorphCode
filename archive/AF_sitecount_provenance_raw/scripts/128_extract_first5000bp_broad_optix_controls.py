#!/usr/bin/env python3

from pathlib import Path
import re

import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

VCF = (
    PROJECT
    / "external/krista_genotype_plots_original/Genotype_Plots"
    / "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

CONTROL_FILE = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "optix_highest_depth_WGS_controls.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis/original_502/checkpoints"
)

OUT_LONG = (
    OUTDIR
    / "optix_first5000bp_focal_and_WGS_controls.long.tsv"
)

OUT_SELECTION = (
    OUTDIR
    / "optix_first5000bp_selection.tsv"
)

OUT_AUDIT = (
    OUTDIR
    / "optix_first5000bp_sample_genotype_audit.tsv"
)

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

    if (
        len(alleles) != 2
        or any(allele not in {"0", "1"} for allele in alleles)
    ):
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


for required in [VCF, CONTROL_FILE]:
    if not required.is_file():
        raise RuntimeError(f"Missing required input: {required}")

controls = pd.read_csv(
    CONTROL_FILE,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

controls["coverage_rank"] = pd.to_numeric(
    controls["coverage_rank"],
    errors="raise",
).astype(int)

controls = controls.sort_values("coverage_rank").head(3)

if len(controls) != 3:
    raise RuntimeError(
        "Expected three selected WGS controls."
    )

control_samples = controls["sample"].tolist()
wanted_samples = FOCAL_SAMPLES + control_samples

rows = []
interval_chrom = None
interval_start = None
interval_end = None

with VCF.open(
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

            missing = set(wanted_samples).difference(sample_index)

            if missing:
                raise RuntimeError(
                    "Samples missing from broad VCF: "
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

        if interval_start is None:
            interval_chrom = chrom
            interval_start = pos
            interval_end = interval_start + 4999

        if chrom != interval_chrom or pos > interval_end:
            break

        # Retain every biallelic SNP in the 500-bp interval.
        if (
            len(ref) != 1
            or len(alt) != 1
            or "," in alt
        ):
            continue

        format_names = fields[8].split(":")
        format_index = {
            name: index
            for index, name in enumerate(format_names)
        }

        for sample_order, sample in enumerate(
            wanted_samples,
            start=1,
        ):
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
                ad_total = ref_depth + alt_depth
                alt_fraction = alt_depth / ad_total
            else:
                ad_total = None
                alt_fraction = None

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
                taxon_name = "gynandromorph"
                phenotype = "Left"
                coverage_rank = 0

            elif sample == "TN_Right":
                sample_type = "focal"
                display_name = "TN_Right | pooled right"
                taxon_name = "gynandromorph"
                phenotype = "Right"
                coverage_rank = 0

            else:
                control = controls.loc[
                    controls["sample"] == sample
                ].iloc[0]

                sample_type = "WGS_control"
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

            rows.append(
                {
                    "region": "Broad_Optix_first5000bp",
                    "source_vcf": str(VCF),
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
                    "sample_order": sample_order,
                }
            )

if interval_start is None:
    raise RuntimeError(
        "No variant records were found in the broad VCF."
    )

data = pd.DataFrame(rows)

if data.empty:
    raise RuntimeError(
        "No biallelic SNPs were found in the first 500-bp span."
    )

data = data.sort_values(
    [
        "sample_order",
        "pos",
        "ref",
        "alt",
    ]
)

data.to_csv(
    OUT_LONG,
    sep="\t",
    index=False,
)

audit_rows = []

for sample in wanted_samples:
    sample_data = data.loc[
        data["sample"] == sample
    ]

    qc_data = sample_data.loc[
        sample_data["QC_pass"]
    ]

    counts = (
        qc_data["GT_class"]
        .value_counts()
        .to_dict()
    )

    audit_rows.append(
        {
            "sample": sample,
            "sample_type": sample_data[
                "sample_type"
            ].iloc[0],
            "biallelic_SNPs": sample_data[
                ["chrom", "pos", "ref", "alt"]
            ].drop_duplicates().shape[0],
            "QC_sites": len(qc_data),
            "GT_0_0": counts.get("0/0", 0),
            "GT_0_1": counts.get("0/1", 0),
            "GT_1_1": counts.get("1/1", 0),
            "median_DP": sample_data["DP"].median(),
        }
    )

audit = pd.DataFrame(audit_rows)

audit.to_csv(
    OUT_AUDIT,
    sep="\t",
    index=False,
)

selection = pd.DataFrame(
    [
        {
            "region": (
                f"{interval_chrom}:"
                f"{interval_start}-{interval_end}"
            ),
            "chrom": interval_chrom,
            "start": interval_start,
            "end": interval_end,
            "span_bp": 5000,
            "biallelic_SNPs": data[
                ["chrom", "pos", "ref", "alt"]
            ].drop_duplicates().shape[0],
            "selection_basis": (
                "first_5000bp_from_first_broad_VCF_record"
            ),
            "controls_selected_by": (
                "highest_median_DP_WGS_only"
            ),
        }
    ]
)

selection.to_csv(
    OUT_SELECTION,
    sep="\t",
    index=False,
)

print("FIRST 5000 BP OF BROAD OPTIX VCF")
print()
print(selection.to_string(index=False))
print()
print("SAMPLE GENOTYPE AUDIT")
print(audit.to_string(index=False))
print()
print("The three controls are WGS only.")
print("TN_Left and TN_Right did not influence control selection.")
print()
print(f"Plot table: {OUT_LONG}")
print(f"Selection: {OUT_SELECTION}")
print(f"Audit: {OUT_AUDIT}")
