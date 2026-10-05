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
    / "optix_fullVCF_focal_and_WGS_control_AF.long.tsv.gz"
)

OUT_SUMMARY = (
    OUTDIR
    / "optix_fullVCF_focal_and_WGS_control_AF.summary.tsv"
)


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

OUTDIR.mkdir(parents=True, exist_ok=True)

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
        "Expected exactly three WGS controls."
    )

wanted_samples = [
    "TN_Left",
    "TN_Right",
    *controls["sample"].tolist(),
]

sample_information = {
    "TN_Left": {
        "sample_type": "focal",
        "display_name": "TN_Left | pooled left",
        "sample_order": 1,
        "phenotype": "Left",
        "taxon_name": "gynandromorph",
    },
    "TN_Right": {
        "sample_type": "focal",
        "display_name": "TN_Right | pooled right",
        "sample_order": 2,
        "phenotype": "Right",
        "taxon_name": "gynandromorph",
    },
}

for row in controls.itertuples(index=False):
    sample_information[row.sample] = {
        "sample_type": "WGS_control",
        "display_name": (
            f"{row.sample} | WGS coverage rank "
            f"{row.coverage_rank} | "
            f"{clean(row.taxon_name)} | "
            f"R={clean(row.phenotype_R)}"
        ),
        "sample_order": int(row.coverage_rank) + 2,
        "phenotype": clean(row.phenotype_R),
        "taxon_name": clean(row.taxon_name),
    }

rows = []
biallelic_snps = 0
first_position = None
last_position = None

with VCF.open(
    "r",
    encoding="utf-8",
    errors="replace",
) as handle:
    sample_index = None

    for line_number, line in enumerate(handle, start=1):
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

        if (
            len(ref) != 1
            or len(alt) != 1
            or "," in alt
        ):
            continue

        biallelic_snps += 1

        if first_position is None:
            first_position = pos

        last_position = pos

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

            information = sample_information[sample]

            rows.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "sample": sample,
                    "sample_type": information["sample_type"],
                    "display_name": information["display_name"],
                    "sample_order": information["sample_order"],
                    "phenotype": information["phenotype"],
                    "taxon_name": information["taxon_name"],
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

        if biallelic_snps % 10000 == 0:
            print(
                f"Processed {biallelic_snps:,} biallelic SNPs",
                flush=True,
            )

data = pd.DataFrame(rows)

data.to_csv(
    OUT_LONG,
    sep="\t",
    index=False,
    compression="gzip",
)

summary_rows = []

for sample in wanted_samples:
    sample_data = data.loc[
        data["sample"] == sample
    ]

    qc = sample_data.loc[
        sample_data["QC_pass"]
    ].copy()

    canonical_distance = qc["alt_fraction"].map(
        lambda value: min(
            abs(value - state)
            for state in (0, 0.5, 1)
        )
    )

    gt_counts = (
        qc["GT_class"]
        .value_counts()
        .to_dict()
    )

    het = qc.loc[
        qc["GT_class"] == "0/1",
        "alt_fraction",
    ]

    summary_rows.append(
        {
            "sample": sample,
            "sample_type": sample_information[
                sample
            ]["sample_type"],
            "display_name": sample_information[
                sample
            ]["display_name"],
            "biallelic_SNPs": biallelic_snps,
            "QC_sites": len(qc),
            "QC_fraction": len(qc) / biallelic_snps,
            "GT_0_0": gt_counts.get("0/0", 0),
            "GT_0_1": gt_counts.get("0/1", 0),
            "GT_1_1": gt_counts.get("1/1", 0),
            "median_DP": qc["DP"].median(),
            "median_heterozygote_AF": (
                het.median() if len(het) else None
            ),
            "heterozygote_AF_Q25": (
                het.quantile(0.25) if len(het) else None
            ),
            "heterozygote_AF_Q75": (
                het.quantile(0.75) if len(het) else None
            ),
            "median_distance_from_0_0.5_1": (
                canonical_distance.median()
            ),
            "fraction_within_0.10_of_0_0.5_1": (
                (canonical_distance <= 0.10).mean()
            ),
        }
    )

summary = pd.DataFrame(summary_rows)

summary.to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False,
)

print()
print("FULL OPTIX VCF ALLELE-FRACTION EXTRACTION")
print()
print(f"Range: Hmr18:{first_position}-{last_position}")
print(f"Biallelic SNPs: {biallelic_snps:,}")
print()
print(summary.to_string(index=False))
print()
print(f"Long table: {OUT_LONG}")
print(f"Summary: {OUT_SUMMARY}")
