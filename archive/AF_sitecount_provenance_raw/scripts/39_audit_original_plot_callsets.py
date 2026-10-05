#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict
import csv
import tarfile


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

SOURCE_DIR = (
    PROJECT
    / "external/krista_genotype_plots_original/"
    "Genotype_Plots"
)

ARCHIVE = PROJECT / "Genotype_Plots.tar"

OUTPUT_SUMMARY = (
    PROJECT
    / "analysis/checkpoints/"
    "original_plot_callset_summary.tsv"
)

OUTPUT_SAMPLES = (
    PROJECT
    / "analysis/checkpoints/"
    "original_plot_callset_samples.tsv"
)

BROAD_MEMBER = (
    "Genotype_Plots/"
    "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

VCFS = {
    "band": SOURCE_DIR / "optix_band-1112651-1114314.vcf",
    "dennis1": SOURCE_DIR / "dennis1-1215210-1216124.vcf",
    "dennis2": SOURCE_DIR / "dennis2-1218557-1221494.vcf",
    "ray1": SOURCE_DIR / "ray1-1156746-1160220.vcf",
    "ray2": SOURCE_DIR / "ray2-1174547-1187833.vcf",
    "ivory": (
        SOURCE_DIR
        / "Ivory_Hmr15-1627000-1637000_biallelic.vcf"
    ),
}

FOCAL_NAMES = [
    "TN2L_S1",
    "TN3L_S3",
    "TN2R_S2",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
]


def parse_vcf(path):
    samples = []
    records = 0
    biallelic_snps = 0
    ranges = defaultdict(
        lambda: [None, None]
    )

    with path.open(errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            records += 1

            if ranges[chrom][0] is None:
                ranges[chrom][0] = pos

            ranges[chrom][1] = pos

            if (
                "," not in alt
                and len(ref) == 1
                and len(alt) == 1
            ):
                biallelic_snps += 1

    range_text = ";".join(
        f"{chrom}:{values[0]}-{values[1]}"
        for chrom, values in sorted(ranges.items())
    )

    return {
        "samples": samples,
        "records": records,
        "biallelic_snps": biallelic_snps,
        "ranges": range_text,
    }


def archived_broad_header():
    samples = []

    with tarfile.open(ARCHIVE, mode="r:") as archive:
        member = archive.getmember(BROAD_MEMBER)
        handle = archive.extractfile(member)

        if handle is None:
            return samples

        for raw_line in handle:
            line = raw_line.decode(
                "utf-8",
                errors="replace",
            )

            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break

    return samples


def focal_list(samples):
    return [
        name for name in FOCAL_NAMES
        if name in samples
    ]


def main():
    missing = [
        path for path in VCFS.values()
        if not path.is_file()
    ]

    if missing:
        print("ERROR: regional VCFs are missing:")
        for path in missing:
            print(path)
        return

    summary_rows = []
    sample_rows = []

    for label, path in VCFS.items():
        result = parse_vcf(path)
        samples = result["samples"]

        summary_rows.append(
            {
                "callset": label,
                "source": "extracted_regional_vcf",
                "basename": path.name,
                "bytes": path.stat().st_size,
                "samples": len(samples),
                "records": result["records"],
                "biallelic_SNPs": (
                    result["biallelic_snps"]
                ),
                "observed_ranges": result["ranges"],
                "focal_present": ",".join(
                    focal_list(samples)
                ),
                "focal_missing": ",".join(
                    name for name in FOCAL_NAMES
                    if name not in samples
                ),
            }
        )

        for order, sample in enumerate(samples, start=1):
            sample_rows.append(
                {
                    "callset": label,
                    "order": order,
                    "sample": sample,
                    "focal_role": (
                        sample
                        if sample in FOCAL_NAMES
                        else ""
                    ),
                }
            )

    broad_samples = archived_broad_header()

    summary_rows.append(
        {
            "callset": "optix_broad_in_tar",
            "source": "header_only_read_from_tar",
            "basename": Path(BROAD_MEMBER).name,
            "bytes": 1997262015,
            "samples": len(broad_samples),
            "records": "",
            "biallelic_SNPs": "",
            "observed_ranges": "",
            "focal_present": ",".join(
                focal_list(broad_samples)
            ),
            "focal_missing": ",".join(
                name for name in FOCAL_NAMES
                if name not in broad_samples
            ),
        }
    )

    for order, sample in enumerate(
        broad_samples,
        start=1,
    ):
        sample_rows.append(
            {
                "callset": "optix_broad_in_tar",
                "order": order,
                "sample": sample,
                "focal_role": (
                    sample
                    if sample in FOCAL_NAMES
                    else ""
                ),
            }
        )

    summary_fields = list(summary_rows[0].keys())
    sample_fields = list(sample_rows[0].keys())

    with OUTPUT_SUMMARY.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=summary_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    with OUTPUT_SAMPLES.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=sample_fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(sample_rows)

    print("ORIGINAL KRISTA PLOT-CALLSET AUDIT")
    print()
    print(
        "callset\tsamples\trecords\tbiallelic_SNPs\t"
        "observed_ranges\tfocal_present"
    )

    for row in summary_rows:
        print(
            f'{row["callset"]}\t'
            f'{row["samples"]}\t'
            f'{row["records"]}\t'
            f'{row["biallelic_SNPs"]}\t'
            f'{row["observed_ranges"]}\t'
            f'{row["focal_present"]}'
        )

    print()
    print(f"Summary: {OUTPUT_SUMMARY}")
    print(f"Sample manifest: {OUTPUT_SAMPLES}")


if __name__ == "__main__":
    main()
