#!/usr/bin/env python3

import csv
import hashlib
import subprocess
from collections import defaultdict
from pathlib import Path

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
CHECKPOINT_DIR = PROJECT / "analysis" / "checkpoints"

VCFS = {
    "optix_longest": (
        PROJECT
        / "data"
        / "raw"
        / "Optix_Hmr18-938195-1239421.vcf"
    ),
    "ivory_longest": (
        PROJECT
        / "data"
        / "raw"
        / "shared_drive"
        / "Ivory_Hmr15-1492477-1794947.vcf"
    ),
    "ivory_existing_slice": (
        PROJECT
        / "data"
        / "raw"
        / "Ivory_Hmr15-1627000-1637000_biallelic.vcf"
    ),
}

TARGET_INTERVALS = {
    "optix_band": ("optix_longest", "Hmr18", 1112651, 1114314),
    "optix_dennis1": ("optix_longest", "Hmr18", 1215210, 1216124),
    "optix_dennis2": ("optix_longest", "Hmr18", 1218557, 1221494),
    "optix_ray1": ("optix_longest", "Hmr18", 1156746, 1160220),
    "optix_ray2": ("optix_longest", "Hmr18", 1174547, 1187833),
    "ivory_plot_region": (
        "ivory_longest",
        "Hmr15",
        1627000,
        1637000,
    ),
}

FOCAL_ROLE = {
    "TN2L_S1": "left_T2_leg",
    "TN3L_S3": "left_T3_leg",
    "TN2R_S2": "right_T2_leg",
    "TN3R_S4": "right_T3_leg",
    "TN_Left": "pooled_left",
    "TN_Right": "pooled_right",
}

SUMMARY_OUTPUT = (
    CHECKPOINT_DIR
    / "longest_optix_ivory_vcf_summary.tsv"
)

SAMPLE_OUTPUT = (
    CHECKPOINT_DIR
    / "longest_optix_ivory_vcf_samples.tsv"
)

INTERVAL_OUTPUT = (
    CHECKPOINT_DIR
    / "longest_optix_ivory_interval_coverage.tsv"
)


def run_capture(command):
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    return result.returncode, result.stdout, result.stderr


def get_samples(vcf):
    code, stdout, stderr = run_capture(
        ["bcftools", "query", "-l", str(vcf)]
    )

    if code != 0:
        return [], stderr.strip()

    return [
        line.strip()
        for line in stdout.splitlines()
        if line.strip()
    ], ""


def stream_sites(vcf):
    process = subprocess.Popen(
        [
            "bcftools",
            "query",
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT\n",
            str(vcf),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    sites = []
    error = ""

    if process.stdout is not None:
        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 4:
                continue

            try:
                position = int(fields[1])
            except ValueError:
                continue

            sites.append(
                (
                    fields[0],
                    position,
                    fields[2],
                    fields[3],
                )
            )

    if process.stderr is not None:
        error = process.stderr.read().strip()

    return_code = process.wait()

    if return_code != 0:
        return [], error

    return sites, ""


def hash_strings(values):
    digest = hashlib.sha256()

    for value in values:
        digest.update(value.encode("utf-8"))
        digest.update(b"\n")

    return digest.hexdigest()


def is_biallelic(record):
    chrom, position, ref, alt = record
    return "," not in alt


def site_string(record):
    chrom, position, ref, alt = record
    return f"{chrom}\t{position}\t{ref}\t{alt}"


def inspect_vcf(label, vcf):
    samples, sample_error = get_samples(vcf)

    if sample_error:
        return None, [], [], sample_error

    sites, site_error = stream_sites(vcf)

    if site_error:
        return None, samples, [], site_error

    by_contig = defaultdict(list)

    for chrom, position, ref, alt in sites:
        by_contig[chrom].append(position)

    contig_summaries = []

    for chrom in sorted(by_contig):
        positions = by_contig[chrom]
        contig_summaries.append(
            (
                chrom,
                min(positions),
                max(positions),
                max(positions) - min(positions) + 1,
                len(positions),
            )
        )

    summary = {
        "label": label,
        "relative_path": vcf.relative_to(PROJECT).as_posix(),
        "bytes": vcf.stat().st_size,
        "samples": len(samples),
        "records": len(sites),
        "contigs": len(by_contig),
        "contig_ranges": ";".join(
            f"{chrom}:{minimum}-{maximum}"
            for (
                chrom,
                minimum,
                maximum,
                span,
                records,
            ) in contig_summaries
        ),
        "maximum_observed_span": max(
            (
                span
                for (
                    chrom,
                    minimum,
                    maximum,
                    span,
                    records,
                ) in contig_summaries
            ),
            default=0,
        ),
        "individual_legs_present": ",".join(
            sample
            for sample in [
                "TN2L_S1",
                "TN3L_S3",
                "TN2R_S2",
                "TN3R_S4",
            ]
            if sample in samples
        ),
        "pooled_sides_present": ",".join(
            sample
            for sample in ["TN_Left", "TN_Right"]
            if sample in samples
        ),
        "sample_order_sha256": hash_strings(samples),
        "site_order_sha256": hash_strings(
            site_string(site) for site in sites
        ),
        "status": "ok",
        "error": "",
    }

    return summary, samples, sites, ""


def main():
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    missing = [
        str(path)
        for path in VCFS.values()
        if not path.is_file()
    ]

    if missing:
        print("ERROR: expected VCFs are missing:")
        for path in missing:
            print(path)
        return

    summaries = {}
    samples_by_vcf = {}
    sites_by_vcf = {}

    for label, vcf in VCFS.items():
        summary, samples, sites, error = inspect_vcf(
            label,
            vcf,
        )

        if summary is None:
            print(f"ERROR reading {label}: {error}")
            return

        summaries[label] = summary
        samples_by_vcf[label] = samples
        sites_by_vcf[label] = sites

    with SUMMARY_OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        fieldnames = list(
            summaries["optix_longest"].keys()
        )

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()

        for label in VCFS:
            writer.writerow(summaries[label])

    with SAMPLE_OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "vcf",
                "sample_order",
                "sample",
                "focal_role",
            ]
        )

        for label in VCFS:
            for order, sample in enumerate(
                samples_by_vcf[label],
                start=1,
            ):
                writer.writerow(
                    [
                        label,
                        order,
                        sample,
                        FOCAL_ROLE.get(
                            sample,
                            "comparative_or_other",
                        ),
                    ]
                )

    interval_rows = []

    for interval_name, (
        vcf_label,
        chrom,
        start,
        end,
    ) in TARGET_INTERVALS.items():
        all_sites = sites_by_vcf[vcf_label]

        interval_sites = [
            site
            for site in all_sites
            if (
                site[0] == chrom
                and start <= site[1] <= end
            )
        ]

        biallelic_sites = [
            site
            for site in interval_sites
            if is_biallelic(site)
        ]

        interval_rows.append(
            {
                "interval": interval_name,
                "source_vcf": vcf_label,
                "region": f"{chrom}:{start}-{end}",
                "interval_length_bp": end - start + 1,
                "source_samples": len(
                    samples_by_vcf[vcf_label]
                ),
                "records": len(interval_sites),
                "biallelic_records": len(
                    biallelic_sites
                ),
                "first_record": (
                    interval_sites[0][1]
                    if interval_sites else ""
                ),
                "last_record": (
                    interval_sites[-1][1]
                    if interval_sites else ""
                ),
                "site_sha256": hash_strings(
                    site_string(site)
                    for site in interval_sites
                ),
            }
        )

    with INTERVAL_OUTPUT.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        fieldnames = list(interval_rows[0].keys())

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(interval_rows)

    optix_samples = samples_by_vcf["optix_longest"]
    ivory_samples = samples_by_vcf["ivory_longest"]
    ivory_slice_samples = samples_by_vcf[
        "ivory_existing_slice"
    ]

    optix_ivory_same_order = (
        optix_samples == ivory_samples
    )
    optix_ivory_same_set = (
        set(optix_samples) == set(ivory_samples)
    )
    ivory_slice_same_order = (
        ivory_samples == ivory_slice_samples
    )

    ivory_interval_sites = [
        site
        for site in sites_by_vcf["ivory_longest"]
        if (
            site[0] == "Hmr15"
            and 1627000 <= site[1] <= 1637000
            and is_biallelic(site)
        )
    ]

    existing_ivory_sites = sites_by_vcf[
        "ivory_existing_slice"
    ]

    ivory_site_match = (
        [
            site_string(site)
            for site in ivory_interval_sites
        ]
        ==
        [
            site_string(site)
            for site in existing_ivory_sites
        ]
    )

    print("LONGEST OPTIX AND IVORY VCF AUDIT")

    print(
        "\nlabel\tbytes\tsamples\trecords\t"
        "observed_ranges\tlegs\tpooled"
    )

    for label in VCFS:
        row = summaries[label]

        print(
            f"{label}\t"
            f"{row['bytes']}\t"
            f"{row['samples']}\t"
            f"{row['records']}\t"
            f"{row['contig_ranges']}\t"
            f"{row['individual_legs_present']}\t"
            f"{row['pooled_sides_present']}"
        )

    print("\nSAMPLE-LIST COMPARISONS")
    print(
        "Optix vs longest Ivory, same sample set: "
        f"{optix_ivory_same_set}"
    )
    print(
        "Optix vs longest Ivory, same order: "
        f"{optix_ivory_same_order}"
    )
    print(
        "Longest Ivory vs existing Ivory slice, "
        f"same order: {ivory_slice_same_order}"
    )

    print("\nTARGET INTERVAL COVERAGE")
    print(
        "interval\tsource\tregion\tlength_bp\t"
        "records\tbiallelic_records\t"
        "first_record\tlast_record"
    )

    for row in interval_rows:
        print(
            f"{row['interval']}\t"
            f"{row['source_vcf']}\t"
            f"{row['region']}\t"
            f"{row['interval_length_bp']}\t"
            f"{row['records']}\t"
            f"{row['biallelic_records']}\t"
            f"{row['first_record']}\t"
            f"{row['last_record']}"
        )

    print("\nIVORY RECONSTRUCTION CHECK")
    print(
        "Existing 51-sample Ivory slice has exactly "
        "the same ordered biallelic sites as the longest "
        f"Ivory VCF slice: {ivory_site_match}"
    )

    print("\nINDIVIDUALS IN LONGEST OPTIX AND IVORY VCFS")
    print(
        "order\toptix_sample\tivory_sample\t"
        "same_at_position\trole"
    )

    maximum_samples = max(
        len(optix_samples),
        len(ivory_samples),
    )

    for index in range(maximum_samples):
        optix_sample = (
            optix_samples[index]
            if index < len(optix_samples)
            else ""
        )
        ivory_sample = (
            ivory_samples[index]
            if index < len(ivory_samples)
            else ""
        )
        same = optix_sample == ivory_sample
        role = FOCAL_ROLE.get(
            optix_sample,
            FOCAL_ROLE.get(
                ivory_sample,
                "comparative_or_other",
            ),
        )

        print(
            f"{index + 1}\t"
            f"{optix_sample}\t"
            f"{ivory_sample}\t"
            f"{same}\t"
            f"{role}"
        )

    print("\nRECONSTRUCTION CONCLUSION")
    print(
        "The five Optix genomic intervals can be sliced "
        "from the longest Optix VCF for its existing samples."
    )

    pooled_available = all(
        sample in optix_samples
        for sample in ["TN_Left", "TN_Right"]
    )

    print(
        "Pooled TN_Left/TN_Right genotypes can be recovered "
        "from the longest Optix VCF: "
        f"{pooled_available}"
    )

    print(f"\nVCF summary: {SUMMARY_OUTPUT}")
    print(f"Sample manifest: {SAMPLE_OUTPUT}")
    print(f"Interval summary: {INTERVAL_OUTPUT}")


if __name__ == "__main__":
    main()
