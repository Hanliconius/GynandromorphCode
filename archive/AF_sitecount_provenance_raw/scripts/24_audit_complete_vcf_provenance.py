#!/usr/bin/env python3

import csv
import hashlib
import subprocess
from collections import defaultdict
from pathlib import Path

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
RAW_DIR = PROJECT / "data" / "raw"
CHECKPOINT_DIR = PROJECT / "analysis" / "checkpoints"

OUTPUT = CHECKPOINT_DIR / "complete_vcf_provenance_manifest.tsv"

FOCAL_SAMPLES = [
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
]

EXPECTED_PLOT_VCFS = [
    "optix_band-1112651-1114314.vcf",
    "dennis1-1215210-1216124.vcf",
    "dennis2-1218557-1221494.vcf",
    "ray1-1156746-1160220.vcf",
    "ray2-1174547-1187833.vcf",
]


def infer_source(path):
    relative = path.relative_to(PROJECT).as_posix()

    if relative.startswith("data/raw/region_panel/"):
        return "vcf_files.zip"

    if relative.startswith("data/raw/ivory_variants/"):
        return "Ivory_VCFs.zip"

    if relative.startswith("data/raw/shared_drive/"):
        return "shared_Drive_individual_file"

    if relative.startswith("data/raw/"):
        return "earlier_shared_Drive_download"

    return "unknown"


def file_sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)

            if not chunk:
                break

            digest.update(chunk)

    return digest.hexdigest()


def run_capture(command):
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )

    return result.returncode, result.stdout, result.stderr


def stream_command_hash(command, count_lines=False):
    digest = hashlib.sha256()
    line_count = 0
    first_chrom = ""
    last_chrom = ""
    minimum_position = None
    maximum_position = None

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=False,
    )

    if process.stdout is not None:
        for raw_line in process.stdout:
            digest.update(raw_line)

            if count_lines:
                line_count += 1

                fields = raw_line.rstrip(b"\n").split(b"\t")

                if len(fields) >= 2:
                    chrom = fields[0].decode(
                        "utf-8", errors="replace"
                    )

                    try:
                        position = int(fields[1])
                    except ValueError:
                        position = None

                    if not first_chrom:
                        first_chrom = chrom

                    last_chrom = chrom

                    if position is not None:
                        if minimum_position is None:
                            minimum_position = position
                        else:
                            minimum_position = min(
                                minimum_position, position
                            )

                        if maximum_position is None:
                            maximum_position = position
                        else:
                            maximum_position = max(
                                maximum_position, position
                            )

    return_code = process.wait()

    return {
        "return_code": return_code,
        "sha256": digest.hexdigest(),
        "line_count": line_count,
        "first_chrom": first_chrom,
        "last_chrom": last_chrom,
        "minimum_position": minimum_position,
        "maximum_position": maximum_position,
    }


def inspect_vcf(path):
    code, sample_text, sample_error = run_capture(
        ["bcftools", "query", "-l", str(path)]
    )

    if code != 0:
        return {
            "status": "bcftools_sample_error",
            "error": sample_error.strip().replace("\t", " ")[:300],
        }

    samples = [
        line.strip()
        for line in sample_text.splitlines()
        if line.strip()
    ]

    sample_digest = hashlib.sha256()
    sample_digest.update(
        ("\n".join(samples) + "\n").encode("utf-8")
    )

    site_result = stream_command_hash(
        [
            "bcftools",
            "query",
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT\n",
            str(path),
        ],
        count_lines=True,
    )

    if site_result["return_code"] != 0:
        return {
            "status": "bcftools_site_error",
            "error": "Could not stream site records",
        }

    genotype_result = stream_command_hash(
        [
            "bcftools",
            "query",
            "-f",
            "%CHROM\t%POS[\t%GT]\n",
            str(path),
        ],
        count_lines=False,
    )

    if genotype_result["return_code"] != 0:
        return {
            "status": "bcftools_genotype_error",
            "error": "Could not stream genotype records",
        }

    sample_set = set(samples)

    focal_present = [
        sample for sample in FOCAL_SAMPLES
        if sample in sample_set
    ]

    return {
        "status": "ok",
        "error": "",
        "samples": len(samples),
        "variants": site_result["line_count"],
        "sample_sha256": sample_digest.hexdigest(),
        "site_sha256": site_result["sha256"],
        "genotype_sha256": genotype_result["sha256"],
        "first_chrom": site_result["first_chrom"],
        "last_chrom": site_result["last_chrom"],
        "minimum_position": site_result["minimum_position"],
        "maximum_position": site_result["maximum_position"],
        "focal_present": ",".join(focal_present),
        "individual_legs_present": sum(
            sample in sample_set
            for sample in FOCAL_SAMPLES[:4]
        ),
        "pooled_sides_present": sum(
            sample in sample_set
            for sample in FOCAL_SAMPLES[4:]
        ),
    }


def main():
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    vcfs = sorted(
        path
        for path in RAW_DIR.rglob("*")
        if path.is_file()
        and not path.name.endswith(".part")
        and (
            path.name.endswith(".vcf")
            or path.name.endswith(".vcf.gz")
            or path.name.endswith(".bcf")
        )
    )

    rows = []

    for path in vcfs:
        result = inspect_vcf(path)

        row = {
            "source": infer_source(path),
            "basename": path.name,
            "relative_path": path.relative_to(PROJECT).as_posix(),
            "bytes": path.stat().st_size,
            "file_sha256": file_sha256(path),
            "status": result.get("status", ""),
            "samples": result.get("samples", ""),
            "variants": result.get("variants", ""),
            "sample_sha256": result.get("sample_sha256", ""),
            "site_sha256": result.get("site_sha256", ""),
            "genotype_sha256": result.get("genotype_sha256", ""),
            "first_chrom": result.get("first_chrom", ""),
            "last_chrom": result.get("last_chrom", ""),
            "minimum_position": result.get("minimum_position", ""),
            "maximum_position": result.get("maximum_position", ""),
            "individual_legs_present": result.get(
                "individual_legs_present", ""
            ),
            "pooled_sides_present": result.get(
                "pooled_sides_present", ""
            ),
            "focal_present": result.get("focal_present", ""),
            "error": result.get("error", ""),
        }

        rows.append(row)

    fieldnames = [
        "source",
        "basename",
        "relative_path",
        "bytes",
        "status",
        "samples",
        "variants",
        "first_chrom",
        "last_chrom",
        "minimum_position",
        "maximum_position",
        "individual_legs_present",
        "pooled_sides_present",
        "focal_present",
        "file_sha256",
        "sample_sha256",
        "site_sha256",
        "genotype_sha256",
        "error",
    ]

    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    print("COMPLETE VCF PROVENANCE AUDIT")
    print(f"VCF/BCF files found: {len(rows)}")

    print("\nVCF CALLSETS")
    print(
        "source\tbasename\tsamples\tvariants\t"
        "legs\tpooled\tfile_hash\tsite_hash\tgenotype_hash"
    )

    for row in rows:
        print(
            f"{row['source']}\t"
            f"{row['basename']}\t"
            f"{row['samples']}\t"
            f"{row['variants']}\t"
            f"{row['individual_legs_present']}\t"
            f"{row['pooled_sides_present']}\t"
            f"{row['file_sha256'][:12]}\t"
            f"{row['site_sha256'][:12]}\t"
            f"{row['genotype_sha256'][:12]}"
        )

    basename_groups = defaultdict(list)
    file_hash_groups = defaultdict(list)
    site_sample_groups = defaultdict(list)
    genotype_groups = defaultdict(list)

    for row in rows:
        basename_groups[row["basename"]].append(row)

        if row["file_sha256"]:
            file_hash_groups[row["file_sha256"]].append(row)

        if row["sample_sha256"] and row["site_sha256"]:
            key = (
                row["sample_sha256"],
                row["site_sha256"],
            )
            site_sample_groups[key].append(row)

        if (
            row["sample_sha256"]
            and row["site_sha256"]
            and row["genotype_sha256"]
        ):
            key = (
                row["sample_sha256"],
                row["site_sha256"],
                row["genotype_sha256"],
            )
            genotype_groups[key].append(row)

    print("\nREPEATED BASENAMES")
    repeated_basenames = {
        key: group
        for key, group in basename_groups.items()
        if len(group) > 1
    }

    if repeated_basenames:
        for basename, group in sorted(repeated_basenames.items()):
            print(f"\n{basename}")
            for row in group:
                print(
                    f"  {row['source']}\t"
                    f"{row['relative_path']}\t"
                    f"{row['bytes']} bytes\t"
                    f"{row['samples']} samples\t"
                    f"{row['variants']} variants\t"
                    f"file={row['file_sha256'][:12]}\t"
                    f"sites={row['site_sha256'][:12]}\t"
                    f"GT={row['genotype_sha256'][:12]}"
                )
    else:
        print("None")

    print("\nEXACT BYTE-DUPLICATE FILES")
    exact_groups = [
        group
        for group in file_hash_groups.values()
        if len(group) > 1
    ]

    if exact_groups:
        for group in exact_groups:
            print(
                " | ".join(
                    row["relative_path"] for row in group
                )
            )
    else:
        print("None")

    print("\nSAME SAMPLE AND SITE SETS")
    logical_groups = [
        group
        for group in site_sample_groups.values()
        if len(group) > 1
    ]

    if logical_groups:
        for group in logical_groups:
            print(
                " | ".join(
                    row["relative_path"] for row in group
                )
            )
    else:
        print("None")

    print("\nIDENTICAL SAMPLE, SITE, AND GT STREAMS")
    gt_groups = [
        group
        for group in genotype_groups.values()
        if len(group) > 1
    ]

    if gt_groups:
        for group in gt_groups:
            print(
                " | ".join(
                    row["relative_path"] for row in group
                )
            )
    else:
        print("None")

    available_basenames = set(basename_groups)

    print("\nPLOT-SPECIFIC OPTIX INPUTS")
    for expected in EXPECTED_PLOT_VCFS:
        status = (
            "present"
            if expected in available_basenames
            else "missing_as_standalone_file"
        )
        print(f"{expected}\t{status}")

    print(f"\nFull manifest: {OUTPUT}")


if __name__ == "__main__":
    main()
