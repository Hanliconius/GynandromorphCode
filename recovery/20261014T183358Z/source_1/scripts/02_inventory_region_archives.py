#!/usr/bin/env python3

import csv
import zipfile
from pathlib import Path, PurePosixPath


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
INCOMING = PROJECT / "data" / "incoming"
CHECKPOINTS = PROJECT / "analysis" / "checkpoints"

ARCHIVES = [
    INCOMING / "vcf_files.zip",
    INCOMING / "Ivory_VCFs.zip",
]

VCF_SUFFIXES = (
    ".vcf",
    ".vcf.gz",
    ".bcf",
    ".tbi",
    ".csi",
)


def unsafe_path(name):
    path = PurePosixPath(name)
    return path.is_absolute() or ".." in path.parts or name.startswith("\\")


def member_type(name):
    lowered = name.lower()

    for suffix in VCF_SUFFIXES:
        if lowered.endswith(suffix):
            return suffix.removeprefix(".")

    suffix = PurePosixPath(name).suffix.lower()
    return suffix.removeprefix(".") if suffix else "no_extension"


def main():
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)

    missing = [archive for archive in ARCHIVES if not archive.is_file()]

    if missing:
        print("MISSING ARCHIVES")
        for archive in missing:
            print(archive)
        return 1

    output = CHECKPOINTS / "region_archive_members.tsv"
    summaries = []
    vcf_members = []

    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "archive",
                "member",
                "type",
                "uncompressed_bytes",
                "compressed_bytes",
                "unsafe_path",
            ]
        )

        for archive_path in ARCHIVES:
            with zipfile.ZipFile(archive_path) as archive:
                members = [
                    member
                    for member in archive.infolist()
                    if not member.is_dir()
                ]

            unsafe_count = 0
            total_uncompressed = 0
            sequence_count = 0

            for member in members:
                unsafe = unsafe_path(member.filename)
                kind = member_type(member.filename)

                unsafe_count += int(unsafe)
                total_uncompressed += member.file_size

                if member.filename.lower().endswith(VCF_SUFFIXES):
                    sequence_count += 1
                    vcf_members.append(
                        (
                            archive_path.name,
                            member.filename,
                            member.file_size,
                            kind,
                        )
                    )

                writer.writerow(
                    [
                        archive_path.name,
                        member.filename,
                        kind,
                        member.file_size,
                        member.compress_size,
                        int(unsafe),
                    ]
                )

            summaries.append(
                (
                    archive_path.name,
                    archive_path.stat().st_size,
                    len(members),
                    total_uncompressed,
                    sequence_count,
                    unsafe_count,
                )
            )

    print("REGION ARCHIVE INVENTORY")
    print(
        "archive\tarchive_bytes\tmembers\t"
        "uncompressed_bytes\tvcf_index_members\tunsafe_paths"
    )

    for row in summaries:
        print("\t".join(str(value) for value in row))

    print("\nVCF/BCF AND INDEX MEMBERS")
    print("archive\tmember\tuncompressed_bytes\ttype")

    for row in vcf_members:
        print("\t".join(str(value) for value in row))

    print(f"\nFull inventory: {output}")
    print("No files were extracted and no VCF records were read.")

    if any(row[-1] > 0 for row in summaries):
        print("ERROR: unsafe archive paths detected; do not extract.")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
