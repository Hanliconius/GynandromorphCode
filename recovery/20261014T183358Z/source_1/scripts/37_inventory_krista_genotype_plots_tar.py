#!/usr/bin/env python3

from pathlib import Path, PurePosixPath
from collections import Counter, defaultdict
import csv
import hashlib
import tarfile


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
ARCHIVE = PROJECT / "Genotype_Plots.tar"

CHECKPOINT_DIR = PROJECT / "analysis/checkpoints"

FULL_INVENTORY = (
    CHECKPOINT_DIR
    / "krista_Genotype_Plots_tar_inventory.tsv"
)

RELEVANT_INVENTORY = (
    CHECKPOINT_DIR
    / "krista_Genotype_Plots_tar_relevant_files.tsv"
)

CHECKSUM_FILE = (
    CHECKPOINT_DIR
    / "krista_Genotype_Plots_tar_sha256.txt"
)

INTERESTING_SUFFIXES = (
    ".vcf",
    ".vcf.gz",
    ".bcf",
    ".bcf.gz",
    ".csi",
    ".tbi",
    ".csv",
    ".tsv",
    ".txt",
    ".r",
    ".rmd",
    ".html",
    ".pdf",
    ".svg",
    ".eps",
    ".ai",
    ".xlsx",
    ".xls",
    ".png",
    ".jpg",
    ".jpeg",
)

RELEVANT_TERMS = (
    "optix",
    "ivory",
    "dennis",
    "ray",
    "band",
    "genotype",
    "metadata",
    "vcf",
    "coverage",
    "gynandro",
)


def human_bytes(number):
    units = ["B", "KiB", "MiB", "GiB", "TiB"]
    value = float(number)

    for unit in units:
        if value < 1024 or unit == units[-1]:
            return f"{value:.2f} {unit}"
        value /= 1024

    return f"{number} B"


def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(1024 * 1024)

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


def member_type(member):
    if member.isfile():
        return "file"
    if member.isdir():
        return "directory"
    if member.issym():
        return "symbolic_link"
    if member.islnk():
        return "hard_link"
    if member.ischr():
        return "character_device"
    if member.isblk():
        return "block_device"
    if member.isfifo():
        return "fifo"
    return "other"


def unsafe_reason(member):
    path = PurePosixPath(member.name)
    reasons = []

    if path.is_absolute():
        reasons.append("absolute_path")

    if ".." in path.parts:
        reasons.append("parent_traversal")

    if member.issym() or member.islnk():
        target = PurePosixPath(member.linkname)

        if target.is_absolute():
            reasons.append("absolute_link_target")

        if ".." in target.parts:
            reasons.append("link_parent_traversal")

    if member.ischr() or member.isblk() or member.isfifo():
        reasons.append("special_file")

    return ",".join(reasons)


def is_relevant(name):
    lower = name.lower()

    return (
        lower.endswith(INTERESTING_SUFFIXES)
        or any(term in lower for term in RELEVANT_TERMS)
    )


def main():
    CHECKPOINT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not ARCHIVE.is_file() or ARCHIVE.stat().st_size == 0:
        print("ERROR: archive is missing or empty:")
        print(ARCHIVE)
        print(
            "Run this only after the originating scp command "
            "has completed."
        )
        return

    archive_hash = sha256(ARCHIVE)

    CHECKSUM_FILE.write_text(
        f"{archive_hash}  {ARCHIVE.name}\n"
    )

    rows = []

    try:
        with tarfile.open(ARCHIVE, mode="r:*") as tar:
            members = tar.getmembers()

            for order, member in enumerate(
                members,
                start=1,
            ):
                member_path = PurePosixPath(member.name)

                rows.append(
                    {
                        "order": order,
                        "member": member.name,
                        "basename": member_path.name,
                        "parent": str(member_path.parent),
                        "top_level": (
                            member_path.parts[0]
                            if member_path.parts
                            else ""
                        ),
                        "type": member_type(member),
                        "size_bytes": member.size,
                        "link_target": (
                            member.linkname
                            if member.issym()
                            or member.islnk()
                            else ""
                        ),
                        "unsafe_reason": unsafe_reason(member),
                        "relevant": is_relevant(member.name),
                    }
                )

    except tarfile.TarError as error:
        print("ERROR: archive could not be parsed:")
        print(error)
        return

    fields = [
        "order",
        "member",
        "basename",
        "parent",
        "top_level",
        "type",
        "size_bytes",
        "link_target",
        "unsafe_reason",
        "relevant",
    ]

    with FULL_INVENTORY.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    relevant_rows = [
        row for row in rows
        if row["relevant"]
        and row["type"] == "file"
    ]

    with RELEVANT_INVENTORY.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(relevant_rows)

    file_rows = [
        row for row in rows
        if row["type"] == "file"
    ]

    unsafe_rows = [
        row for row in rows
        if row["unsafe_reason"]
    ]

    type_counts = Counter(
        row["type"] for row in rows
    )

    top_counts = Counter(
        row["top_level"] for row in rows
    )

    basename_paths = defaultdict(list)

    for row in file_rows:
        basename_paths[row["basename"]].append(
            row["member"]
        )

    duplicate_basenames = {
        basename: paths
        for basename, paths in basename_paths.items()
        if len(paths) > 1
    }

    print("KRISTA GENOTYPE_PLOTS TAR AUDIT")
    print()
    print(f"Archive: {ARCHIVE}")
    print(
        f"Archive size: "
        f"{human_bytes(ARCHIVE.stat().st_size)}"
    )
    print(f"SHA256: {archive_hash}")
    print(f"Members: {len(rows)}")
    print(f"Regular files: {len(file_rows)}")
    print(
        "Total uncompressed regular-file size: "
        f"{human_bytes(sum(row['size_bytes'] for row in file_rows))}"
    )
    print(f"Relevant files: {len(relevant_rows)}")
    print(f"Unsafe members: {len(unsafe_rows)}")
    print(
        "Repeated basenames: "
        f"{len(duplicate_basenames)}"
    )

    print()
    print("MEMBER TYPES")
    print("type\tmembers")

    for kind, count in sorted(type_counts.items()):
        print(f"{kind}\t{count}")

    print()
    print("TOP-LEVEL PATHS")
    print("path\tmembers")

    for name, count in sorted(
        top_counts.items(),
        key=lambda item: (-item[1], item[0]),
    ):
        print(f"{name}\t{count}")

    print()
    print("RELEVANT FILES")
    print("size_bytes\tmember")

    for row in relevant_rows:
        print(
            f'{row["size_bytes"]}\t{row["member"]}'
        )

    if duplicate_basenames:
        print()
        print("REPEATED BASENAMES")

        for basename in sorted(duplicate_basenames):
            print()
            print(basename)

            for path in duplicate_basenames[basename]:
                print(f"  {path}")

    if unsafe_rows:
        print()
        print("UNSAFE MEMBERS")

        for row in unsafe_rows:
            print(
                f'{row["unsafe_reason"]}\t'
                f'{row["member"]}\t'
                f'{row["link_target"]}'
            )

        print()
        print(
            "The archive was not extracted because unsafe "
            "members were detected."
        )
    else:
        print()
        print(
            "No unsafe paths, unsafe links, or special files "
            "were detected."
        )
        print("No files were extracted.")

    print()
    print(f"Full inventory: {FULL_INVENTORY}")
    print(f"Relevant inventory: {RELEVANT_INVENTORY}")
    print(f"Checksum: {CHECKSUM_FILE}")


if __name__ == "__main__":
    main()
