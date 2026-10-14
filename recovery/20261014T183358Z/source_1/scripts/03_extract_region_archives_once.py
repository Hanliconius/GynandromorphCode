#!/usr/bin/env python3

import shutil
import zlib
import zipfile
from pathlib import Path, PurePosixPath


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

ARCHIVE_TARGETS = {
    PROJECT / "data/incoming/vcf_files.zip":
        PROJECT / "data/raw/region_panel",
    PROJECT / "data/incoming/Ivory_VCFs.zip":
        PROJECT / "data/raw/ivory_variants",
}


def crc32_file(path):
    checksum = 0

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            checksum = zlib.crc32(block, checksum)

    return checksum & 0xFFFFFFFF


def safe_member(name):
    path = PurePosixPath(name)
    return (
        not path.is_absolute()
        and ".." not in path.parts
        and not name.startswith("\\")
    )


def main():
    failed = False

    for archive_path, target_dir in ARCHIVE_TARGETS.items():
        if not archive_path.is_file():
            print(f"ERROR\tmissing archive\t{archive_path}")
            failed = True
            continue

        target_dir.mkdir(parents=True, exist_ok=True)

        with zipfile.ZipFile(archive_path) as archive:
            for member in archive.infolist():
                if member.is_dir():
                    continue

                if not safe_member(member.filename):
                    print(
                        f"ERROR\tunsafe archive member\t"
                        f"{archive_path.name}\t{member.filename}"
                    )
                    failed = True
                    continue

                relative = Path(*PurePosixPath(member.filename).parts)
                destination = target_dir / relative
                temporary = destination.with_name(destination.name + ".part")

                destination.parent.mkdir(parents=True, exist_ok=True)

                if destination.is_file():
                    size_matches = destination.stat().st_size == member.file_size
                    crc_matches = (
                        crc32_file(destination) == member.CRC
                        if size_matches
                        else False
                    )

                    if size_matches and crc_matches:
                        print(f"REUSE\t{destination}")
                    else:
                        print(
                            f"ERROR\texisting file differs from archive\t"
                            f"{destination}"
                        )
                        failed = True

                    continue

                if temporary.exists():
                    temporary.unlink()

                with archive.open(member) as source, temporary.open("wb") as output:
                    shutil.copyfileobj(source, output)

                size_matches = temporary.stat().st_size == member.file_size
                crc_matches = (
                    crc32_file(temporary) == member.CRC
                    if size_matches
                    else False
                )

                if size_matches and crc_matches:
                    temporary.replace(destination)
                    print(f"EXTRACT\t{destination}")
                else:
                    print(f"ERROR\tfailed integrity check\t{destination}")
                    temporary.unlink(missing_ok=True)
                    failed = True

    if failed:
        print("\nExtraction encountered one or more errors.")
        return 1

    print("\nExtraction completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
