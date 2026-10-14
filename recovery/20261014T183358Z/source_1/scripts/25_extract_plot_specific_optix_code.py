#!/usr/bin/env python3

import hashlib
import re
from pathlib import Path

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
INPUT = (
    PROJECT
    / "analysis"
    / "checkpoints"
    / "original_plotting_code_blocks.txt"
)
OUTPUT = (
    PROJECT
    / "analysis"
    / "checkpoints"
    / "plot_specific_optix_code.txt"
)

TARGETS = [
    "Optix_Hmr18-938195-1239421.vcf",
    "optix_band-1112651-1114314.vcf",
    "dennis1-1215210-1216124.vcf",
    "dennis2-1218557-1221494.vcf",
    "ray1-1156746-1160220.vcf",
    "ray2-1174547-1187833.vcf",
    "Optix_Band_metadata.csv",
    "Optix_Dennis_metadata.csv",
    "Optix_Ray_metadata.csv",
    "Optix_metadata.csv",
]

HEADER_PATTERN = re.compile(
    r"^===== SOURCE: (.*?); BLOCK: ([0-9]+) =====$",
    re.MULTILINE,
)


def normalize_block(text):
    lines = [
        line.rstrip()
        for line in text.strip().splitlines()
    ]
    return "\n".join(lines).strip()


def main():
    if not INPUT.is_file():
        print(f"ERROR: missing input: {INPUT}")
        return

    text = INPUT.read_text(
        encoding="utf-8",
        errors="replace",
    )

    matches = list(HEADER_PATTERN.finditer(text))
    blocks = []

    for index, match in enumerate(matches):
        start = match.end()

        if index + 1 < len(matches):
            end = matches[index + 1].start()
        else:
            end = len(text)

        source = match.group(1)
        block_number = match.group(2)
        block = normalize_block(text[start:end])

        if any(target in block for target in TARGETS):
            blocks.append(
                {
                    "source": source,
                    "block_number": block_number,
                    "text": block,
                }
            )

    unique_blocks = []
    seen_hashes = set()

    for block in blocks:
        digest = hashlib.sha256(
            block["text"].encode("utf-8")
        ).hexdigest()

        if digest not in seen_hashes:
            seen_hashes.add(digest)
            block["sha256"] = digest
            unique_blocks.append(block)

    with OUTPUT.open("w", encoding="utf-8") as handle:
        for output_number, block in enumerate(
            unique_blocks,
            start=1,
        ):
            handle.write(
                "===== UNIQUE BLOCK "
                f"{output_number}; "
                f"SOURCE: {block['source']}; "
                f"ORIGINAL BLOCK: {block['block_number']}; "
                f"SHA256: {block['sha256']} =====\n"
            )
            handle.write(block["text"])
            handle.write("\n\n")

    print("PLOT-SPECIFIC OPTIX CODE EXTRACTION")
    print(f"Matching source blocks: {len(blocks)}")
    print(f"Unique matching blocks: {len(unique_blocks)}")
    print(f"Output: {OUTPUT}")

    print("\nEXTRACTED CODE\n")

    if unique_blocks:
        print(OUTPUT.read_text(
            encoding="utf-8",
            errors="replace",
        ))
    else:
        print("No matching code blocks were found.")


if __name__ == "__main__":
    main()
