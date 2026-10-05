#!/usr/bin/env python3

import html
import re
from html.parser import HTMLParser
from pathlib import Path

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
SOURCE_DIR = PROJECT / "external" / "krista_analysis"
CHECKPOINT_DIR = PROJECT / "analysis" / "checkpoints"

RMD = SOURCE_DIR / "Optix_and_Ivory_Analysis.Rmd"
HTML_FILES = [
    SOURCE_DIR / "Optix_and_Ivory_Analysis.html",
    SOURCE_DIR / "gynandromorph_genotype_plots.html",
]

CODE_OUTPUT = CHECKPOINT_DIR / "original_plotting_code_blocks.txt"
REFERENCE_OUTPUT = CHECKPOINT_DIR / "original_plotting_file_references.tsv"
RELEVANT_OUTPUT = CHECKPOINT_DIR / "original_plotting_relevant_lines.tsv"

FILE_PATTERN = re.compile(
    r"""(?ix)
    (?:
        [A-Za-z0-9_~@./\\:-]+
    )
    \.
    (?:
        vcf(?:\.gz)?
        |bcf
        |csv
        |tsv
        |txt
        |bed
        |bam
        |rds
        |rdata
        |xlsx?
        |pdf
        |svg
    )
    """
)

RELEVANT_PATTERN = re.compile(
    r"""(?ix)
    \bvcf\b
    |\.vcf
    |read[._]?vcf
    |readVcf
    |vcfR
    |extract\.gt
    |genotype
    |optix
    |ivory
    |band
    |dennis
    |ray
    |TN_Left
    |TN_Right
    |TN2L
    |TN2R
    |TN3L
    |TN3R
    |gynandromorph
    |pool
    |metadata
    |accession
    |sample
    """
)


class CodeBlockParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.pre_depth = 0
        self.code_depth = 0
        self.current = []
        self.blocks = []

    def handle_starttag(self, tag, attrs):
        if tag == "pre":
            self.pre_depth += 1
            if self.pre_depth == 1:
                self.current = []
        elif tag == "code" and self.pre_depth:
            self.code_depth += 1

    def handle_endtag(self, tag):
        if tag == "code" and self.code_depth:
            self.code_depth -= 1
        elif tag == "pre" and self.pre_depth:
            self.pre_depth -= 1
            if self.pre_depth == 0:
                block = "".join(self.current).strip()
                if block:
                    self.blocks.append(block)
                self.current = []

    def handle_data(self, data):
        if self.pre_depth:
            self.current.append(data)


def clean_line(line):
    line = html.unescape(line)
    line = line.replace("\t", " ")
    line = re.sub(r"\s+", " ", line)
    return line.strip()


def extract_rmd_code(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    blocks = []

    fenced = False
    current = []

    for line in text.splitlines():
        if line.startswith("```"):
            if fenced:
                block = "\n".join(current).strip()
                if block:
                    blocks.append(block)
                current = []
                fenced = False
            else:
                fenced = True
                current = []
        elif fenced:
            current.append(line)

    return text, blocks


def extract_html_code(path):
    parser = CodeBlockParser()

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            parser.feed(chunk)

    parser.close()
    return parser.blocks


def main():
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    missing = [
        str(path)
        for path in [RMD, *HTML_FILES]
        if not path.is_file()
    ]

    if missing:
        print("ERROR: expected analysis source files are missing:")
        for path in missing:
            print(path)
        return 1

    rmd_text, rmd_blocks = extract_rmd_code(RMD)

    source_blocks = [(RMD.name, block) for block in rmd_blocks]

    for path in HTML_FILES:
        for block in extract_html_code(path):
            source_blocks.append((path.name, block))

    with CODE_OUTPUT.open("w", encoding="utf-8") as handle:
        for block_number, (source, block) in enumerate(
            source_blocks, start=1
        ):
            handle.write(
                f"===== SOURCE: {source}; BLOCK: {block_number} =====\n"
            )
            handle.write(block.rstrip())
            handle.write("\n\n")

    references = set()
    relevant_rows = []

    # Scan the original Rmd prose and code.
    for line_number, line in enumerate(rmd_text.splitlines(), start=1):
        cleaned = clean_line(line)

        if not cleaned:
            continue

        for match in FILE_PATTERN.findall(cleaned):
            references.add((RMD.name, match))

        if RELEVANT_PATTERN.search(cleaned):
            relevant_rows.append(
                (RMD.name, f"line_{line_number}", cleaned)
            )

    # Scan extracted HTML code blocks, not the embedded image data.
    for block_number, (source, block) in enumerate(source_blocks, start=1):
        if source == RMD.name:
            continue

        for line_number, line in enumerate(block.splitlines(), start=1):
            cleaned = clean_line(line)

            if not cleaned:
                continue

            for match in FILE_PATTERN.findall(cleaned):
                references.add((source, match))

            if RELEVANT_PATTERN.search(cleaned):
                relevant_rows.append(
                    (
                        source,
                        f"block_{block_number}_line_{line_number}",
                        cleaned,
                    )
                )

    with REFERENCE_OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write("source\treferenced_file\n")
        for source, reference in sorted(references):
            handle.write(f"{source}\t{reference}\n")

    with RELEVANT_OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write("source\tlocation\ttext\n")
        for source, location, text in relevant_rows:
            handle.write(f"{source}\t{location}\t{text}\n")

    print("ORIGINAL PLOTTING SOURCE AUDIT")
    print(f"Rmd code blocks: {len(rmd_blocks)}")

    for path in HTML_FILES:
        count = sum(
            source == path.name for source, block in source_blocks
        )
        print(f"{path.name} code blocks: {count}")

    print(f"Distinct file references: {len(references)}")
    print(f"Relevant code/provenance lines: {len(relevant_rows)}")

    print("\nREFERENCED INPUT/OUTPUT FILES")
    if references:
        print("source\treferenced_file")
        for source, reference in sorted(references):
            print(f"{source}\t{reference}")
    else:
        print("No file references detected.")

    print("\nRELEVANT LINES FROM ORIGINAL RMD")
    rmd_rows = [
        row for row in relevant_rows
        if row[0] == RMD.name
    ]

    for source, location, text in rmd_rows:
        print(f"{location}\t{text}")

    print("\nFull extracted code:")
    print(CODE_OUTPUT)
    print("File-reference table:")
    print(REFERENCE_OUTPUT)
    print("Relevant-line table:")
    print(RELEVANT_OUTPUT)

    return 0


if __name__ == "__main__":
    main()
