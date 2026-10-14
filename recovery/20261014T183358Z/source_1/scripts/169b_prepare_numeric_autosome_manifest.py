#!/usr/bin/env python3

from pathlib import Path


project = Path(
    "/gpfs/automountdir/gpfs/scratch/"
    "martinlab/joe/GynandroCalling"
)

reference_dir = (
    project
    / "analysis/genomewide_mixing_input/reference/"
    "GCF_056956725.1"
)

report = (
    reference_dir
    / "GCF_056956725.1_ASM5695672v1_assembly_report.txt"
)

fai = (
    reference_dir
    / "GCF_056956725.1_ASM5695672v1_genomic.fna.fai"
)

output = (
    project
    / "analysis/genomewide_mixing_input/metadata/"
    "autosome_calling_manifest.tsv"
)


fai_lengths = {}

with fai.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        fai_lengths[fields[0]] = int(fields[1])


autosomes = []
z_rows = 0

with report.open() as handle:
    for line in handle:
        if line.startswith("#") or not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 9:
            continue

        role = fields[1]
        molecule = fields[2]

        if role != "assembled-molecule":
            continue

        if molecule.upper() == "Z":
            z_rows += 1
            continue

        # Exclude non-autosomes before resolving accessions.
        if not molecule.isdigit():
            continue

        number = int(molecule)

        if number < 1 or number > 20:
            continue

        candidates = [
            fields[6],
            fields[4],
            fields[0],
        ]

        if len(fields) > 9:
            candidates.append(fields[9])

        contig = next(
            (
                candidate
                for candidate in candidates
                if (
                    candidate
                    and candidate.lower() != "na"
                    and candidate in fai_lengths
                )
            ),
            None,
        )

        if contig is None:
            raise RuntimeError(
                f"Cannot resolve autosome {number}; "
                f"candidates={candidates}"
            )

        report_length = int(fields[8])
        fasta_length = fai_lengths[contig]

        if report_length != fasta_length:
            raise RuntimeError(
                f"Length mismatch for autosome {number}: "
                f"report={report_length}; "
                f"FASTA={fasta_length}"
            )

        autosomes.append(
            (
                number,
                contig,
                fasta_length,
            )
        )


autosomes.sort(key=lambda row: row[0])

numbers = [
    row[0]
    for row in autosomes
]

if numbers != list(range(1, 21)):
    raise RuntimeError(
        "Expected one copy of autosomes 1–20; "
        f"found {numbers}"
    )

if z_rows != 1:
    raise RuntimeError(
        "Expected exactly one excluded Z chromosome; "
        f"found {z_rows}"
    )


temporary = Path(
    str(output) + ".new"
)

with temporary.open("w") as handle:
    handle.write(
        "array_id\t"
        "assigned_molecule\t"
        "contig\t"
        "length\n"
    )

    for array_id, row in enumerate(autosomes):
        number, contig, length = row

        handle.write(
            f"{array_id}\t"
            f"{number}\t"
            f"{contig}\t"
            f"{length}\n"
        )

temporary.replace(output)

print("NUMERIC AUTOSOME MANIFEST COMPLETE")
print(f"autosomes={len(autosomes)}")
print(f"excluded_Z={z_rows}")
print(f"output={output}")
print()
print(output.read_text(), end="")
