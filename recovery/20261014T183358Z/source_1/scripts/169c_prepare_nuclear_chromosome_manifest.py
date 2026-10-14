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
    "chromosome_calling_manifest.tsv"
)


fai_lengths = {}

with fai.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        fai_lengths[fields[0]] = int(fields[1])


chromosomes = []

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

        # Nuclear chromosomes are numbered 1–21.
        # Chromosome 21 is Z.
        if not molecule.isdigit():
            continue

        chromosome = int(molecule)

        if chromosome < 1 or chromosome > 21:
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
                f"Cannot resolve chromosome {chromosome}; "
                f"candidates={candidates}"
            )

        report_length = int(fields[8])
        fasta_length = fai_lengths[contig]

        if report_length != fasta_length:
            raise RuntimeError(
                f"Length mismatch for chromosome {chromosome}: "
                f"report={report_length}; "
                f"FASTA={fasta_length}; "
                f"contig={contig}"
            )

        chromosomes.append(
            (
                chromosome,
                contig,
                fasta_length,
            )
        )


chromosomes.sort(key=lambda row: row[0])

observed = [
    row[0]
    for row in chromosomes
]

expected = list(range(1, 22))

if observed != expected:
    raise RuntimeError(
        "Expected chromosomes 1 through 21; "
        f"found {observed}"
    )


temporary = Path(
    str(output) + ".new"
)

with temporary.open("w") as handle:
    handle.write(
        "array_id\t"
        "assigned_molecule\t"
        "contig\t"
        "length\t"
        "chromosome_type\n"
    )

    for array_id, row in enumerate(chromosomes):
        chromosome, contig, length = row

        chromosome_type = (
            "Z"
            if chromosome == 21
            else "autosome"
        )

        handle.write(
            f"{array_id}\t"
            f"{chromosome}\t"
            f"{contig}\t"
            f"{length}\t"
            f"{chromosome_type}\n"
        )

temporary.replace(output)

print("NUCLEAR CHROMOSOME MANIFEST COMPLETE")
print(f"chromosomes={len(chromosomes)}")
print("autosomes=20")
print("Z_chromosomes=1")
print("Z_assigned_molecule=21")
print(f"output={output}")
print()
print(output.read_text(), end="")
