#!/usr/bin/env python3

from pathlib import Path


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/"
    "martinlab/joe/GynandroCalling"
)

META = (
    PROJECT
    / "analysis/genomewide_mixing_input/metadata"
)

REF_DIR = (
    PROJECT
    / "analysis/genomewide_mixing_input/reference/"
    "GCF_056956725.1"
)

REPORT = (
    REF_DIR
    / "GCF_056956725.1_ASM5695672v1_assembly_report.txt"
)

FAI = (
    REF_DIR
    / "GCF_056956725.1_ASM5695672v1_genomic.fna.fai"
)

OUTPUT = META / "autosome_calling_manifest.tsv"
AUDIT = META / "autosome_identifier_resolution.tsv"


def main():
    if not REPORT.is_file():
        raise RuntimeError(
            f"Missing assembly report: {REPORT}"
        )

    if not FAI.is_file():
        raise RuntimeError(
            f"Missing FASTA index: {FAI}"
        )

    fai_lengths = {}

    with FAI.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            fai_lengths[fields[0]] = int(fields[1])

    autosomes = []
    audit_rows = []
    z_count = 0

    with REPORT.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            sequence_name = fields[0]
            sequence_role = fields[1]
            molecule = fields[2]
            genbank = fields[4]
            refseq = fields[6]
            report_length = int(fields[8])

            if sequence_role != "assembled-molecule":
                continue

            if molecule.upper() == "Z":
                z_count += 1
                continue

            if not molecule.isdigit():
                continue

            candidates = [
                refseq,
                genbank,
                sequence_name,
            ]

            if len(fields) > 9:
                candidates.append(fields[9])

            candidates = [
                value
                for value in candidates
                if value and value.lower() != "na"
            ]

            matches = [
                value
                for value in candidates
                if value in fai_lengths
            ]

            if not matches:
                raise RuntimeError(
                    "No FASTA identifier matched autosome "
                    f"{molecule}. Candidates: {candidates}"
                )

            contig = matches[0]
            fasta_length = fai_lengths[contig]

            if fasta_length != report_length:
                raise RuntimeError(
                    f"Length mismatch for autosome {molecule}: "
                    f"report={report_length}; "
                    f"FASTA={fasta_length}; "
                    f"contig={contig}"
                )

            autosomes.append(
                (
                    int(molecule),
                    contig,
                    fasta_length,
                )
            )

            audit_rows.append(
                (
                    molecule,
                    sequence_name,
                    genbank,
                    refseq,
                    contig,
                    fasta_length,
                )
            )

    autosomes.sort(key=lambda row: row[0])

    observed = [
        row[0]
        for row in autosomes
    ]

    expected = list(range(1, 21))

    if observed != expected:
        raise RuntimeError(
            "Expected autosomes 1 through 20; "
            f"found {observed}"
        )

    if z_count != 1:
        raise RuntimeError(
            "Expected exactly one assembled Z chromosome; "
            f"found {z_count}"
        )

    META.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_new = Path(
        str(OUTPUT) + ".new"
    )

    with output_new.open("w") as handle:
        handle.write(
            "array_id\t"
            "assigned_molecule\t"
            "contig\t"
            "length\n"
        )

        for array_id, row in enumerate(autosomes):
            molecule, contig, length = row

            handle.write(
                f"{array_id}\t"
                f"{molecule}\t"
                f"{contig}\t"
                f"{length}\n"
            )

    audit_new = Path(
        str(AUDIT) + ".new"
    )

    with audit_new.open("w") as handle:
        handle.write(
            "assigned_molecule\t"
            "sequence_name\t"
            "genbank_accession\t"
            "refseq_accession\t"
            "selected_fasta_contig\t"
            "length\n"
        )

        for row in audit_rows:
            handle.write(
                "\t".join(map(str, row)) + "\n"
            )

    output_new.replace(OUTPUT)
    audit_new.replace(AUDIT)

    print("AUTOSOME CALLING MANIFEST COMPLETE")
    print(f"autosomes={len(autosomes)}")
    print(f"excluded_Z={z_count}")
    print(f"output={OUTPUT}")
    print(f"audit={AUDIT}")

    print("\nAUTOSOME MANIFEST")
    print(OUTPUT.read_text(), end="")


if __name__ == "__main__":
    main()
